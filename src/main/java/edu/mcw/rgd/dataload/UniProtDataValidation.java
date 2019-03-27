package edu.mcw.rgd.dataload;

import java.util.*;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.GenomicElement;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 */
public class UniProtDataValidation {
    UniProtDAO dbDao;
    java.util.Map<String,Integer> activeXdbIdMap;
    int speciesTypeKey;

    private Map<Integer, UniProtRecord> records = new HashMap<>();
    private Map<Integer, List<ProteinDomain>> domainsMap = new HashMap<>();

    /**
     * assign incoming data to matching rgd ids;
     * note: a set of incoming data could be assigned to multiple active rgd ids
     * @param incomingRecords
     */
    public void mergeIncomingRecords(List<UniProtRatRecord> incomingRecords) throws Exception {

        int insertedDomains = 0;
        int upToDateDomains = 0;

        for( UniProtRatRecord incomingRec: incomingRecords ) {

            for( Integer activeRgdId: incomingRec.matchingRgdIds ) {

                UniProtRecord rec = records.get(activeRgdId);
                if( rec==null ) {
                    rec = new UniProtRecord(activeXdbIdMap);
                    records.put(activeRgdId, rec);
                }

                rec.init(activeRgdId, incomingRec);
            }

            // merge protein domains
            for (ProteinDomain pd : incomingRec.domains) {

                if (pd.geInRgd == null) {
                    GenomicElement ge = dbDao.insertDomainName(pd.getDomainName());
                    // if ge has been inserted, its objectStatus property will be null
                    if (ge.getObjectStatus() == null) {
                        insertedDomains++;
                    } else {
                        upToDateDomains++;
                    }
                    pd.geInRgd = ge;
                }

                List<ProteinDomain> list = domainsMap.get(pd.geInRgd.getRgdId());
                if (list == null) {
                    list = new ArrayList<>();
                    domainsMap.put(pd.geInRgd.getRgdId(), list);
                }
                list.add(pd);
            }
        }

        System.out.println("DOMAINS INSERTED: "+insertedDomains);
        System.out.println("DOMAINS UP-TO-DATE: "+upToDateDomains);
    }

    public void load() throws Exception {

        int primaryMapKey = MapManager.getInstance().getReferenceAssembly(speciesTypeKey).getKey();

        System.out.println(" loading uniprot xdb ids ...");
        for( UniProtRecord rec: records.values() ) {
            load(rec);
        }

        System.out.println(" loading protein domains ...");
        for (Map.Entry<Integer, List<ProteinDomain>> entry : domainsMap.entrySet()) {
            System.out.println("processing PD " + entry.getKey());

            List<MapData> domainLoci = new ArrayList<>();
            for (ProteinDomain pd : entry.getValue()) {
                if( pd.loci!=null ) {
                    for (MapData md : pd.loci) {
                        addDomainLoci(domainLoci, md);
                    }
                }
            }
            updateDomainLociInDb(entry.getKey(), primaryMapKey, "UniProtKB", domainLoci);
        }

        System.out.println(" loading OK!");
    }

    void addDomainLoci(List<MapData> loci, MapData md) {

        // see for duplicate loci, with same chr, strand, start and stop
        for( MapData mdLoci: loci ) {
            if( mdLoci.getChromosome().equals(md.getChromosome())
                    && mdLoci.getStrand().equals(md.getStrand())
                    && mdLoci.getStartPos().equals(md.getStartPos())
                    && mdLoci.getStopPos().equals(md.getStopPos()) ) {

                System.out.println("-- protein-domain: merging duplicate loci");
                if( Utils.isStringEmpty(mdLoci.getNotes()) ) {
                    mdLoci.setNotes(md.getNotes());
                } else if( !Utils.isStringEmpty(md.getNotes()) ) {
                    mdLoci.setNotes( mdLoci.getNotes()+"; "+md.getNotes() );
                }
            }
        }
    }

    void updateDomainLociInDb( int domainRgdId, int mapKey, String srcPipeline, List<MapData> loci ) throws Exception {
        MapDAO mdao = new MapDAO();
        List<MapData> mdsInRgd = mdao.getMapData(domainRgdId, mapKey, srcPipeline);

        // if incoming locus has a match in RGD, remove them from the lists
        List<MapData> mdsUpToDate = new ArrayList<>(mdsInRgd.size());

        Iterator<MapData> it = loci.iterator();
        while( it.hasNext() ) {
            MapData mdIncoming = it.next();

            // find a match in RGD
            Iterator<MapData> itInRgd = mdsInRgd.iterator();
            while( itInRgd.hasNext() ) {
                MapData mdInRgd = itInRgd.next();
                if( mdInRgd.equalsByGenomicCoords(mdIncoming) ) {
                    mdsUpToDate.add(mdInRgd);
                    itInRgd.remove();
                    it.remove();
                    break;
                }
            }
        }

        if( !mdsInRgd.isEmpty() ) {
            System.out.println("LOCI to be removed from RGD");
        }
        if( !loci.isEmpty() ) {
            dbDao.insertMapData(loci);
        }
    }

    class DomainLociComparator implements Comparator<MapData> {

        @Override
        public int compare(MapData o1, MapData o2) {
            int r = o1.getStartPos() - o2.getStartPos();
            if( r!=0 ) {
                return r;
            }
            r = o1.getStopPos() - o2.getStopPos();
            if( r!=0 ) {
                return r;
            }
            r = o1.getChromosome().compareTo(o2.getChromosome());
            if( r!=0 ) {
                return r;
            }
            return o1.getStrand().compareTo(o2.getStrand());
        }
    }

    public void load(UniProtRecord rec) throws Exception {

        handleXdbIds(rec);
    }

    void handleXdbIds(UniProtRecord rec) throws Exception {
        // from incoming data remove those xdb acc ids from trembl that are the same as xdb acc ids from sprot
        int removedTremblDuplicates = rec.removeSprotTremblDuplicates();
        if( removedTremblDuplicates>0 ) {
            UniProtDataLoadManager.getInstance().incrementCounter("skipped TREMBL duplicates", removedTremblDuplicates);
        }
        List<XdbId> xdbIdsIncoming = new ArrayList<>(rec.getXdbIds());

        // get uniprot xdb ids from RGD
        List<XdbId> rgdXdbIds = dbDao.getUniProtXdbIds(rec.getRgdId(), speciesTypeKey);

        // construct a set of duplicates: xdbids that are present both in rgd and in incoming list
        List<XdbId> dupXdbs = new ArrayList<>(rgdXdbIds);
        dupXdbs.retainAll(xdbIdsIncoming);
        // remove from incoming xdbids those that are already in rgd
        xdbIdsIncoming.removeAll(rgdXdbIds);
        // remove from rgdids the duplicates (present on incoming xdbid list)
        // but remove only one copy if there are more identical duplicates
        removeDuplicates(rgdXdbIds, dupXdbs);

        // real work
        if( xdbIdsIncoming.size()>0 ) {
            dbDao.insertXdbIds(xdbIdsIncoming, RgdId.OBJECT_KEY_GENES);
        }
        if( rgdXdbIds.size() > 0 ) {
            dbDao.deleteXdbIds(rgdXdbIds, RgdId.OBJECT_KEY_GENES);
        }

        // update nr of untouched rows
        dbDao.updateLastModificationDate(dupXdbs);
    }

    // similar to rgdXdbIds.removeAll(dupXdbs) except that it removes only
    // one copy from rgdXdbIds table if duplicate is found
    // return nr of removed duplicates
    int removeDuplicates(List<XdbId> rgdXdbIds, List<XdbId> dupXdbs) {

        int removeCount = 0;
        for( XdbId xdbId: dupXdbs ) {
            if( rgdXdbIds.remove(xdbId) )
                removeCount++;
        }
        return removeCount;
    }

    /**
     * @return Returns the dbCheckDao.
     */
    public UniProtDAO getDbDao() {
        return dbDao;
    }
    /**
     * @param dbDao The dbCheckDao to set.
     */
    public void setDbDao(UniProtDAO dbDao) {
        this.dbDao = dbDao;
    }

    public Map<String,Integer> getActiveXdbIdMap() {
        return activeXdbIdMap;
    }

    public void setActiveXdbIdMap(Map<String,Integer> activeXdbIdMap) {
        this.activeXdbIdMap = activeXdbIdMap;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }
}
