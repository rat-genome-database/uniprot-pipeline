package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.ProteinDAO;
import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.PipelineLogger;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;

import java.util.*;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author mtutaj
 * @since 6/18/15
 * QC module: matches incoming records with gene rgd ids
 */
public class UniProtQC {

    private PipelineLogger dbLog = PipelineLogger.getInstance();
    private UniProtDAO dao;
    private int speciesTypeKey;
    private int unMatched;
    private int inActiveGene;
    private int newActiveGene;
    private List<String> uniprotSources = new ArrayList<>(2);
    private Map<String, Integer> matchXdbCount = new TreeMap<>();

    public UniProtQC() {
        uniprotSources.add(PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.SWISSPROT);
        uniprotSources.add(PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.TREMBL);
    }

    public void qc(List<UniProtRatRecord> incomingRecords) throws Exception {
        for( UniProtRatRecord rec: incomingRecords ) {
            if( qc(rec) ) {
                qcProteinDomains(rec);
                qcProteinSequence(rec);
            }
        }
    }

    /** perform qc on the incoming data; the goal is to determine a list of active rgd ids
     * that will be associated with the incoming data
     * 
     * @param data incoming data
     * @throws Exception
     */
    public boolean qc(UniProtRatRecord data) throws Exception {

        // continue until given stage of processing will give at least one active RGD_ID;

        // Does the record have any NCBI Gene ID?
        if( matchByXdbId(data, "GeneID", XdbId.XDB_KEY_NCBI_GENE, "GeneID", " 1. ") )
            return true;

        // Does the record have any HGNC ID? (will work for human only)
        if( matchByXdbId(data, "HGNC", XdbId.XDB_KEY_HGNC, "HGNC", " 2. ") )
            return true;

        // Does the record have any MGI ID? (will work for mouse only)
        if( matchByXdbId(data, "MGI", XdbId.XDB_KEY_MGD, "MGI", " 3. ") )
            return true;

        // Does the record have any RefSeq ID?
        if( matchByXdbId(data, "RefSeq", XdbId.XDB_KEY_GENEBANKPROT, "RefSeq", " 4. ") )
            return true;

        // does the record have any UniProtKB/Swiss-Prot ID?
        if( matchByXdbId(data, "UniProt", XdbId.XDB_KEY_UNIPROT, "UniProt", " 5. ") )
            return true;

        // does the record have any UniProtKB/Swiss-Prot ID imported by Gene pipelines as GeneBank protein?
        if( matchByXdbId(data, "UniProt", XdbId.XDB_KEY_GENEBANKPROT, "GeneBankProt", " 6. ") )
            return true;

        // Does the record have any Ensembl ID?
        if( matchByXdbId(data, "EnsemblP", XdbId.XDB_KEY_ENSEMBL_PROTEIN, "EnsemblP", " 7. ") )
            return true;

        // Does the record have any Ensembl ID?
        if( matchByXdbId(data, "EnsemblG", XdbId.XDB_KEY_ENSEMBL_GENES, "EnsemblG", " 8. ") )
            return true;

        // Does the record have any Ensembl ID?
        if( matchByXdbId(data, "EnsemblT", XdbId.XDB_KEY_ENSEMBL_TRANSCRIPT, "EnsemblT", " 9. ") )
            return true;


        // no match on NCBI Gene ID, no match on UniProt ID and no match on RefSeq ID
        unMatched++;
        incrementMatchCount("10. unmatched");
        dbLog.addLogProp("no match for NCBI Gene ID, HGNC ID, MGI ID, UniProt ID, RefSeq ID or Ensembl ID", "UNMATCHED", data.getRecNo(), PipelineLogger.REC_FLAG, data.makeRecAsString());
        return false;
    }

    /**
     * match a record of rat data against given xdb; return true if matching active rgd_id is found
      *@param rec incoming data record
     * @param dataXdbName name of xdb in rat data record
     * @param xdbKey xdb key
     * @param xdbName xdb short name to be used in logs
     * @return true if a matching active rgd_is is found
     */
    boolean matchByXdbId(UniProtRatRecord rec, String dataXdbName, int xdbKey, String xdbName, String matchPriority) throws Exception {

        int activeRgdIdCount = 0;

        // Does the record have any xdb ids of given name?
        List<String> xdbInfo = rec.getXdbInfo(dataXdbName);
        // if no active rgdid's found, run match against XDB IDs
        for( int i=0; i<xdbInfo.size(); i+=2 ) {
            String accId = xdbInfo.get(i);

            // does incoming xdb id has a single match with RGD?
            for( Integer rgdId: dao.matchXdbId(accId, xdbKey, getSpeciesTypeKey(), uniprotSources) ) {
                int activeRgdId = getActiveRgdId(rgdId, false, rec);
                if( activeRgdId!=0 ) {
                    rec.matchingRgdIds.add(activeRgdId);
                    activeRgdIdCount++;
                }
                dbLog.addLogProp(xdbName+":"+accId+"\tRGDID:"+activeRgdId, accId, rec.getRecNo(), PipelineLogger.REC_XDBID);
            }
        }

        if( activeRgdIdCount>1 ) {
            dbLog.addLogProp("multiple "+xdbName+"s: "+ Utils.concatenate(rec.matchingRgdIds,","), "MULTI-"+xdbName, rec.getRecNo(), PipelineLogger.REC_FLAG, rec.makeRecAsString());
        }

        if( activeRgdIdCount>0 ) {
            incrementMatchCount(matchPriority + xdbName);
            return true;
        } else {
            return false;
        }
    }

    /**
     * check if given rgd_id is active, or if it has active rgd_id in history
     * @param rgdid RGD_ID to check
     * @param isHistory true, if rgdid being checked was pulled from RGD_ID_HISTORY
     * @return value of active rgd id or 0 otherwise
     * @throws Exception
     */
    public int getActiveRgdId(int rgdid, boolean isHistory, UniProtRatRecord rec) throws Exception {

        // check if rgd id active
        Boolean isActive = dao.isGeneActive(rgdid, getSpeciesTypeKey());
        if( isActive==null ) {
            System.out.println("CONFLICT: RGD ID "+rgdid+" is for different species!");
            return 0;
        }

        if( isActive ) {
            // LOG file: NEW_ACTIVE_RGDID
            if( isHistory ) {
                newActiveGene++;
                this.dbLog.addLogProp(rgdid+" is new active RGD_ID retrieved from RGD_ID_HISTORY", "NEW_ACTIVE_GENE", rec.getRecNo(), PipelineLogger.REC_FLAG, rec.makeRecAsString());
            }

            return rgdid;
        }

        // rgd-id is not active -- check the history
        // Has the inactive gene been replaced by an new one?
        int newRgdId = dao.checkHistory(rgdid);
        // no active rgd_id found in rgd_id_history table
        if( newRgdId==0 ) {
            inActiveGene++;
            this.dbLog.addLogProp(rgdid+" is not active and it does not have rgd id history", "INACTIVE_GENE", rec.getRecNo(), PipelineLogger.REC_FLAG, rec.makeRecAsString());
            return 0;
        }
        else {
            // check if rgdid pulled from history is active, and if not active check its history as well...
            this.dbLog.addLogProp("inactive "+rgdid+" has rgd id history "+newRgdId, "GENE_WITH_HISTORY", rec.getRecNo(), PipelineLogger.REC_FLAG, rec.makeRecAsString());
            return getActiveRgdId(newRgdId, true, rec);
        }
    }

    void incrementMatchCount(String xdbName) {
        Integer count = matchXdbCount.get(xdbName);
        if( count == null ) {
            count = 1;
        } else {
            count++;
        }
        matchXdbCount.put(xdbName, count);
    }

    void qcProteinSequence( UniProtRatRecord data) throws Exception {
        if( data.proteinSequence!=null ) {
            if( !data.proteinSequence.matches("[a-zA-Z]+") ) {
                throw new Exception("ERROR: Unexpected protein sequence: "+data.proteinSequence);
            }
        }
    }

    static Pattern domainNameCounterPattern = Pattern.compile(" \\d+$");

    void qcProteinDomains( UniProtRatRecord data) throws Exception {
        // many proteins have the same domain appearing multiple times
        // f.e. P58365 protein has Cadherin domain appearing 27 times!
        //   so the domain names are like this: 'Cadherin 1', 'Cadherin 2', ... 'Cadherin 27'
        // therefore we have to strip the last part and keep only the domain name f.e. 'Cadherin'
        for (ProteinDomain pd : data.domains) {
            Matcher m = domainNameCounterPattern.matcher(pd.getDomainName());
            int counterPos = -1;
            while( m.find() ) {
                counterPos = m.start();
            }
            if( counterPos>0 ) {
                pd.setDomainName( pd.getDomainName().substring(0, counterPos));
            }
            //logDomain.info(pd.getDomainName() +" =" +data.getUniProtAccId());

            pd.geInRgd = dao.getProteinDomainObject(pd.getDomainName());

            if( pd.geInRgd!=null ) {
                pd.loci = positionProteinDomain(pd, data.uniProtAccId);
            }
        }
    }

    List<MapData> positionProteinDomain(ProteinDomain pd, String uniProtAccId) throws Exception {

        List<MapData> results = new ArrayList<>();

        // convert 1-based aaPos into 0-based nucPos
        int nucStartPos = (pd.aaStartPos-1)*3; // including
        int lenToGo = (pd.aaStopPos-pd.aaStartPos+1)*3;

        // determine primary assembly
        TranscriptDAO tdao = new TranscriptDAO();
        XdbIdDAO xdao = new XdbIdDAO();
        ProteinDAO pdao = new ProteinDAO();

        int primaryMapKey = MapManager.getInstance().getReferenceAssembly(speciesTypeKey).getKey();
        CdsUtils utils = new CdsUtils(tdao, primaryMapKey);

        // map protein domain to protein
        Protein p = pdao.getProteinByUniProtId(uniProtAccId);

        // map protein to ncbi protein acc ids
        List<XdbId> ncbiProtAccIds = xdao.getXdbIdsByRgdId(XdbId.XDB_KEY_GENEBANKPROT, p.getRgdId());
        // get NCBI transcript given ncbi protein acc ids
        for( XdbId xid: ncbiProtAccIds ) {
            List<Transcript> transcripts = tdao.getTranscriptsByProteinAccId(xid.getAccId());
            for (Transcript tr : transcripts) {
                for (MapData md : tr.getGenomicPositions()) {
                    if( md.getMapKey()!=primaryMapKey ) {
                        continue;
                    }
                    List<CodingFeature> cfs = utils.buildCfList(md);

                    List<MapData> mds;
                    if( md.getStrand().equals("-") ) {
                        mds = handleNegativeStrand(cfs, nucStartPos, lenToGo, pd.geInRgd.getRgdId(), uniProtAccId);
                    } else {
                        mds = handlePositiveStrand(cfs, nucStartPos, lenToGo, pd.geInRgd.getRgdId(), uniProtAccId);
                    }

                    // merge all cds chunks of protein domain position into one range
                    MapData mdRange = null;
                    for( MapData m: mds ) {
                        if( mdRange==null ) {
                            mdRange = m.clone();
                        } else {
                            if( m.getStartPos() < mdRange.getStartPos() ) {
                                mdRange.setStartPos( m.getStartPos() );
                            }
                            if( m.getStopPos() > mdRange.getStopPos() ) {
                                mdRange.setStopPos( m.getStopPos() );
                            }
                        }
                    }
                    if( mdRange!=null ) {
                        results.add(mdRange);
                    }
                }
            }
        }
        return results;
    }

    List<MapData> handlePositiveStrand(List<CodingFeature> cfs, int nucStartPos, int lenToGo, int domainRgdId, String uniProtAccId) {

        List<MapData> results = new ArrayList<>();

        // iterate over CDS features
        for (CodingFeature cf : cfs) {
            if (cf.getFeatureType() == TranscriptFeature.FeatureType.CDS) {
                // we found a CDS
                int cdsLen = cf.getStopPos() - cf.getStartPos() + 1;
                if (nucStartPos >= cdsLen) {
                    // nucStartPos is outside of this CDS
                    nucStartPos -= cdsLen;
                    continue;
                }

                // nucStartPos is within this cds
                //
                MapData mdDomain = new MapData();
                mdDomain.setStartPos(cf.getStartPos() + nucStartPos);
                mdDomain.setChromosome(cf.getChromosome());
                mdDomain.setMapKey(cf.getMapKey());
                mdDomain.setSrcPipeline("UniProtKB");
                mdDomain.setStrand(cf.getStrand());
                mdDomain.setRgdId(domainRgdId);
                results.add(mdDomain);
                mdDomain.setNotes(uniProtAccId+" part "+results.size());

                // is nucStopPos entirely within this CDS
                int nucStopPos = nucStartPos + lenToGo;
                if (nucStopPos < cdsLen) {
                    mdDomain.setStopPos(cf.getStartPos() + nucStopPos);
                    return results;
                } else {
                    // add this part of protein
                    mdDomain.setStopPos(cf.getStopPos());

                    int domainPartLen = mdDomain.getStopPos() - mdDomain.getStartPos() + 1;
                    nucStartPos = 0;
                    lenToGo -= domainPartLen;
                }
            }
        }
        System.out.println("plus strand problem "+uniProtAccId+" "+lenToGo);
        return results;
    }

    List<MapData> handleNegativeStrand(List<CodingFeature> cfs, int nucStartPos, int lenToGo, int domainRgdId, String uniProtAccId) {

        List<MapData> results = new ArrayList<>();

        // iterate over CDS features
        for( int i=cfs.size()-1; i>=0; i-- ) {
            CodingFeature cf = cfs.get(i);
            if (cf.getFeatureType() != TranscriptFeature.FeatureType.CDS ) {
                continue;
            }

            // we found a CDS
            int cdsLen = cf.getStopPos() - cf.getStartPos() + 1;
            if (nucStartPos >= cdsLen) {
                // nucStartPos is outside of this CDS
                nucStartPos -= cdsLen;
                continue;
            }

            MapData mdDomain = new MapData();
            mdDomain.setChromosome(cf.getChromosome());
            mdDomain.setMapKey(cf.getMapKey());
            mdDomain.setSrcPipeline("UniProtKB");
            mdDomain.setStrand(cf.getStrand());
            mdDomain.setRgdId(domainRgdId);
            results.add(mdDomain);
            mdDomain.setNotes(uniProtAccId+" part "+results.size());

            // nucStartPos is within this cds
            // is nucStopPos entirely within this CDS
            int stopPos = cf.getStopPos() - nucStartPos + 1;
            int startPos = stopPos - lenToGo + 1;
            if (startPos >= cf.getStartPos()) {
                // whole part is within this CDS (i.e. part of this CDS
                mdDomain.setStartPos(startPos);
                mdDomain.setStopPos(stopPos);
                return results;
            } else {
                // add this part of protein
                mdDomain.setStartPos(cf.getStartPos());
                mdDomain.setStopPos(stopPos);

                int domainPartLen = mdDomain.getStopPos() - mdDomain.getStartPos() + 1;
                nucStartPos = 0;
                lenToGo -= domainPartLen;
            }
        }
        System.out.println("minus strand problem "+uniProtAccId+" "+lenToGo);
        return results;
    }

    public void dumpMatchSummary() {
        System.out.println("== MATCH SUMMARY (by match priority) ==");
        for( Map.Entry<String, Integer> entry: matchXdbCount.entrySet() ) {
            System.out.println("  "+entry.getKey()+": "+entry.getValue() );
        }
        System.out.println("== MATCH SUMMARY ==");
    }

    public UniProtDAO getDao() {
        return dao;
    }

    public void setDao(UniProtDAO dao) {
        this.dao = dao;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public int getUnMatched() {
        return unMatched;
    }

    public void setUnMatched(int unMatched) {
        this.unMatched = unMatched;
    }

    public int getInActiveGene() {
        return inActiveGene;
    }

    public void setInActiveGene(int inActiveGene) {
        this.inActiveGene = inActiveGene;
    }

    public int getNewActiveGene() {
        return newActiveGene;
    }

    public void setNewActiveGene(int newActiveGene) {
        this.newActiveGene = newActiveGene;
    }
}
