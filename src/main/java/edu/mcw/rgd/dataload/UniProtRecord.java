package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.XdbId;

import java.util.*;

/**
 * @author mtutaj
 * @since 6/16/15
 * <p>
 * represents external db information, downloaded from UniProtKB, and associated with given gene
 */
public class UniProtRecord {

    private Map<String,Integer> activeXdbIdMap;
    private int rgdId;
    private List<XdbId> xdbIds = new ArrayList<>();
    List<ProteinDomain> domains = new ArrayList<>();

    public UniProtRecord(Map<String,Integer> activeXdbIdMap) {
        this.activeXdbIdMap = activeXdbIdMap;
    }

    public void init(int rgdId, UniProtRatRecord rec) throws Exception {
        setRgdId(rgdId);

        // build list of XdbIds (not duplicates) ready for insertion
        Map<String, List<String>> xdbMap = rec.getXdbMap();
        for( Map.Entry<String, List<String>> entry: xdbMap.entrySet() ) {
            String xdbName = entry.getKey();
            int xdbKey = getXdbKeyFromXdbName(xdbName);
            if( xdbKey>0 ) {
                List<String> xdbList = entry.getValue();
                // on list, process by pairs, first goes the xid, then link-info
                for( int i=0; i<xdbList.size(); i+=2 ) {
                    XdbId xdbId = addXdbId(rgdId, xdbKey, xdbList.get(i), xdbList.get(i+1), rec.getUniProtAccId(), rec.getSrcPipeline());
                    if( !xdbIds.contains(xdbId) ) {
                        xdbIds.add(xdbId);
                    }
                }
            }
        }

        domains.addAll(rec.domains);
    }

    // convert external database name, as found in the source file, to XDB_KEY value, as found in RGD_XDB table
    int getXdbKeyFromXdbName(String xdbName) throws Exception {

        Integer xdbKey = activeXdbIdMap.get(xdbName);
        if( xdbKey!=null )
            return xdbKey;
        else {
            return -1; // xbd name not found in mappings
        }
    }

    // set up properties of new XdbId object (never null)
    XdbId addXdbId(int rgdId, int xdbkey, String accid, String textlink, String uniProtAccId, String srcPipeline) throws Exception {

        XdbId xdbId = new XdbId();
        xdbId.setAccId(accid.trim());
        xdbId.setRgdId(rgdId);
        xdbId.setXdbKey(xdbkey);
        if(textlink==null)
            xdbId.setLinkText(accid);
        else
            xdbId.setLinkText(textlink);
        xdbId.setSrcPipeline(srcPipeline);

        Date currentDate = new Date();
        xdbId.setCreationDate(currentDate);
        xdbId.setModificationDate(currentDate);
        xdbId.setNotes(uniProtAccId);
        return xdbId;
    }

    // remove from incoming uniprot id those trembl ids that have the corresponding sprot ids
    // return count of removed trembl entries that were identical as sprot entries
    public int removeSprotTremblDuplicates() {

        int tremblSprotDuplicates = 0;

        // sort entries by xdb_key, acc_id and src_pipeline for easy detection of duplicates
        // note: TrEMBL entries will be always *after* Swiss-Prot entries
        Collections.sort(xdbIds, new Comparator<XdbId>() {
            @Override
            public int compare(XdbId o1, XdbId o2) {
                int r = o1.getXdbKey() - o2.getXdbKey();
                if( r!=0 )
                    return r;
                r = o1.getAccId().compareTo(o2.getAccId());
                if( r!=0 )
                    return r;
                return o1.getSrcPipeline().compareTo(o2.getSrcPipeline());
            }
        });

        // remove duplicates, starting from last but one entry
        for( int i=xdbIds.size()-2; i>=0; i-- ) {
            XdbId x1 = xdbIds.get(i);
            XdbId x2 = xdbIds.get(i+1);

            if( x1.getXdbKey()==x2.getXdbKey()
               && x1.getAccId().equals(x2.getAccId())
               && x1.getSrcPipeline().endsWith(UniProtDAO.SWISSPROT)
               && x2.getSrcPipeline().endsWith(UniProtDAO.TREMBL) ) {

                // duplicate trembl entry detetected
                xdbIds.remove(i+1);
                tremblSprotDuplicates++;
            }
        }
        return tremblSprotDuplicates;
    }

    public int getRgdId() {
        return rgdId;
    }

    public void setRgdId(int rgdId) {
        this.rgdId = rgdId;
    }

    public List<XdbId> getXdbIds() {
        return xdbIds;
    }
}
