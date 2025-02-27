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

    public UniProtRecord(Map<String,Integer> activeXdbIdMap) {
        this.activeXdbIdMap = activeXdbIdMap;
    }

    public void init(int rgdId, UniProtRatRecord rec) {
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
    }

    // convert external database name, as found in the source file, to XDB_KEY value, as found in RGD_XDB table
    int getXdbKeyFromXdbName(String xdbName) {

        Integer xdbKey = activeXdbIdMap.get(xdbName);
        if( xdbKey!=null )
            return xdbKey;
        else {
            return -1; // xbd name not found in mappings
        }
    }

    // set up properties of new XdbId object (never null)
    XdbId addXdbId(int rgdId, int xdbkey, String accid, String textlink, String uniProtAccId, String srcPipeline) {

        XdbId xdbId = new XdbId();
        xdbId.setAccId(accid.trim());
        xdbId.setRgdId(rgdId);
        xdbId.setXdbKey(xdbkey);
        if( textlink!=null ) {
            xdbId.setLinkText(textlink);
        }
        xdbId.setSrcPipeline(srcPipeline);

        Date currentDate = new Date();
        xdbId.setCreationDate(currentDate);
        xdbId.setModificationDate(currentDate);
        xdbId.setNotes(uniProtAccId);
        return xdbId;
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

    public String fullDump() {
        StringBuffer buf = new StringBuffer();
        buf.append("   RGD:"+rgdId+"\n");
        buf.append("   activeXdbIdMap\n");
        buf.append(activeXdbIdMap.toString());
        buf.append("   xdb ids\n");
        for( XdbId id: xdbIds ) {
            buf.append("     "+id.dump("|")+"\n");
        }
        return buf.toString();
    }
}
