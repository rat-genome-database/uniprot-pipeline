package edu.mcw.rgd.dataload;

import java.util.*;

import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.XdbId;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 */
public class UniProtDataValidation {
    UniProtDAO dbDao;
    java.util.Map<String,Integer> activeXdbIdMap;
    int speciesTypeKey;
    Logger logMain = LogManager.getLogger("main");
    private Map<Integer, UniProtRecord> records = new HashMap<>();

    /**
     * assign incoming data to matching rgd ids;
     * note: a set of incoming data could be assigned to multiple active rgd ids
     * @param incomingRecords
     */
    public void mergeIncomingRecords(List<UniProtRatRecord> incomingRecords) throws Exception {

        for( UniProtRatRecord incomingRec: incomingRecords ) {

            for( Integer activeRgdId: incomingRec.matchingRgdIds ) {

                UniProtRecord rec = records.get(activeRgdId);
                if( rec==null ) {
                    rec = new UniProtRecord(activeXdbIdMap);
                    records.put(activeRgdId, rec);
                }

                rec.init(activeRgdId, incomingRec);
            }
        }
    }

    public void load() throws Exception {

        logMain.debug(" loading uniprot xdb ids ...");
        for( UniProtRecord rec: records.values() ) {
            load(rec);
        }

        logMain.debug(" loading OK!");
    }

    public void load(UniProtRecord rec) throws Exception {

        handleXdbIds(rec);
    }

    void handleXdbIds(UniProtRecord rec) throws Exception {
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
