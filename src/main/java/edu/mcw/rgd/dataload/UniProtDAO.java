package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.*;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.IntStringMapQuery;
import edu.mcw.rgd.dao.spring.StringListQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.PipelineLogger;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 * All code dealing with database is found here -- to centralize database access
 */
public class UniProtDAO extends AbstractDAO {

    public final static String TREMBL = "/TrEMBL";
    public final static String SWISSPROT = "/Swiss-Prot";
    public final static int OBJECT_KEY_PROTEIN_DOMAINS = 23;

    static Logger logInsertedIds = Logger.getLogger("ids_inserted");
    static Logger logDeletedIds = Logger.getLogger("ids_deleted");
    static Logger logInsertedProteins = Logger.getLogger("proteins_inserted");
    static Logger logUpdatedProteins = Logger.getLogger("proteins_updated");
    static Logger logProteinsXdbIds = Logger.getLogger("proteins_xdb_ids");
    static Logger logAssociations = Logger.getLogger("associations");
    static Logger logAliases = Logger.getLogger("aliases");
    static Logger logSequences = Logger.getLogger("sequences");

    private AliasDAO aliasDAO = new AliasDAO();
    private AssociationDAO associationDAO = new AssociationDAO();
    private GenomicElementDAO geDAO = new GenomicElementDAO();
    private ProteinDAO proteinDAO = new ProteinDAO();
    private RGDManagementDAO rgdidsDAO = new RGDManagementDAO();
    private SequenceDAO seqDAO = new SequenceDAO();
    private XdbIdDAO xdbidDAO = new XdbIdDAO();

    private int rowsDeleted;
    private int rowsInserted;
    private int rowsMatched;
    private Date processingStartTime;
    private int staleRowsDeleted;
    private int aliasesInserted;
    private int aliasesUpToDate;

    /**
     * return true if RGD_ID is active; false otherwise
     * @param rgdId  gene rgd id
     * @return true if RGD_ID is active; false otherwise
     * @throws Exception
     */
    public Boolean isGeneActive(int rgdId, int speciesTypeKey) throws Exception {

        String cacheKey = rgdId+"|"+speciesTypeKey;
        Boolean isActive = _cacheObjectStatus.get(cacheKey);
        if( isActive!=null )
            return isActive;

        RgdId id = rgdidsDAO.getRgdId2(rgdId);
        if( id.getSpeciesTypeKey()!=speciesTypeKey )
            return null;
        isActive = id.getObjectStatus().equalsIgnoreCase("ACTIVE");
        _cacheObjectStatus.put(cacheKey, isActive);
        return isActive;
    }

    // cache of gene statuses
    static Map<String,Boolean> _cacheObjectStatus = new HashMap<>();

    /**
     * return new RGD_ID from rgd_id_history, or 0, if there is no history
     * @param rgdid rgd id of object to check for history
     * @return new RGD_ID from rgd_id_history, or 0, if there is no history
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int checkHistory(int rgdid) throws Exception {

        return rgdidsDAO.getRgdIdFromHistory(rgdid);
    }

    // return list of RGD_IDs matching given XDB_KEY and ACC_ID; exclude genes of type 'splice' or 'allele'
    // exclude matching by xdb ids imported by UniProtKB pipeline
    public List<Integer> matchXdbId(String accId, int xdbKey, int speciesTypeKey, List<String> uniprotSources) throws Exception {

        return xdbidDAO.getGeneRgdIdsByXdbId(xdbKey, accId, speciesTypeKey, uniprotSources);
    }

    /**
     * get xdb ids for given rgd id that were loaded by UniProtKB pipeline
     * @param rgdId rgd id
     * @param speciesTypeKey species results are limited to
     * @return a list, possibly empty, of XdbIds
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<XdbId> getUniProtXdbIds(int rgdId, int speciesTypeKey) throws Exception {

        if( uniProtXdbIds.isEmpty() ) {
            loadCacheOfUniProtXdbIds(PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.SWISSPROT, speciesTypeKey);
            loadCacheOfUniProtXdbIds(PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.TREMBL, speciesTypeKey);
        }

        List<XdbId> xdbIds = uniProtXdbIds.get(rgdId);
        return xdbIds!=null ? xdbIds : Collections.<XdbId>emptyList();
    }
    static Map<Integer, List<XdbId>> uniProtXdbIds = new HashMap<>();


    private void loadCacheOfUniProtXdbIds(String srcPipeline, int speciesTypeKey) throws Exception {
        XdbId filter = new XdbId();
        filter.setSrcPipeline(srcPipeline);
        for( XdbId xdbId: xdbidDAO.getXdbIds(filter, speciesTypeKey) ) {
            List<XdbId> xdbIds = uniProtXdbIds.get(xdbId.getRgdId());
            if( xdbIds==null ) {
                xdbIds = new ArrayList<>();
                uniProtXdbIds.put(xdbId.getRgdId(), xdbIds);
            }
            xdbIds.add(xdbId);
        }
    }

    public Map<Integer, List<XdbId>> getUniProtSecondaryIds(int speciesType, int objectKey) throws Exception {
        XdbId filter = new XdbId();
        filter.setXdbKey(XdbId.XDB_KEY_UNIPROT_SECONDARY);

        Map<Integer, List<XdbId>> secondaryIdMap = new HashMap<>();
        for( XdbId id: xdbidDAO.getXdbIds(filter, speciesType, objectKey) ){
            List<XdbId> list = secondaryIdMap.get(id.getRgdId());
            if( list==null ) {
                list = new ArrayList<>();
                secondaryIdMap.put(id.getRgdId(), list);
            }
            list.add(id);
        }
        return secondaryIdMap;
    }

    public List<XdbId> getRefSeqIdsForProteins(String srcPipeline) throws Exception {
        XdbId filter = new XdbId();
        filter.setXdbKey(XdbId.XDB_KEY_GENEBANKPROT);
        filter.setSrcPipeline(srcPipeline);

        return xdbidDAO.getXdbIds(filter, SpeciesType.ALL, RgdId.OBJECT_KEY_PROTEINS);
    }

    /**
     * insert a number of XdbIds into database
     * @param xdbIds a list of XdbId objects to insert
     * @param objectKey object key of object whose XdbIds are being inserted
     * @return count of rows affected
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int insertXdbIds(List<XdbId> xdbIds, int objectKey) throws Exception {

        // insert uniprot ids
        int rowsAffected = xdbIds.size();
        if( rowsAffected>0 ) {
            xdbidDAO.insertXdbs(xdbIds);
            this.rowsInserted += rowsAffected;
        }

        String species = "";
        if( UniProtDataLoadManager.getInstance()!=null ) {
            species = "|species:"+UniProtDataLoadManager.getInstance().getSpeciesTypeKey();
        }

        // dump inserted uniprot ids
        for( XdbId xdbId: xdbIds ) {

            String msg = "INSERT|"+xdbId.dump("|")+species;
            if( objectKey==RgdId.OBJECT_KEY_GENES ) {
                logInsertedIds.info(msg);
                if( _insertedRgdIds.add(xdbId.getRgdId()) ) {
                    logInsertedIds.debug("INSERTED "+xdbId.getRgdId());
                }
            } else if( objectKey==RgdId.OBJECT_KEY_PROTEINS ) {
                logProteinsXdbIds.info(msg);
            }
        }

        return rowsAffected;
    }
    static Set<Integer> _insertedRgdIds = new HashSet<>();

    /**
     * delete a number of XdbIds from database
     * @param xdbIds a list of XdbId objects to delete
     * @param objectKey object key of object whose XdbIds are being deleted
     * @return count of rows affected
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int deleteXdbIds(List<XdbId> xdbIds, int objectKey) throws Exception {

        // spring framework returns wrong numbers of rows affected when doing inserts in batches;
        // could not figure it out why
        int rowsAffected = xdbIds.size();
        if( rowsAffected>0 ) {
            // dump deleted uniprot ids
            for( XdbId xdbId: xdbIds ) {
                String msg = "DELETE|" + xdbId.dump("|") + "|species:" + UniProtDataLoadManager.getInstance().getSpeciesTypeKey();

                if( xdbId.getXdbKey()==XdbId.XDB_KEY_UNIPROT || xdbId.getXdbKey()==XdbId.XDB_KEY_UNIPROT_SECONDARY ) {
                    insertOldProteinIdAlias(xdbId.getRgdId(), xdbId.getAccId());
                }

                if( objectKey==RgdId.OBJECT_KEY_GENES ) {
                    logDeletedIds.info(msg);
                    if (_deletedRgdIds.add(xdbId.getRgdId())) {
                        logDeletedIds.info("DELETED RGD:" + xdbId.getRgdId());
                    }
                } else if( objectKey==RgdId.OBJECT_KEY_PROTEINS ) {
                    logProteinsXdbIds.info(msg);
                }
            }

            deleteXdbIds(xdbIds);
            rowsDeleted += rowsAffected;
        }
        return rowsAffected;
    }
    static Set<Integer> _deletedRgdIds = new HashSet<>();

    void deleteXdbIds(List<XdbId> xdbIds) throws Exception {
        // dump deleted uniprot ids
        for( XdbId xdbId: xdbIds ) {
            if( xdbId.getXdbKey()==XdbId.XDB_KEY_UNIPROT || xdbId.getXdbKey()==XdbId.XDB_KEY_UNIPROT_SECONDARY ) {
                insertOldProteinIdAlias(xdbId.getRgdId(), xdbId.getAccId());
            }
        }

        xdbidDAO.deleteXdbIds(xdbIds);
    }

    synchronized void insertOldProteinIdAlias(int rgdId, String accId) throws Exception {

        // check if the alias to be deleted is already in RGD
        List<Alias> aliasesInRgd = aliasDAO.getAliases(rgdId, "old_protein_id");
        for( Alias alias: aliasesInRgd ) {
            if( alias.getValue().equals(accId) ) {
                logAliases.info("UP_TO_DATE "+alias.dump("|"));
                aliasesUpToDate++;
                return;
            }
        }

        Alias alias = new Alias();
        alias.setRgdId(rgdId);
        alias.setTypeName("old_protein_id");
        alias.setValue(accId);
        alias.setNotes("created by UniProtKB pipeline");
        aliasDAO.insertAlias(alias);
        logAliases.info("INSERTED "+alias.dump("|"));
        aliasesInserted++;
    }

    /**
     * update last modification date for a number of XdbIds
     * @param xdbIds a list of XdbId objects with last modification date to be updates
     * @return count of rows affected
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int updateLastModificationDate(List<XdbId> xdbIds) throws Exception {

        List<Integer> accXdbKeys = new ArrayList<>(xdbIds.size());
        for( XdbId id: xdbIds ) {
            // do not update rows that already have been processed by the pipeline during this run
            if( id.getModificationDate()==null || id.getModificationDate().compareTo(this.getProcessingStartTime()) < 0 )
                accXdbKeys.add(id.getKey());
        }

        int rowsAffected = accXdbKeys.size();
        if( rowsAffected>0 ) {
            xdbidDAO.updateModificationDate(accXdbKeys);
            rowsMatched += rowsAffected;
        }
        return rowsAffected;
    }

    /**
     * Delete all stale rows from RGD_ACC_XDB table.
     * <p>
     * Note: stale rows are rows that had not been 'touched' by the pipeline during the current run,
     * but that had been inserted by previous runs of the pipeline.
     * @param speciesTypeKey species type key
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void deleteStaleRows(int speciesTypeKey) throws Exception {

        // delete stale rows for SPROT pipeline
        int rows = deleteStaleRows(PipelineLogger.PIPELINE_UNIPROT+SWISSPROT, speciesTypeKey);
        UniProtDataLoadManager.getInstance().incrementCounter("stale rows deleted for SWISSPROT", rows);

        // delete stale rows for TREMBL pipeline
        rows = deleteStaleRows(PipelineLogger.PIPELINE_UNIPROT+TREMBL, speciesTypeKey);
        UniProtDataLoadManager.getInstance().incrementCounter("stale rows deleted for TREMBL", rows);
    }

    private int deleteStaleRows(String srcPipeline, int speciesTypeKey) throws Exception {

        List<XdbId> staleXdbIds = xdbidDAO.getXdbIdsModifiedBefore(this.getProcessingStartTime(), srcPipeline, speciesTypeKey);
        for( XdbId xdbId: staleXdbIds ) {
            logDeletedIds.info("DELETE_STALE|"+xdbId.dump("|")+"|species:"+UniProtDataLoadManager.getInstance().getSpeciesTypeKey());
            if( _deletedStaleRgdIds.add(xdbId.getRgdId()) ) {
                logDeletedIds.debug("DELETED STALE RGD:"+xdbId.getRgdId());
            }
        }

        // spring framework returns wrong numbers of rows affected when doing inserts in batches;
        // could not figure it out why
        int rowsAffected = staleXdbIds.size();
        if( rowsAffected>0 ) {
            deleteXdbIds(staleXdbIds);
            staleRowsDeleted += rowsAffected;
        }
        return rowsAffected;
    }
    static Set<Integer> _deletedStaleRgdIds = new HashSet<>();

    public int getRowsDeleted() {
        return rowsDeleted;
    }

    public void setRowsDeleted(int rowsDeleted) {
        this.rowsDeleted = rowsDeleted;
    }

    public int getRowsInserted() {
        return rowsInserted;
    }

    public void setRowsInserted(int rowsInserted) {
        this.rowsInserted = rowsInserted;
    }

    public int getRowsMatched() {
        return rowsMatched;
    }

    public void setRowsMatched(int rowsMatched) {
        this.rowsMatched = rowsMatched;
    }

    public Date getProcessingStartTime() {
        return processingStartTime;
    }

    public void setProcessingStartTime(Date processingStartTime) {
        this.processingStartTime = processingStartTime;
    }

    public int getStaleRowsDeleted() {
        return staleRowsDeleted;
    }

    public int getAliasesInserted() {
        return aliasesInserted;
    }

    public int getAliasesUpToDate() {
        return aliasesUpToDate;
    }

    /////////////////
    /// ASSOCIATION DAO
    /////////////////

    /**
     * return list of all association for given master rgd id and association type
     * @param masterRgdId master rgd id
     * @param assocType association type
     * @return list of Association objects; never null, but returned list could be empty
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<Association> getAssociationsForMasterRgdId(int masterRgdId, String assocType) throws Exception {
        return associationDAO.getAssociationsForMasterRgdId(masterRgdId, assocType);
    }

    public void insertAssociations(Collection<Association> assocsForInsert) throws Exception {

        for( Association assoc: assocsForInsert ) {
            associationDAO.insertAssociation(assoc);
            logAssociations.info("INSERT|"+assoc.dump("|"));
        }
    }

    public void deleteAssociations(Collection<Association> assocsForDelete) throws Exception {

        for( Association assoc: assocsForDelete ) {
            associationDAO.deleteAssociationByKey(assoc.getAssocKey());
            logAssociations.info("DELETE|"+assoc.dump("|"));
        }
    }

    /////////////////
    /// PROTEIN DAO
    /////////////////

    /**
     * create a new RGD_ID, and insert PROTEIN object
     * @param protein Protein object
     * @return Protein object with newly generated RGD_ID
     * @throws Exception
     */
    public Protein insertProtein(Protein protein) throws Exception {
        RgdId id = rgdidsDAO.createRgdId(RgdId.OBJECT_KEY_PROTEINS, "ACTIVE", "", protein.getSpeciesTypeKey());
        protein.setRgdId(id.getRgdId());

        proteinDAO.insertProtein(protein);
        logInsertedProteins.info(protein.dump("|"));
        return protein;
    }

    public void updateProteinUniProtId(Protein protein, String oldUniProtId) throws Exception {
        logUpdatedProteins.info("RGD:"+protein.getRgdId()+"|OLD_UNIPROT_ID:"+oldUniProtId+"|NEW_UNIPROT_ID:"+protein.getUniprotId());
        proteinDAO.updateProtein(protein);
    }

    public Protein getProteinByUniProtId(String uniProtId) throws Exception {
        return proteinDAO.getProteinByUniProtId(uniProtId);
    }

    public GenomicElement getProteinDomainObject(String proteinDomain) throws Exception {
        List<GenomicElement> list = geDAO.getElementsBySymbol(proteinDomain, OBJECT_KEY_PROTEIN_DOMAINS);
        return list.isEmpty() ? null : list.get(0);
    }

    public synchronized GenomicElement insertDomainName(String proteinDomain) throws Exception {
        // ensure the domain name is not in db
        GenomicElement ge = getProteinDomainObject(proteinDomain);
        if( ge!=null ) {
            return ge;
        }

        RgdId id = rgdidsDAO.createRgdId(OBJECT_KEY_PROTEIN_DOMAINS, "ACTIVE", null, 0);
        ge = new GenomicElement();
        ge.setRgdId(id.getRgdId());
        ge.setSymbol(proteinDomain);
        ge.setName(proteinDomain);
        ge.setSoAccId("SO:0000417");
        ge.setObjectType("protein domain");
        ge.setSource("UniProtKB");
        geDAO.insertElement(ge);
        return ge;
    }

    public List<Sequence2> getObjectSequences(int rgdId, String seqType) throws Exception {
        return seqDAO.getObjectSequences2(rgdId, seqType);
    }

    public int insertSequence(Sequence2 seq) throws Exception {
        int r = seqDAO.insertSequence(seq);
        logSequences.info("INSERTED "+seq.dump("|"));
        return r;
    }


    public String getMD5ForObjectSequences(int rgdId) throws Exception {
        return _rgdId2md5.get(rgdId);
    }

    public void loadMD5ForProteinSequences(int speciesTypeKey, String seqType) throws Exception {

        _rgdId2md5.clear();

        String query = "SELECT s.rgd_id,s.seq_data_md5 FROM rgd_sequences s,rgd_ids r "+
                "WHERE s.rgd_id=r.rgd_id AND s.seq_type=? AND r.object_key="+RgdId.OBJECT_KEY_PROTEINS+" AND r.species_type_key=?";
        for( IntStringMapQuery.MapPair pair: IntStringMapQuery.execute(this, query, seqType, speciesTypeKey) ) {

            String prevMD5 = _rgdId2md5.put(pair.keyValue, pair.stringValue);
            if( prevMD5!=null ) {
                throw new Exception("ERROR: multiple sequences in RGD for protein RGDID:"+pair.keyValue);
            }
        }
        System.out.println("  --loaded md5 for "+SpeciesType.getCommonName(speciesTypeKey)+" and seq_type="+seqType+": "+_rgdId2md5.size());
    }
    static Map<Integer, String> _rgdId2md5 = new HashMap<>();
}
