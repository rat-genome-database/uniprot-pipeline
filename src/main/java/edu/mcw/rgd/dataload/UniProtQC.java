package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since 6/18/15
 * QC module: matches incoming records with gene rgd ids
 */
public class UniProtQC {

    private UniProtDAO dao;
    private int speciesTypeKey;
    private int unMatched;
    private int inActiveGene;
    private int newActiveGene;
    private List<String> uniprotSources = new ArrayList<>(2);
    private Map<String, Integer> matchXdbCount = new TreeMap<>();

    static Logger logMain = LogManager.getLogger("main");

    public UniProtQC() {
        uniprotSources.add(UniProtDAO.SWISSPROT);
        uniprotSources.add(UniProtDAO.TREMBL);
    }

    public void qc(List<UniProtRatRecord> incomingRecords) throws Exception {
        for( UniProtRatRecord rec: incomingRecords ) {
            if( qc(rec) ) {

                if(rec.getUniProtAccId().equals("A0A0H2UH92") ) {
                    // in Oct 2019, Stan agreed, that protein A0A0H2UH92 has been mistakenly associated with gene Slc25a16 RGD:1311311
                    // instead of Dna2 RGD:1306791
                    // so we override this association here
                    if( rec.matchingRgdIds.contains(1311311) ) {
                        int index = rec.matchingRgdIds.indexOf(1311311);
                        rec.matchingRgdIds.set(index, 1306791);
                        System.out.println("A0A0H2UH92 override: replaced gene association RGD:1311311 -> RGD:1306791");
                    }
                }

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

        // Does the record have a match by gene symbol?
        if( matchByGeneSymbol(data, "10. ") ) {
            return true;
        }

        // no match by NCBI Gene ID, no match on UniProt ID and no match on RefSeq ID
        unMatched++;
        incrementMatchCount("11. unmatched");
        // no match for NCBI Gene ID, HGNC ID, MGI ID, UniProt ID, RefSeq ID, Ensembl ID or UniProt Gene Name
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
            }
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
            logMain.warn("CONFLICT: RGD ID "+rgdid+" is for different species!");
            return 0;
        }

        if( isActive ) {
            // LOG file: NEW_ACTIVE_RGDID
            if( isHistory ) {
                newActiveGene++;
                // rgdid+" is new active RGD_ID retrieved from RGD_ID_HISTORY
            }

            return rgdid;
        }

        // rgd-id is not active -- check the history
        // Has the inactive gene been replaced by an new one?
        int newRgdId = dao.checkHistory(rgdid);
        // no active rgd_id found in rgd_id_history table
        if( newRgdId==0 ) {
            inActiveGene++;
            // rgdid+" is not active and it does not have rgd id history
            return 0;
        }
        else {
            // check if rgdid pulled from history is active, and if not active check its history as well...
            // "inactive "+rgdid+" has rgd id history "+newRgdId
            return getActiveRgdId(newRgdId, true, rec);
        }
    }

    boolean matchByGeneSymbol(UniProtRatRecord rec, String matchPriority) throws Exception {

        if( rec.getGeneName()==null ) {
            return false;
        }

        int activeRgdIdCount = 0;

        for( Gene gene: dao.getGenesBySymbol(rec.getGeneName(), getSpeciesTypeKey()) ) {
            int activeRgdId = getActiveRgdId(gene.getRgdId(), false, rec);
            if( activeRgdId!=0 ) {
                rec.matchingRgdIds.add(activeRgdId);
                activeRgdIdCount++;
            }
        }

        if( activeRgdIdCount>1 ) {
            //dbLog.addLogProp("multiple GNs: "+ Utils.concatenate(rec.matchingRgdIds,","), "MULTI-"+rec.getGeneName(), rec.getRecNo(), PipelineLogger.REC_FLAG, rec.makeRecAsString());
        }

        if( activeRgdIdCount>0 ) {
            incrementMatchCount(matchPriority + "UniProtGeneName");
            return true;
        } else {
            return false;
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

    public void dumpMatchSummary() {
        logMain.info("== MATCH SUMMARY (by match priority) ==");
        for( Map.Entry<String, Integer> entry: matchXdbCount.entrySet() ) {
            logMain.info("  "+entry.getKey()+": "+entry.getValue() );
        }
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
