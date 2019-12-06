package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Created by mtutaj on 2/11/2016.
 * <p>
 * loads incoming proteins into RGD table PROTEINS including:
 *  protein-to-gene associations, protein sequences, uniprot secondary accession ids
 */
public class ProteinLoader {

    UniProtDataLoadManager manager;
    Logger logMain = Logger.getLogger("main");

    private UniProtDAO dao;
    private int speciesTypeKey;
    private String assocType;
    private String sequenceType;
    private String oldSequenceType;

    public void run(List<UniProtRatRecord> records) throws Exception {

        long time0 = System.currentTimeMillis();

        manager = UniProtDataLoadManager.getInstance();

        getDao().loadMD5ForProteinSequences(getSpeciesTypeKey(), getSequenceType());

        Set<XdbId> incomingSecondaryIds = new HashSet<>();

        for( UniProtRatRecord rec: records ) {

            // get Protein object in RGD, by primary or secondary uniprot id
            Protein protein = getProteinByUniProtId(rec);

            if( protein==null ) {
                // insert new protein object if it was not in RGD yet
                protein = insertProtein(rec);
                manager.incrementCounter("  proteins inserted: ");
            } else {
                manager.incrementCounter("  proteins matched : ");
            }

            handleProteinSequences(rec, protein);

            // handle protein_to_gene associations
            handleProteinToGeneAssociations(rec, protein);

            // handle uniprot secondary accession ids
            handleXdbIds(rec, protein, incomingSecondaryIds);
        }

        handleSecondaryUniProtIds(incomingSecondaryIds);

        logMain.info("PROTEIN LOADER OK,   elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    Protein getProteinByUniProtId(UniProtRatRecord rec) throws Exception {

        // try primary uniprot accession id
        Protein protein = getDao().getProteinByUniProtId(rec.getUniProtAccId());
        if( protein!=null ) {
            return protein;
        }

        // UniProt could have replaced uniprot ids between loads -- try secondary acc ids
        List<String> secondaries = rec.getXdbInfo("UniProtSecondary");
        for( int i=0; i<secondaries.size(); i+=2 ) {
            String uniProtSecondaryId = secondaries.get(i);
            protein = getDao().getProteinByUniProtId(uniProtSecondaryId);
            if( protein!=null ) {
                // protein by secondary accession id -- update uniprot id
                protein.setUniprotId(rec.getUniProtAccId());
                getDao().updateProteinUniProtId(protein, uniProtSecondaryId);
                manager.incrementCounter("  proteins updated : ");
                return protein;
            }
        }
        return null;
    }

    Protein insertProtein(UniProtRatRecord rec) throws Exception {

        // determine protein symbol
        String symbol;
        List<String> uniProtInfo = rec.getXdbInfo("UniProt");
        if( uniProtInfo!=null && uniProtInfo.size()>1 ) {
            symbol = uniProtInfo.get(1);
        } else {
            symbol = rec.getUniProtAccId();
        }

        // create a new RGD id for protein object
        Protein protein = new Protein();
        protein.setSrcPipeline(rec.getSrcPipeline());
        protein.setUniprotId(rec.getUniProtAccId());
        protein.setSymbol(symbol);
        protein.setName(rec.getProteinName());
        protein.setSpeciesTypeKey(getSpeciesTypeKey());
        return getDao().insertProtein(protein);
    }

    void handleProteinSequences(UniProtRatRecord rec, Protein protein) throws Exception {

        Sequence seqIncoming = new Sequence();
        seqIncoming.setRgdId(protein.getRgdId());
        seqIncoming.setSeqData(rec.proteinSequence);
        seqIncoming.setSeqType(getSequenceType());

        String seqInRgdMD5 = getDao().getMD5ForObjectSequences(protein.getRgdId(), getSequenceType());
        if( seqInRgdMD5==null ) {
            getDao().insertSequence(seqIncoming);
            manager.incrementCounter("  protein-sequences inserted", 1);
        } else {
            // see if the protein sequence is the same
            String seqIncomingMD5 = Utils.generateMD5(seqIncoming.getSeqData());
            if( !seqIncomingMD5.equals(seqInRgdMD5) ) {

                // add new sequence
                getDao().insertSequence(seqIncoming);

                // incoming sequence differs from sequence in RGD!
                // downgrade the old sequence to 'old_uniprot_seq'
                getDao().changeSequenceType(seqIncoming, seqInRgdMD5, getOldSequenceType());

                manager.incrementCounter("  protein-sequences updated, old_uniprot_seq created", 1);
            } else {
                manager.incrementCounter("  protein-sequences up-to-date", 1);
            }
        }
    }

    void handleProteinToGeneAssociations(UniProtRatRecord rec, Protein protein) throws Exception {

        // determine type of association
        if( rec.matchingRgdIds.isEmpty() ) {
            manager.incrementCounter("  protein-to-gene assoc match none ");
        } else if( rec.matchingRgdIds.size()==1 ) {
            manager.incrementCounter("  protein-to-gene assoc match single ");
        } else {
            manager.incrementCounter("  protein-to-gene assoc match multiple ");
        }

        // create incoming associations
        List<Association> assocIncoming = new ArrayList<>();
        for( Integer geneRgdId: rec.matchingRgdIds ) {
            Association assoc = new Association();
            assoc.setSrcPipeline(rec.getSrcPipeline());
            assoc.setAssocType(getAssocType()); // 'protein_to_gene'
            assoc.setMasterRgdId(protein.getRgdId());
            assoc.setDetailRgdId(geneRgdId);
            assocIncoming.add(assoc);
        }

        // get existing associations in RGD
        List<Association> assocInRgd = getDao().getAssociationsForMasterRgdId(protein.getRgdId(), getAssocType());

        Collection<Association> assocMatched = CollectionUtils.intersection(assocIncoming, assocInRgd);
        Collection<Association> assocToBeInserted = CollectionUtils.subtract(assocIncoming, assocInRgd);
        Collection<Association> assocToBeDeleted = CollectionUtils.subtract(assocInRgd, assocIncoming);

        getDao().insertAssociations(assocToBeInserted);
        getDao().deleteAssociations(assocToBeDeleted);

        manager.incrementCounter("  protein-to-gene assoc matched ", assocMatched.size());
        manager.incrementCounter("  protein-to-gene assoc inserted ", assocToBeInserted.size());
        manager.incrementCounter("  protein-to-gene assoc deleted ", assocToBeDeleted.size());
    }

    void handleXdbIds(UniProtRatRecord rec, Protein protein, Set<XdbId> incomingSecondaryIds) throws Exception {

        // create incoming xdb ids
        List<String> secondaries = rec.getXdbInfo("UniProtSecondary");
        for( int i=0; i<secondaries.size(); i+=2 ) {
            String uniProtSecondaryId = secondaries.get(i);
            XdbId xdbId = new XdbId();
            xdbId.setAccId(uniProtSecondaryId);
            xdbId.setSrcPipeline(rec.getSrcPipeline());
            xdbId.setRgdId(protein.getRgdId());
            xdbId.setXdbKey(XdbId.XDB_KEY_UNIPROT_SECONDARY);
            incomingSecondaryIds.add(xdbId);
        }
    }

    void handleSecondaryUniProtIds(Collection<XdbId> incomingSecondaryIds) throws Exception {

        List<XdbId> secondaryIdsInRgd = getDao().getUniProtSecondaryIds(getSpeciesTypeKey(), RgdId.OBJECT_KEY_PROTEINS);

        Collection<XdbId> xdbIdsMatched = CollectionUtils.intersection(incomingSecondaryIds, secondaryIdsInRgd);
        Collection<XdbId> xdbIdsToBeInserted = CollectionUtils.subtract(incomingSecondaryIds, secondaryIdsInRgd);
        Collection<XdbId> xdbIdsToBeDeleted = CollectionUtils.subtract(secondaryIdsInRgd, incomingSecondaryIds);

        if( !xdbIdsToBeInserted.isEmpty() ) {
            getDao().insertXdbIds(new ArrayList<>(xdbIdsToBeInserted), RgdId.OBJECT_KEY_PROTEINS);
            manager.incrementCounter("  protein secondary uniprot ids inserted ", xdbIdsToBeInserted.size());
        }
        if( !xdbIdsToBeDeleted.isEmpty() ) {
            getDao().deleteXdbIds(new ArrayList<>(xdbIdsToBeDeleted), RgdId.OBJECT_KEY_PROTEINS);
            manager.incrementCounter("  protein secondary uniprot ids deleted ", xdbIdsToBeDeleted.size());
        }
        manager.incrementCounter("  protein secondary uniprot ids matched ", xdbIdsMatched.size());
    }

    public void setDao(UniProtDAO dao) {
        this.dao = dao;
    }

    public UniProtDAO getDao() {
        return dao;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public void setAssocType(String assocType) {
        this.assocType = assocType;
    }

    public String getAssocType() {
        return assocType;
    }

    public void setSequenceType(String sequenceType) {
        this.sequenceType = sequenceType;
    }

    public String getSequenceType() {
        return sequenceType;
    }

    public void setOldSequenceType(String oldSequenceType) {
        this.oldSequenceType = oldSequenceType;
    }

    public String getOldSequenceType() {
        return oldSequenceType;
    }
}
