package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.ProteinDAO;
import edu.mcw.rgd.datamodel.Protein;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.commons.collections.CollectionUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by mtutaj on 3/21/2018.
 */
public class Refseq2UniprotLoader {

    UniProtDAO dao = new UniProtDAO();

    private String mappingsFile;
    private String srcPipeline;

    public void run() throws Exception {

        dao.setProcessingStartTime(new Date());

        System.out.println("Refseq2UniprotLoader: load RefSeq protein acc ids from NCBI to protein objects");

        // #NCBI_protein_accession	UniProtKB_protein_accession
        // NP_659298	Q34879
        String localFile = downloadFileWithMappings();

        // load uniProtAcc to protein-rgd-id
        Map<String, Integer> uniProtMap = loadProteinAccId2RgdIds();

        // generate incoming data: NCBI protein accession ids mapped to existing UniProtKB acc ids
        List<XdbId> ncbiProtIds = generateIncomingData(localFile, uniProtMap);
        System.out.println("incoming xdb ids   ="+ncbiProtIds.size());

        // QC
        List<XdbId> ncbiProtIdsInRgd = dao.getRefSeqIdsForProteins(getSrcPipeline());
        System.out.println("xdb ids in RGD     ="+ncbiProtIdsInRgd.size());

        Collection<XdbId> ncbiProtIdsForInsert = CollectionUtils.subtract(ncbiProtIds, ncbiProtIdsInRgd);
        System.out.println("xdb ids for insert ="+ncbiProtIdsForInsert.size());
        Collection<XdbId> ncbiProtIdsForDelete = CollectionUtils.subtract(ncbiProtIdsInRgd, ncbiProtIds);
        System.out.println("xdb ids for delete ="+ncbiProtIdsForDelete.size());
        Collection<XdbId> ncbiProtIdsUpToDate = CollectionUtils.intersection(ncbiProtIdsInRgd, ncbiProtIds);
        System.out.println("xdb ids up-to-date ="+ncbiProtIdsUpToDate.size());

        // DATA LOAD
        if( !ncbiProtIdsForInsert.isEmpty() ) {
            List<XdbId> forInsert = new ArrayList<>(ncbiProtIdsForInsert);
            dao.insertXdbIds(forInsert, RgdId.OBJECT_KEY_PROTEINS);
            System.out.println("xdb ids inserted ="+forInsert.size());
        }

        if( !ncbiProtIdsForDelete.isEmpty() ) {
            List<XdbId> forDelete = new ArrayList<>(ncbiProtIdsForDelete);
            dao.deleteXdbIds(forDelete, RgdId.OBJECT_KEY_PROTEINS);
            System.out.println("xdb ids deleted ="+forDelete.size());
        }

        if( !ncbiProtIdsUpToDate.isEmpty() ) {
            List<XdbId> forUpdate = new ArrayList<>(ncbiProtIdsUpToDate);
            dao.updateLastModificationDate(forUpdate);
            System.out.println("xdb ids up-to-date ="+forUpdate.size());
        }

        System.out.println("=== ELAPSED TIME: "+Utils.formatElapsedTime(System.currentTimeMillis(), dao.getProcessingStartTime().getTime()));
    }

    String downloadFileWithMappings() throws Exception {
        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(getMappingsFile());
        downloader.setLocalFile("data/gene_refseq_uniprotkb.txt.gz");
        downloader.setAppendDateStamp(true);
        return downloader.downloadNew();
    }

    Map<String, Integer> loadProteinAccId2RgdIds() throws Exception {
        System.out.println("loading proteins in RGD");

        // load uniProtAcc to protein-rgd-id
        ConcurrentHashMap<String, Integer> uniProtMap = new ConcurrentHashMap<>();
        List<Protein> proteins = new ProteinDAO().getProteins();
        proteins.parallelStream().forEach(p -> {
            Integer rgdId = uniProtMap.get(p.getUniprotId());
            if( rgdId!=null ) {
                System.out.println("conflict: "+p.getUniprotId()+" "+p.getSymbol()+" "+p.getName());
            }
            uniProtMap.put(p.getUniprotId(), p.getRgdId());
        });
        System.out.println("loaded proteins: "+uniProtMap.size());
        return uniProtMap;
    }

    List<XdbId> generateIncomingData(String mappingsFileName, Map<String, Integer> uniProtMap) throws IOException {
        // generate incoming data: NCBI protein accession ids mapped to existing UniProtKB acc ids
        int matchingLines = 0;
        int skippedLines = 0;
        List<XdbId> ncbiProtAccs = new ArrayList<>();

        BufferedReader in = Utils.openReader(mappingsFileName);
        String line;
        while( (line=in.readLine())!=null ) {
            if( line.startsWith("#") ) {
                continue;
            }
            String[] cols = line.split("[\\t]", -1);
            if( cols.length<2 ) {
                System.out.println("unexpected line format: "+line);
                continue;
            }
            String ncbiProtAcc = cols[0];
            String uniprotAcc = cols[1];

            Integer proteinRgdId = uniProtMap.get(uniprotAcc);
            if( proteinRgdId==null ) {
                skippedLines++;
                continue;
            }

            XdbId id = new XdbId();
            id.setRgdId(proteinRgdId);
            id.setXdbKey(XdbId.XDB_KEY_GENEBANKPROT);
            id.setAccId(ncbiProtAcc);
            id.setSrcPipeline(getSrcPipeline());
            ncbiProtAccs.add(id);
            matchingLines++;
        }

        in.close();

        System.out.println("skipped lines (unknown UniProtKB acc): "+skippedLines);
        System.out.println("matching lines  (known UniProtKB acc): "+matchingLines);

        return ncbiProtAccs;
    }

    public void setMappingsFile(String mappingsFile) {
        this.mappingsFile = mappingsFile;
    }

    public String getMappingsFile() {
        return mappingsFile;
    }

    public void setSrcPipeline(String srcPipeline) {
        this.srcPipeline = srcPipeline;
    }

    public String getSrcPipeline() {
        return srcPipeline;
    }
}
