package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Protein;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * @author mtutaj
 * @since 3/27/2019
 * Downloads files with deleted accessions from UniProtKB and withdraws corresponding proteins in RGD.
 */
public class DeletedAccessions {

    private String deletedSprotAccessionsFile;
    private String deletedTremblAccessionsFile;

    UniProtDAO dao = new UniProtDAO();

    Logger log = Logger.getLogger("deleted_acc");

    public void run() throws Exception {
        try {
            _run();
        } catch( Exception e ) {
            Utils.printStackTrace(e, log);  // log the exception and rethrow it
            throw e;
        }
    }

    void _run() throws Exception {

        long time0 = System.currentTimeMillis();
        log.info("=== Starting module: Protein Deleted Accessions");

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        log.info("    "+sdt.format(new Date(time0)));

        List<String> delAccList = new ArrayList<>();

        // download files
        FileDownloader fd = new FileDownloader();
        fd.setUseCompression(true);
        fd.setPrependDateStamp(true);

        fd.setExternalFile(getDeletedSprotAccessionsFile());
        fd.setLocalFile("data/delacc_sp.txt.gz");
        String localFile = fd.downloadNew();
        log.info("  downloaded 'delacc_sp.txt'");
        parseFile(localFile, delAccList);

        fd.setExternalFile(getDeletedTremblAccessionsFile());
        fd.setLocalFile("data/delacc_tr.txt.gz");
        localFile = fd.downloadNew();
        log.info("  downloaded 'delacc_tr.txt'");
        parseFile(localFile, delAccList);

        // process all species
        log.info("    total accessions read:  "+delAccList.size());

        Collection<Integer> speciesTypeKeys = SpeciesType.getSpeciesTypeKeys();
        speciesTypeKeys.parallelStream().forEach( speciesTypeKey -> {
            processForSpecies(speciesTypeKey, delAccList);
        });

        log.info("=== Ending module: Protein Deleted Accessions --- elapsed "+ Utils.formatElapsedTime(time0, System.currentTimeMillis()));
        log.info("");
    }

    void parseFile(String fileName, List<String> delAccList) throws IOException {

        // a valid accession contains only uppercase letters and numbers
        int linesRead = 0;
        int accRead = 0;

        BufferedReader in = Utils.openReader(fileName);
        String line;
        while( (line=in.readLine())!=null ) {

            linesRead++;
            if( line.matches("[A-Z][A-Z0-9]+") ) {
                delAccList.add(line);
                accRead++;
            }
        }
        in.close();

        log.info("    file lines read:  "+linesRead);
        log.info("    accessions read:  "+accRead);
    }

    void processForSpecies(int speciesTypeKey, List<String> delAccList) {

        Map<String, Integer> accIdsInRgd = getAccIdsInRgd(speciesTypeKey);
        int activeProteinsInRgd = accIdsInRgd.size();
        if( activeProteinsInRgd==0 ) {
            return;
        }

        String speciesName = SpeciesType.getCommonName(speciesTypeKey);

        int toBeWithdrawnCount = 0;
        for( String acc: delAccList ) {

            Integer proteinRgdId = accIdsInRgd.get(acc);
            if( proteinRgdId!=null ) {
                toBeWithdrawnCount++;

                try {
                    dao.withdrawProtein(proteinRgdId);
                } catch(Exception e) {
                    throw new RuntimeException(e);
                }

                log.debug("withdrawing protein "+acc+", RGD:"+proteinRgdId+", for "+speciesName);
            }
        }

        log.info(speciesName+": active proteins in RGD: "+activeProteinsInRgd);
        log.info(speciesName+": withdrawn proteins: "+toBeWithdrawnCount);
    }

    Map<String, Integer> getAccIdsInRgd(int speciesTypeKey) {

        if( speciesTypeKey==0 ) {
            return Collections.EMPTY_MAP;
        }

        try {
            List<Protein> activeProteinsInRgd = dao.getActiveProteins(speciesTypeKey);
            if( activeProteinsInRgd.isEmpty() ) {
                return Collections.EMPTY_MAP;
            }
            Map<String, Integer> accIdsInRgd = new HashMap<>();

            for (Protein p : activeProteinsInRgd) {
                accIdsInRgd.put(p.getUniprotId(), p.getRgdId());
            }
            return accIdsInRgd;

        } catch(Exception e) {
            throw new RuntimeException(e);
        }
    }

    public void setDeletedSprotAccessionsFile(String deletedSprotAccessionsFile) {
        this.deletedSprotAccessionsFile = deletedSprotAccessionsFile;
    }

    public String getDeletedSprotAccessionsFile() {
        return deletedSprotAccessionsFile;
    }

    public void setDeletedTremblAccessionsFile(String deletedTremblAccessionsFile) {
        this.deletedTremblAccessionsFile = deletedTremblAccessionsFile;
    }

    public String getDeletedTremblAccessionsFile() {
        return deletedTremblAccessionsFile;
    }
}
