package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.Protein;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CanonicalProteins {

    private Map<String,String> canonicalFiles;

    private Logger log = LogManager.getLogger("canonical_proteins");
    private UniProtDAO dao = new UniProtDAO();

    public void run(int speciesTypeKey) throws Exception {

        if( speciesTypeKey<=0 ) {
            for( String species: getCanonicalFiles().keySet() ) {
                run(SpeciesType.parse(species));
            }
        } else {
            String species = SpeciesType.getShortName(speciesTypeKey);
            String fileUrl = getCanonicalFiles().get(species);
            if( fileUrl==null ) {
                log.warn("WARNING: cannot find a file mapping for "+species);
                return;
            }
            run(speciesTypeKey, species, fileUrl);
        }
    }

    void run(int speciesTypeKey, String species, String fileUrl) throws Exception {

        log.info("START for "+species);

        String localFile = downloadFile(species, fileUrl);

        List<String> incomingCanonicalProteins = getIncomingCanonicalProteins(localFile);

        CounterPool counters = new CounterPool();

        List<Protein> activeProteins = dao.getActiveProteins(speciesTypeKey);
        Map<String, Protein> proteinMap = new HashMap<>();
        for( Protein p: activeProteins ) {
            proteinMap.put(p.getUniprotId(), p);
        }
        log.info("   proteins in RGD: "+Utils.formatThousands(proteinMap.size()));

        for( String uniprotId: incomingCanonicalProteins ) {

            // is the protein in rgd?
            Protein p = proteinMap.get(uniprotId);
            if( p==null ) {
                counters.increment("canonical proteins not in RGD");
            } else if( p.isCanonical() ){
                counters.increment("canonical proteins already in RGD");
                proteinMap.remove(uniprotId);
            } else {
                boolean oldIsCanonical = p.isCanonical();
                p.setCanonical(true);
                dao.updateProteinCanonicalStatus(p, oldIsCanonical);
                counters.increment("proteins in RGD upgraded to canonical status");
                proteinMap.remove(uniprotId);
            }
        }

        // remaining proteins in RGD should not be canonical
        for( Protein p: proteinMap.values() ) {
            if( p.isCanonical() ) {
                boolean oldIsCanonical = p.isCanonical();
                p.setCanonical(false);
                dao.updateProteinCanonicalStatus(p, oldIsCanonical);
                counters.increment("proteins in RGD downgraded to non-canonical status");
            } else {
                counters.increment("proteins in RGD with retained non-canonical status");
            }
        }

        log.info(counters.dumpAlphabetically());
    }

    String downloadFile(String species, String fileUrl) throws Exception {

        FileDownloader fd = new FileDownloader();
        fd.setPrependDateStamp(true);
        fd.setExternalFile(fileUrl);
        fd.setLocalFile("data/"+species+"_canonical_proteins.fa");
        fd.setUseCompression(true);
        String localFile = fd.downloadNew();

        log.info("downloaded file "+localFile);
        return localFile;
    }

    List<String> getIncomingCanonicalProteins(String localFile) throws IOException {

        List<String> incomingCanonicalProteins = new ArrayList<>();

        BufferedReader in = Utils.openReader(localFile);
        String line;
        while( (line=in.readLine())!=null ) {

            // sample line to parse:
            // >sp|O75179|ANR17_HUMAN Ankyrin repeat domain-containing protein 17 OS=Homo sapiens OX=9606 GN=ANKRD17 PE=1 SV=3
            int barPos1 = line.indexOf('|');
            if( barPos1<0 ) {
                continue;
            }
            int barPos2 = line.indexOf('|', barPos1+1);
            if( barPos2<0 ) {
                continue;
            }
            String accId = line.substring(barPos1+1, barPos2);
            incomingCanonicalProteins.add(accId);
        }
        in.close();

        log.info("   incoming canonical proteins: "+Utils.formatThousands(incomingCanonicalProteins.size()));

        return incomingCanonicalProteins;
    }

    public Map<String, String> getCanonicalFiles() {
        return canonicalFiles;
    }

    public void setCanonicalFiles(Map<String, String> canonicalFiles) {
        this.canonicalFiles = canonicalFiles;
    }
}
