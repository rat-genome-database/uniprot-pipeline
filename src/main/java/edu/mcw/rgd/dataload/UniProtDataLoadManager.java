package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.log.RGDSpringLogger;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.Date;
import java.util.TreeMap;
import java.util.List;
import java.util.Map;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 */
public class UniProtDataLoadManager {
    UniProtDAO dao;
    UniProtFileParser fileParser;
    UniProtQC qc = new UniProtQC();
    UniProtDataValidation dataValidation;
    static long startMilisec=System.currentTimeMillis();
    static long runSec=0;
    RGDSpringLogger rgdLogger;
    Logger logMain = Logger.getLogger("main");

    int speciesTypeKey;
    private String version;
    private Map<String,Integer> counters = new TreeMap<>();

    static private UniProtDataLoadManager manager;
    private ProteinLoader proteinLoader;

    protected UniProtDataLoadManager() {
    }

    static public UniProtDataLoadManager getInstance() {
        return manager;
    }

    /**
     * Download SPROT and TREMBL files from FTP site of UniProtKB consortium, parse them and load into RGD database.
     */
    public static void main(String[] args) throws Exception {

        if( args.length <= 0) {
            usage();
            return;
        }

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        UniProtDataLoadManager dataManager = (UniProtDataLoadManager) (bf.getBean("dataLoadManager"));
        manager = dataManager;

        // ensure this text matches the 'appver' property from build.xml file
        manager.logMain.info(dataManager.getVersion());

        manager.logMain.info("   "+dataManager.dao.getConnectionInfo());

        boolean loadRefSeq2UniprotMappings = false;

        try {
            String fileName = null, fileName2=null;
            boolean downloadOnly = false;
            int speciesTypeKey = SpeciesType.ALL;
            boolean loadProteinDomains = false;

            for( int i=0; i<args.length; i++ ) {
                // optional: only download source files
                switch (args[i]) {
                    case "-download":
                        downloadOnly = true;
                        break;
                    // optional full path to both files
                    case "-file":
                        if (i + 2 >= args.length) {
                            usage();
                            return;
                        }
                        // override source data files
                        fileName = args[++i];
                        fileName2 = args[++i];
                        break;
                    case "-species":
                        if (i + 1 >= args.length) {
                            usage();
                            return;
                        }
                        // try to parse species
                        speciesTypeKey = SpeciesType.parse(args[++i]);
                        break;
                    case "-skipProteinLoader":
                        dataManager.proteinLoader = null;
                        break;
                    case "--loadRefSeq2UniProt":
                        loadRefSeq2UniprotMappings = true;
                        break;

                    case "--loadProteinDomains":
                        loadProteinDomains = true;
                        break;

                    case "--deletedAccessions":
                        DeletedAccessions module = (DeletedAccessions) (bf.getBean("deletedAccessions"));
                        module.run();
                        return;
                }
            }

            if( loadRefSeq2UniprotMappings ) {
                Refseq2UniprotLoader loader = (Refseq2UniprotLoader) (bf.getBean("RefSeq2UniProtLoader"));
                loader.run();
                return;
            }

            // species must be set
            if( speciesTypeKey==SpeciesType.ALL ) {
                usage();
                return;
            }
            dataManager.speciesTypeKey = speciesTypeKey;
            dataManager.fileParser.setSpecies(speciesTypeKey);

            // optional:
            if( downloadOnly ) {
                dataManager.download();
                return;
            }

            if( loadProteinDomains ) {
                ProteinDomainLoader module = (ProteinDomainLoader) (bf.getBean("proteinDomainLoader"));
                module.run(dataManager.fileParser);
                return;
            }
            dataManager.startPipeline(fileName, fileName2, speciesTypeKey);
        }
        catch( Exception e ) {
            e.printStackTrace();
            throw e;
        }
    }

    static private void usage() {

        System.out.println(
                "UniProtKB pipeline command line arguments:\n" +
                        "-species rat|mouse|human|1|2|3|...  -- [mandatory] specify species to process\n" +
                        "[-download]  -- [optional] download source files and exit\n" +
                        "[-file SwissProtFile TremblFile]  -- [optional] override input files\n"
        );
    }

    public void download() throws Exception {
        // download only
        fileParser.download(fileParser.getFileName());
        fileParser.download(fileParser.getFileName2());
    }

    /**
     * start the pipeline - both in uniprot_sprot and uniprot_trembl mode
     * @param fileName optional fileName pointing to file containing the sprot data to be analyzed
     * @param fileName2 optional fileName2 pointing to another file containing the trembl data to be analyzed
     * @param speciesTypeKey species type key
     * @throws Exception when unexpected error in spring framework occurs
     */
    public void startPipeline(String fileName, String fileName2, int speciesTypeKey) throws Exception {

        logMain.info("===== start "+SpeciesType.getCommonName(speciesTypeKey)+" UniProtKB pipeline =====");

        qc.setSpeciesTypeKey(speciesTypeKey);
        qc.setDao(dao);

        dao.setProcessingStartTime(new Date());
        Thread.sleep(555);

        dataValidation.setSpeciesTypeKey(speciesTypeKey);
        fileParser.processFile1(fileName, UniProtDAO.SWISSPROT);

        // take snapshot of the statistics after analyzing the first file
        int relationshipsIncoming = fileParser.getTotRecord();
        int relationshipsSkipped = fileParser.getSkipRecord();
        int relationshipsUpdated = dao.getRowsInserted();
        int genesInactive = qc.getInActiveGene();
        int genesActive = qc.getNewActiveGene();
        int externalRefDuplicate = dao.getRowsMatched();
        int externalRefDeleted = dao.getRowsDeleted();
        int externalRefUnmatched = qc.getUnMatched();
        // print statistics to database
        writeDbSummary();

        // reset database counter
        fileParser.getMapXdbCount().clear();

        fileParser.processFile2(fileName2, UniProtDAO.TREMBL);

        // QC
        List<UniProtRatRecord> incomingRecords = fileParser.getIncomingRecords();
        qc.qc(incomingRecords);
        qc.dumpMatchSummary();

        // RUN protein loader
        if( proteinLoader!=null ) {
            proteinLoader.setDao(dao);
            proteinLoader.setSpeciesTypeKey(speciesTypeKey);
            proteinLoader.run(incomingRecords);
            setProteinLoader(null); // to release memory
        }

        // DATA LOAD
        dataValidation.mergeIncomingRecords(incomingRecords);
        dataValidation.load();

        dao.deleteStaleRows(speciesTypeKey);

        writeSummary();

        // update the stats: remove from totals the statistics from 1st file
        // so the totals will contain numbers from second file only
        fileParser.setTotRecord(-relationshipsIncoming + fileParser.getTotRecord());
        fileParser.setSkipRecord(-relationshipsSkipped + fileParser.getSkipRecord());
        dao.setRowsInserted(-relationshipsUpdated + dao.getRowsInserted());
        qc.setInActiveGene(-genesInactive + qc.getInActiveGene());
        qc.setNewActiveGene(-genesActive + qc.getNewActiveGene());
        dao.setRowsMatched(-externalRefDuplicate + dao.getRowsMatched());
        dao.setRowsDeleted(-externalRefDeleted + dao.getRowsDeleted());
        qc.setUnMatched(-externalRefUnmatched + qc.getUnMatched());
        writeDbSummary();

        logMain.info("===== done "+SpeciesType.getCommonName(speciesTypeKey)+" UniProtKB pipeline =====");
        logMain.info("");
    }

    // write summary to db pipeline logs
    public void writeDbSummary() throws Exception {
        fileParser.dumpXdbCounts(dataValidation.getActiveXdbIdMap(), logMain);

        long endMilisec=System.currentTimeMillis();
        runSec=endMilisec-startMilisec;
    }

    public void writeSummary() throws Exception {
        long endMilisec=System.currentTimeMillis();
        runSec=(endMilisec-startMilisec)/1000;
        logMain.info("Protein records found in the source file   : "+fileParser.getTotRecord());
        if( fileParser.getSkipRecord()>0 )
            logMain.info("   records skipped (different species): "+fileParser.getSkipRecord());
        if( dao.getRowsInserted()>0 )
            logMain.info("Xdb ids loaded into RGD  : "+dao.getRowsInserted());
        if( dao.getRowsDeleted()>0 )
            logMain.info("  xdb ids deleted        : "+dao.getRowsDeleted());
        if( dao.getStaleRowsDeleted()>0 )
            logMain.info("  xdb ids deleted (stale): "+dao.getStaleRowsDeleted());
        if( dao.getRowsMatched()>0 )
            logMain.info("  xdb ids matching RGD   : "+dao.getRowsMatched());

        if( dao.getAliasesUpToDate()>0 )
            logMain.info("  old_protein_id aliases up-to-date: "+dao.getAliasesUpToDate());
        if( dao.getAliasesInserted()>0 )
            logMain.info("  old_protein_id aliases inserted  : "+dao.getAliasesInserted());

        if( qc.getInActiveGene()>0 )
            logMain.info("inActiveGene : "+qc.getInActiveGene());
        if( qc.getNewActiveGene()>0 )
            logMain.info("newActiveGene: "+qc.getNewActiveGene());
        if( qc.getUnMatched()>0 )
            logMain.info("unMatched    : "+qc.getUnMatched());

        dumpGlobalCounters();

        logMain.info("Process completed in "+Utils.formatElapsedTime(startMilisec, endMilisec));

        // record logs in database
        String subSystem = "UniProtKB"+SpeciesType.getCommonName(this.speciesTypeKey);
        rgdLogger.log(subSystem,"proteinRecordsTotal",fileParser.getTotRecord());
        rgdLogger.log(subSystem,"proteinRecordsSkipped",fileParser.getSkipRecord());
        rgdLogger.log(subSystem,"relationsDeleted",dao.getRowsDeleted());
        rgdLogger.log(subSystem,"relationsUpdated",dao.getRowsInserted());
        rgdLogger.log(subSystem,"genesInactive",qc.getInActiveGene());
        rgdLogger.log(subSystem,"genesActive",qc.getNewActiveGene());
        rgdLogger.log(subSystem,"relationsDuplicate",dao.getRowsMatched());
        rgdLogger.log(subSystem,"relationsUnmatched",qc.getUnMatched());
        rgdLogger.log(subSystem,"timeToRun",runSec);

    }

    void dumpGlobalCounters() {

        logMain.info("=== global counters ===");
        for( Map.Entry<String,Integer> entry: counters.entrySet() ) {
            if( entry.getValue()!=0 ) {
                logMain.info("  "+entry.getKey()+":  "+entry.getValue());
            }
        }
    }

    public int incrementCounter(String counterName, int delta) {

        Integer count = counters.get(counterName);
        if( count==null )
            count = delta;
        else
            count += delta;
        counters.put(counterName, count);
        return count;
    }

    public int incrementCounter(String counterName) {
        return incrementCounter(counterName, 1);
    }

    public RGDSpringLogger getRgdLogger() {
        return rgdLogger;
    }
    public void setRgdLogger(RGDSpringLogger rgdLogger) {
        this.rgdLogger = rgdLogger;
    }

    public UniProtDAO getUniProtDAO() {
        return dao;
    }

    public void setUniProtDAO(UniProtDAO uniProtDAO) {
        this.dao = uniProtDAO;
    }

    public UniProtFileParser getFileParser() {
        return fileParser;
    }

    public void setFileParser(UniProtFileParser fileParser) {
        this.fileParser = fileParser;
    }

    public UniProtDataValidation getDataValidation() {
        return dataValidation;
    }

    public void setDataValidation(UniProtDataValidation dataValidation) {
        this.dataValidation = dataValidation;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setProteinLoader(ProteinLoader proteinLoader) {
        this.proteinLoader = proteinLoader;
    }

    public ProteinLoader getProteinLoader() {
        return proteinLoader;
    }
}
