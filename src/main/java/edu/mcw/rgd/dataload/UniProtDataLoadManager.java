package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.PipelineLog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.log.RGDSpringLogger;
import edu.mcw.rgd.process.PipelineLogger;
import edu.mcw.rgd.process.Utils;
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
    PipelineLogger dbLogger = PipelineLogger.getInstance();
    UniProtDAO dao;
    UniProtFileParser fileParser;
    UniProtQC qc = new UniProtQC();
    UniProtDataValidation dataValidation;
    static long startMilisec=System.currentTimeMillis();
    static long runSec=0;
    RGDSpringLogger rgdLogger;

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
        System.out.println(dataManager.getVersion());

        System.out.println(dataManager.dao.getConnectionInfo());

        boolean loadRefSeq2UniprotMappings = false;

        try {
            String fileName = null, fileName2=null;
            boolean downloadOnly = false;
            int speciesTypeKey = SpeciesType.ALL;

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

        System.out.println("===== start "+SpeciesType.getCommonName(speciesTypeKey)+" UniProtKB pipeline =====");

        qc.setSpeciesTypeKey(speciesTypeKey);
        qc.setDao(dao);

        dao.setProcessingStartTime(new Date());
        Thread.sleep(555);

        // initialize pipeline logging to database for first file
        dbLogger.init(speciesTypeKey, "uniprot_sprot", PipelineLogger.PIPELINE_UNIPROT);

        dataValidation.setSpeciesTypeKey(speciesTypeKey);
        fileParser.processFile1(fileName, PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.SWISSPROT);

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

        dbLogger.close(true);

        // initialize pipeline logging to database for second file
        dbLogger.restart("uniprot_trembl");

        // reset database counter
        fileParser.getMapXdbCount().clear();

        fileParser.processFile2(fileName2, PipelineLogger.PIPELINE_UNIPROT+UniProtDAO.TREMBL);

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

        dbLogger.close(true);

        System.out.println("===== done "+SpeciesType.getCommonName(speciesTypeKey)+" UniProtKB pipeline =====");
        System.out.println();
    }

    // write summary to db pipeline logs
    public void writeDbSummary() throws Exception {
        dbLogger.log("Number of protein records in the source file", Integer.toString(fileParser.getTotRecord()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of protein records skipped (different species)", Integer.toString(fileParser.getSkipRecord()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of new xdb ids loaded into RGD", Integer.toString(dao.getRowsInserted()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of xdb ids deleted from RGD", Integer.toString(dao.getRowsDeleted()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of stale xdb ids deleted from RGD", Integer.toString(dao.getStaleRowsDeleted()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of matching xdb ids", Integer.toString(dao.getRowsMatched()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of inActiveGene", Integer.toString(qc.getInActiveGene()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of newActiveGene", Integer.toString(qc.getNewActiveGene()), PipelineLog.LOGPROP_TOTAL);
        dbLogger.log("Number of unMatched", Integer.toString(qc.getUnMatched()), PipelineLog.LOGPROP_TOTAL);

        dbLogger.log("Total Number of protein records in the source file: ", Integer.toString(fileParser.getTotRecord()), PipelineLog.LOGPROP_RECCOUNT);

        fileParser.dumpXdbCounts(this.dataValidation.getActiveXdbIdMap());

        long endMilisec=System.currentTimeMillis();
        runSec=endMilisec-startMilisec;
        dbLogger.log("Process runtime length: ", Utils.formatElapsedTime(startMilisec, endMilisec), PipelineLog.LOGPROP_TOTAL);
    }

    public void writeSummary() throws Exception {
        long endMilisec=System.currentTimeMillis();
        runSec=(endMilisec-startMilisec)/1000;
        System.out.println("Protein records found in the source file   : "+fileParser.getTotRecord());
        if( fileParser.getSkipRecord()>0 )
            System.out.println("   records skipped (different species): "+fileParser.getSkipRecord());
        if( dao.getRowsInserted()>0 )
            System.out.println("Xdb ids loaded into RGD  : "+dao.getRowsInserted());
        if( dao.getRowsDeleted()>0 )
            System.out.println("  xdb ids deleted        : "+dao.getRowsDeleted());
        if( dao.getStaleRowsDeleted()>0 )
            System.out.println("  xdb ids deleted (stale): "+dao.getStaleRowsDeleted());
        if( dao.getRowsMatched()>0 )
            System.out.println("  xdb ids matching RGD   : "+dao.getRowsMatched());

        if( dao.getAliasesUpToDate()>0 )
            System.out.println("  old_protein_id aliases up-to-date: "+dao.getAliasesUpToDate());
        if( dao.getAliasesInserted()>0 )
            System.out.println("  old_protein_id aliases inserted  : "+dao.getAliasesInserted());

        if( qc.getInActiveGene()>0 )
            System.out.println("inActiveGene : "+qc.getInActiveGene());
        if( qc.getNewActiveGene()>0 )
            System.out.println("newActiveGene: "+qc.getNewActiveGene());
        if( qc.getUnMatched()>0 )
            System.out.println("unMatched    : "+qc.getUnMatched());

        dumpGlobalCounters();

        System.out.println("Process completed in "+Utils.formatElapsedTime(startMilisec, endMilisec));

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

        System.out.println("=== global counters ===");
        for( Map.Entry<String,Integer> entry: counters.entrySet() ) {
            if( entry.getValue()!=0 ) {
                System.out.println("  "+entry.getKey()+":  "+entry.getValue());
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
