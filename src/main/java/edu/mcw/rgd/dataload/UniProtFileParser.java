package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.PipelineLog;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.*;

import java.io.*;
import java.util.*;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 */
public class UniProtFileParser {
    String ftpServer; // ftp server name
    String ftpDirectory; // ftp directory with source files
    String fileName;  // original sprot file
    String fileName2; // alternate trembl file
    int totRecord=0;   // total number of records found in the file
    int skipRecord=0;  // total number of records skipped (species not rat)
    int downloadMaxRetryCount; // how many times we should try to download the file
    int downloadRetryInterval; // interval in seconds between download retrials
    private int speciesTypeKey = SpeciesType.ALL;
    private int taxonid;

    private AccIdDumper accIdDumper;

    // map of counts of how many lines of particular database appears in the data
    Map<String, Integer> mapXdbCount = new HashMap<>(127);
    List<UniProtRatRecord> incomingRecords = new ArrayList<>(1000);

    UniProtDataValidation dataValidation;
    PipelineLogger dblog = PipelineLogger.getInstance();
    private String accIdDumpFile;
    private Map<String,String> swissProtFileNames;
    private Map<String,String> tremblFileNames;

    public void setSpecies(int speciesTypeKey) {

        String speciesCommonName = SpeciesType.getCommonName(speciesTypeKey);
        if( speciesCommonName!=null  ) {
            setFileName(getSwissProtFileNames().get(speciesCommonName));
            setFileName2(getTremblFileNames().get(speciesCommonName));
        } else {
            setFileName(null);
            setFileName2(null);
        }

        this.speciesTypeKey = speciesTypeKey;
        this.taxonid = SpeciesType.getTaxonomicId(speciesTypeKey);
    }

    public void processFile1(String fileName1, String srcPipeline) throws Exception {
        if( fileName1==null )
            fileName1 = download(this.fileName);

        // create a file for dumping primary and secondary SPROT acc ids
        accIdDumper = new AccIdDumper();
        accIdDumper.writeFileHeader(this.getAccIdDumpFile()+"_"+SpeciesType.getCommonName(this.speciesTypeKey));

        processFile(fileName1, srcPipeline);

        accIdDumper.close();
    }

    public void processFile2(String fileName2, String srcPipeline) throws Exception {
        if( fileName2==null )
            fileName2 = download(this.fileName2);

        // no file for dumping TREMBL acc ids
        accIdDumper = null;

        processFile(fileName2, srcPipeline);
    }

    /**
     * download the remote file and store it locally; current data should be appended to name of local file
     * @param fname file name to download
     * @return name of local file
     * @throws Exception if file download for any reason fails
     */
    public String download(String fname) throws Exception {

        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile("ftp://"+ftpServer+ftpDirectory+"/"+fname);
        downloader.setLocalFile("data/"+fname);
        downloader.setAppendDateStamp(true);
        return downloader.downloadNew();
    }

    public void processFile(String fileName, String srcPipeline) throws Exception {
        dblog.log("starting processing of ", fileName, PipelineLogger.INFO);

        // open input stream
        BufferedReader reader = Utils.openReader(fileName);

        String line;
        String uniprotID = null;
        String uniprotAC = null;
        UniProtRatRecord rec = new UniProtRatRecord(srcPipeline);
        int recNo = 0;
        boolean speciesIsCorrect = true;
        String lineOX = ""; // can span multiple lines

        while ((line=reader.readLine()) != null) {

            // handle new record of lines
            if (line.startsWith("//")) {
                //increment counter
                totRecord++;

                // add the UniProt entry
                if( uniprotID!=null && uniprotAC!=null ) {
                    rec.addEntry("UniProt", uniprotAC, uniprotID, null);
                    rec.setUniProtAccId(uniprotAC);
                }

                // keep track of RECNO
                recNo++;
                rec.setRecNo(recNo);

                // check organism
                if( speciesIsCorrect ) {
                    // fix domain names
                    for( ProteinDomain domain: rec.domains ) {
                        domain.cleanupDomainName();
                    }
                    incomingRecords.add(rec);
                } else
                    skipRecord++; // count records of species different from rat

                dblog.writeLogProps();
                dblog.removeAllLogProps();

                // prepare for next record
                rec = new UniProtRatRecord(srcPipeline);
                speciesIsCorrect = false;
                uniprotID = uniprotAC = null;
                lineOX = "";
                continue;
            }

            // process species
            if( line.startsWith("OX") ) {
                lineOX += line; // OX can span multiple lines ...
                if( parseTaxonId(lineOX) == this.taxonid ) {
                    speciesIsCorrect = true;
                }
                continue;
            }
            // extract uniprot id
            if( line.startsWith("ID   ") ) {
                uniprotID = extractWord(line, 5);
                continue;
            }
            // extract uniprot accession id
            if( line.startsWith("AC   ") ) {
                // for some proteins there are multiple AC lines
                // the primary accession id is the first accId from the first line
                if( uniprotAC==null ) {
                    uniprotAC = extractWord(line, 5);

                    // AC   Q9CQV8; O70455; Q3TY33; Q3UAN6;
                    //      primary|secondary ids
                    String lineWithSecondaryAccIds = line.substring(5+uniprotAC.length()+1).trim();
                    rec.parseSecondaryAccessionIds(lineWithSecondaryAccIds);
                    if( accIdDumper!=null ) {
                        accIdDumper.writeData(uniprotAC, lineWithSecondaryAccIds);
                    }
                } else {
                    rec.parseSecondaryAccessionIds(line.substring(5).trim());
                }
                continue;
            }
            // extract protein name
            if( line.startsWith("DE   ") ) {
                parseProteinName(line, rec);
                continue;
            }
            // extract protein feature
            if( line.startsWith("FT   ") ) {
                parseProteinFeature(line, rec);
                continue;
            }
            // extract protein sequence
            if( line.startsWith("     ") ) {
                parseProteinSequence(line, rec);
                continue;
            }

            // skip all other lines not beginning with "DR"
            if( !line.startsWith("DR   ") )
                continue;
            if( line.endsWith(".") )
                line = line.substring(5, line.length()-1); // skip the terminating '.' character
            else
                line = line.substring(5); // skip beginning of the line : "DR   "

            // split the line by "; "
            String[] fields = line.split("; ");

            // first line always contains database name string
            if( fields.length<2 )
                continue; // invalid line
            String xdbName = fields[0];
            String accId = fields[1];
            String link1 = fields.length>2 ? fields[2] : null;
            String link2 = fields.length>3 ? fields[3] : null;

            // customized pre-processing for some of entries
            switch (xdbName) {
                case "Ensembl":
                    // in one line there are several ids that in RGD fall into separate XDB_KEYS
                    // the source line looks like this: ENSRNOT00000016981; ENSRNOP00000016981; ENSRNOG00000010945
                    // or this:               ENSMUST00000018470; ENSMUSP00000018470; ENSMUSG00000018326. [Q9CQV8-1]
                    rec.addEnsemblEntry(accId);
                    rec.addEnsemblEntry(link1);
                    rec.addEnsemblEntry(extractWord(link2, 0));
                    break;
                case "RefSeq":
                    parseRefSeq(xdbName, accId, link1, rec);
                    break;
                case "SUPFAM":
                    // convert 'SSF50353' ACC_ID into '50353'
                    if (accId != null && accId.length() > 3)
                        rec.addEntry(xdbName, accId.substring(3), link1, link2);
                    break;
                case "GermOnline":
                    // source data puts species name, like 'Mus musculus' in link1 column; that's useless
                    // much better to use acc id as link text
                    rec.addEntry(xdbName, accId, accId, null);
                    break;
                default:
                    rec.addEntry(xdbName, accId, link1, link2);
                    break;
            }

            // update count for every database name
            Integer count = mapXdbCount.get(xdbName);
            if( count==null )
                count = 0;
            count++;
            mapXdbCount.put(xdbName, count);
        }
        reader.close();

        dblog.log("finished processing of ", fileName, PipelineLogger.INFO);
        System.out.println("finished processing of "+ fileName);
    }

    /// OX   NCBI_TaxID=442598 {ECO:0000313|EMBL:AXQ06062.1};
    int parseTaxonId(String lineOX) {

        String pattern = "NCBI_TaxID=";
        int posStart = lineOX.indexOf(pattern);
        if( posStart<0 ) {
            return 0;
        }
        posStart += pattern.length();

        // continue as long as there are digits
        int taxonId = 0;
        for( int pos=posStart; pos<lineOX.length(); pos++ ) {
            char c = lineOX.charAt(pos);
            if( !Character.isDigit(c) ) {
                break;
            }
            taxonId *= 10;
            taxonId += Character.digit(c, 10);
        }
        return taxonId;
    }

    void parseRefSeq(String xdbName, String accId, String nucId, UniProtRatRecord ratData) {
        // DR   RefSeq; NP_001258682.1; NM_001271753.1. [O88508-1]

        // strip refseq versions from accession:
        // 'NP_001019465.2','NM_001024294.1' => 'NP_001019465','NM_001024294'
        if( accId!=null ) {
            int dotPos = accId.indexOf('.');
            if( dotPos>0 ) // 'NP_001019465.2' => 'NP_001019465'
                accId = accId.substring(0, dotPos);
        }

        // do not use NM_xxx id as LINK_TEXT; use NP_xxx id
        ratData.addEntry(xdbName, accId, accId, null);
    }

    // analyze the given string starting from 'startingPos' and extract the longest string
    // composed entirely of letters, digits and '_'
    String extractWord(String str, int startingPos) {
        StringBuilder buf = new StringBuilder(str.length());
        for( int i=startingPos; i<str.length(); i++ ) {
            char c = str.charAt(i);
            if( Character.isJavaIdentifierPart(c) )
                buf.append(c);
            else
                break;
        }
        return buf.toString();
    }

    void parseProteinName(String line, UniProtRatRecord rec) {
        final String RECNAME = "RecName: Full=";
        final String SUBNAME = "SubName: Full=";

        int recNamePos = line.indexOf(RECNAME);
        if( recNamePos>0 ) {
            recNamePos += RECNAME.length();
            int recNameEnd = line.indexOf(";", recNamePos);
            if( recNameEnd>recNamePos ) {
                String proteinName = cleanUpProteinName(line.substring(recNamePos, recNameEnd));
                rec.setProteinName(proteinName);
            } else {
                System.out.println(" problems with recommended name");
            }
        } else {
            // parse submitted name only if recommended name is not present
            if( Utils.isStringEmpty(rec.getProteinName()) ) {
                int subNamePos = line.indexOf(SUBNAME);
                if (subNamePos > 0) {
                    subNamePos += SUBNAME.length();
                    int subNameEnd = line.indexOf(";", subNamePos);
                    if (subNameEnd > subNamePos) {
                        String proteinName = cleanUpProteinName(line.substring(subNamePos, subNameEnd));
                        rec.setProteinName(proteinName);
                    } else {
                        System.out.println(" problems with submitted name");
                    }
                }
            }
        }
    }

    void parseProteinSequence(String line, UniProtRatRecord rec) {
        if( line.startsWith("     ") ) {
            String seqFragment = line.replaceAll("\\s+", "");
            if( rec.proteinSequence==null ) {
                rec.proteinSequence = seqFragment;
            } else {
                rec.proteinSequence += seqFragment;
            }
        }
    }

    void parseProteinFeature(String line, UniProtRatRecord rec) {
        if( line.startsWith("FT   DOMAIN") ) {
            inProteinDomain = true;
            // new protein domain record
            ProteinDomain domain = new ProteinDomain();
            domain.aaStartPos = parseFeaturePos(line.substring(14, 20));
            domain.aaStopPos = parseFeaturePos(line.substring(20, 27));
            String text = line.substring(34);
            int dotPos = text.indexOf('.');
            if( dotPos>0 ) {
                domain.setDomainName(text.substring(0, dotPos));
                domain.notes = text.substring(dotPos+1);
            } else {
                domain.setDomainName(text);
            }

            rec.domains.add(domain);
        } else if( !line.startsWith("FT         ") ) {
            // feature different than domain
            inProteinDomain = false;
        } else { // feature continuation
            if( inProteinDomain ) {
                ProteinDomain domain = rec.domains.get(rec.domains.size()-1);
                String text = line.substring(34);
                // is the domain name complete?
                if( domain.getDomainName().endsWith(".") ) {
                    // domain name is complete -- append notes
                    domain.notes += text;
                } else {
                    // domain name is not complete
                    int dotPos = text.indexOf('.');
                    if( dotPos>0 ) {
                        domain.setDomainName(domain.getDomainName() + " "+text.substring(0, dotPos));
                        domain.notes = text.substring(dotPos+1);
                    } else {
                        domain.setDomainName(domain.getDomainName() + " "+text);
                    }
                }
            }
        }
    }
    boolean inProteinDomain;

    int parseFeaturePos(String str) {
        String pos = str.trim();
        if( pos.startsWith("<") || pos.startsWith(">") || pos.startsWith("?") ) {
            pos = pos.substring(1);
        }
        return Integer.parseInt(pos);
    }

    // remove any contents in braces
    String cleanUpProteinName(String proteinName) {
        int braceStart = proteinName.indexOf('{');
        int braceStop = proteinName.indexOf('}');
        if( braceStart>=0 && braceStop>braceStart ) {
            return (proteinName.substring(0, braceStart) + proteinName.substring(braceStop+1)).trim();
        } else {
            return proteinName;
        }
    }

    public void setAccIdDumpFile(String accIdDumpFile) {
        this.accIdDumpFile = accIdDumpFile;
    }

    public String getAccIdDumpFile() {
        return accIdDumpFile;
    }

    public void setSwissProtFileNames(Map<String,String> swissProtFileNames) {
        this.swissProtFileNames = swissProtFileNames;
    }

    public Map<String,String> getSwissProtFileNames() {
        return swissProtFileNames;
    }

    public void setTremblFileNames(Map<String,String> tremblFileNames) {
        this.tremblFileNames = tremblFileNames;
    }

    public Map<String,String> getTremblFileNames() {
        return tremblFileNames;
    }

    // helper class to hold counts of external database ids processed by program
    class XdbIdCount implements Comparable<XdbIdCount> {
        public String xdbName; // name of external database
        public int count; // count of lines with the external database name

        public XdbIdCount(String xdbName, int count) {
            this.xdbName = xdbName;
            this.count = count;
        }

        public int compareTo(XdbIdCount o) {
            // the highest counts go first
            int r = o.count - this.count;
            if( r != 0 )
                return r;
            // if counts are equal, match by xdb name
            return this.xdbName.compareTo(o.xdbName);
        }
    }

    public void dumpXdbCounts(Map<String, Integer> activeXdbIdMap) throws Exception {
        // create a sorted list of XdbIdCount objects
        List<XdbIdCount> slist = new ArrayList<>(this.mapXdbCount.size());
        for( Map.Entry<String, Integer> entry: mapXdbCount.entrySet() ) {
            slist.add(new XdbIdCount(entry.getKey(), entry.getValue()));
        }
        Collections.sort(slist);

        // dump the counts
        int inactiveDbCount = 0;
        int itemCountForInactiveDatabases = 0;
        for( XdbIdCount item: slist ) {
            if( activeXdbIdMap.containsKey(item.xdbName))
                dblog.log("# lines for [ACTIVE] database "+item.xdbName, Integer.toString(item.count), PipelineLog.LOGPROP_TOTAL);
            else {
                inactiveDbCount++;
                itemCountForInactiveDatabases += item.count;
            }
        }
        dblog.log("# lines for [IGNORED] "+inactiveDbCount+" databases", Integer.toString(itemCountForInactiveDatabases), PipelineLog.LOGPROP_TOTAL);
    }



    /**
     * @return Returns the fileName.
     */
    public String getFileName() {
        return fileName;
    }
    /**
     * @param fileName The fileName to set.
     */
    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    /**
     * @return Returns the totRecord.
     */
    public int getTotRecord() {
        return totRecord;
    }
    /**
     * @param totRecord The totRecord to set.
     */
    public void setTotRecord(int totRecord) {
        this.totRecord = totRecord;
    }

    public String getFileName2() {
        return fileName2;
    }

    public void setFileName2(String fileName2) {
        this.fileName2 = fileName2;
    }

    public int getSkipRecord() {
        return skipRecord;
    }

    public void setSkipRecord(int skipRecord) {
        this.skipRecord = skipRecord;
    }

    public Map<String, Integer> getMapXdbCount() {
        return mapXdbCount;
    }

    public void setMapXdbCount(Map<String, Integer> mapXdbCount) {
        this.mapXdbCount = mapXdbCount;
    }

    public String getFtpServer() {
        return ftpServer;
    }

    public void setFtpServer(String ftpServer) {
        this.ftpServer = ftpServer;
    }

    public String getFtpDirectory() {
        return ftpDirectory;
    }

    public void setFtpDirectory(String ftpDirectory) {
        this.ftpDirectory = ftpDirectory;
    }

    public int getDownloadMaxRetryCount() {
        return downloadMaxRetryCount;
    }

    public void setDownloadMaxRetryCount(int downloadMaxRetryCount) {
        this.downloadMaxRetryCount = downloadMaxRetryCount;
    }

    public int getDownloadRetryInterval() {
        return downloadRetryInterval;
    }

    public void setDownloadRetryInterval(int downloadRetryInterval) {
        this.downloadRetryInterval = downloadRetryInterval;
    }

    public UniProtDataValidation getDataValidation() {
        return dataValidation;
    }

    public void setDataValidation(UniProtDataValidation dataValidation) {
        this.dataValidation = dataValidation;
    }

    public List<UniProtRatRecord> getIncomingRecords() {
        return incomingRecords;
    }
}
