package edu.mcw.rgd.dataload;

import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * load protein domain objects
 */
public class ProteinDomainLoader {

    private List<Integer> processedMapKeys;
    private String srcPipeline;

    private UniProtDAO dao = new UniProtDAO();

    TranscriptDAO tdao = new TranscriptDAO();

    private edu.mcw.rgd.datamodel.Map assembly; // current assembly being processed
    private int taxonid; // current taxonid

    private Logger log = Logger.getLogger("domains");
    private Logger logStrandProblem = Logger.getLogger("strand_problem");

    private java.util.Map<Integer, List<ProteinDomain>> domainsMap = new HashMap<>();
    CounterPool counters;

    public void run(UniProtFileParser fileParser) throws Exception {
        List<Integer> mapKeys = new ArrayList<>(getProcessedMapKeys());
        Collections.shuffle(mapKeys);

        for( int mapKey: mapKeys ) {
            this.assembly = MapManager.getInstance().getMap(mapKey);
            this.taxonid = SpeciesType.getTaxonomicId(assembly.getSpeciesTypeKey());
            this.domainsMap.clear();
            this.counters = new CounterPool();
            fileParser.setSpecies(assembly.getSpeciesTypeKey());

            run(fileParser.download(fileParser.getFileName()), fileParser.download(fileParser.getFileName2()));
        }
    }

    public void run(String fileName1, String fileName2) throws Exception {

        String speciesName = SpeciesType.getCommonName(assembly.getSpeciesTypeKey());
        log.info("START for "+speciesName +" "+assembly.getName());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyyMMdd");
        String datePrefix = sdt.format(new Date());


        // parse sprot and trembl files and extract lines with primary accession ids and protein domain info
        // store this as a lean info file
        String outFileName = "data/"+datePrefix+"_"+speciesName.toLowerCase()+"_domains.txt";
        BufferedWriter out = new BufferedWriter(new FileWriter(outFileName));
        processFile(fileName1, out);
        processFile(fileName2, out);
        out.write("//\n");
        out.close();

        // parse the lean source file and create ProteinDomain objects
        List<ProteinDomain> domains = parse(outFileName);

        qc(domains);
        load(domains);

        log.info(counters.dumpAlphabetically());

        log.info("END for "+speciesName +" "+assembly.getName());
        log.info("");
    }

    void processFile(String fileName, BufferedWriter out) throws Exception {
        log.info(" processing "+fileName);

        // open input stream
        BufferedReader reader = Utils.openReader(fileName);

        String line;
        boolean speciesIsCorrect = false;
        String lineOX = ""; // can span multiple lines
        List<String> lines = new ArrayList<>();
        String prevLine = null;

        while ((line=reader.readLine()) != null) {

            // handle new record
            if (line.startsWith("//")) {

                // flush all valid lines to the file
                if( speciesIsCorrect && lines.size()>1 ) {
                    out.write("//\n");
                    int acCount = 0;
                    for (String l : lines) {
                        if( l.startsWith("AC") ) {
                            acCount++;
                        }
                        // for some proteins there are multiple AC lines
                        // the primary accession id is the first accId from the first line
                        // so we skip 2nd and subsequent AC lines
                        if( acCount>1 && l.startsWith("AC") ) {
                            continue;
                        }

                        out.write(l);
                        out.write("\n");
                    }
                }

                // init for new protein
                speciesIsCorrect = false;
                lineOX = "";
                lines.clear();

                prevLine = line;
                continue;
            }

            // process species
            if (line.startsWith("OX")) {
                lineOX += line; // OX can span multiple lines ...
                if (UniProtFileParser.parseTaxonId(lineOX) == this.taxonid) {
                    speciesIsCorrect = true;
                }
            } else if (line.startsWith("AC   ")) {
                lines.add(line);
            } else if (line.startsWith("FT   DOMAIN")) {
                lines.add(line);
            } else if (line.startsWith("FT         ")) {
                if (prevLine.startsWith("FT   DOMAIN")) {
                    // strip 'FT    '
                    // -- find 1st non-space character
                    int i;
                    for( i=2; i<line.length(); i++ ) {
                        if( line.charAt(i)!=' ') {
                            break;
                        }
                    }
                    String oldLine = lines.get(lines.size()-1);
                    String newLine = oldLine + line.substring(i);
                    lines.set(lines.size()-1, newLine);
                }
            }
            prevLine = line;
        }
        reader.close();
    }

    List<ProteinDomain> parse(String ftFileName) throws Exception {
        List<ProteinDomain> domains = new ArrayList<>();

        BufferedReader in = Utils.openReader(ftFileName);

        String line;
        String acc = null;
        while( (line=in.readLine())!=null ) {

            // new record
            if( line.startsWith("//") ) {
                acc = null;
            }

            // extract uniprot accession id
            if( line.startsWith("AC   ") ) {
                acc = UniProtFileParser.extractWord(line, 5);
                continue;
            }

            // extract protein feature
            ProteinDomain pd = parseProteinFeature(line, acc);
            if( pd!=null ) {
                domains.add(pd);
            }
        }
        in.close();

        log.info("PROTEIN DOMAIN RECORDS PARSED: "+domains.size());
        return domains;
    }

    ProteinDomain parseProteinFeature(String line, String acc) throws Exception {
        ProteinDomain domain = null;

        if( line.startsWith("FT   DOMAIN") ) {
            // new protein domain record
            domain = new ProteinDomain();
            domain.uniprotAcc = acc;

            parseDomainPos(line, domain);

            // Example line:
            // FT   DOMAIN          525..758/note="ABC transporter 1"
            int startPos = line.indexOf("/note=\"") + 7;
            int endPos = line.indexOf("\"", startPos);
            if( endPos < 0 ) { // note: there are some lines without terminating double quotes!!! we must handle them
                endPos = line.length();
            }
            if( endPos > startPos ) {
                domain.setDomainName(line.substring(startPos, endPos));
            } else {
                throw new Exception("domain parsing problem: "+line);
            }

            domain.qcDomainName();
        }
        return domain;
    }

    // the position is given like that:
    // FT   DOMAIN          231..303/note="RRM 3"
    boolean parseDomainPos(String line, ProteinDomain pd) {
        int startPos = 21;
        int doubleDotPos = line.indexOf("..", startPos);
        int slashPos = line.indexOf("/", doubleDotPos);
        if( doubleDotPos<0 || slashPos<0 ) {
            return false;
        }
        String pos1 = line.substring(startPos, doubleDotPos);
        pd.aaStartPos = parseFeaturePos(pos1);

        String pos2 = line.substring(doubleDotPos+2, slashPos);
        pd.aaStopPos = parseFeaturePos(pos2);
        return true;
    }

    int parseFeaturePos(String str) {
        String pos = str.trim();
        if( pos.startsWith("<") || pos.startsWith(">") || pos.startsWith("?") ) {
            pos = pos.substring(1);
        }
        return Integer.parseInt(pos);
    }

    void qc(List<ProteinDomain> domains) throws Exception {

        Set<String> domainNames = new HashSet<>();

        for( ProteinDomain pd: domains ) {
            pd.geInRgd = dao.getProteinDomainObject(pd.getDomainName());

            if( pd.geInRgd!=null ) {
                pd.loci = positionProteinDomain(pd, pd.uniprotAcc);
            }

            domainNames.add(pd.getDomainName());
        }

        log.info("PROTEIN DOMAIN COUNT: "+domainNames.size());
    }

    List<MapData> positionProteinDomain(ProteinDomain pd, String uniProtAccId) throws Exception {

        List<MapData> results = new ArrayList<>();

        // convert 1-based aaPos into 0-based nucPos
        int nucStartPos = (pd.aaStartPos-1)*3; // including
        int lenToGo = (pd.aaStopPos-pd.aaStartPos+1)*3;

        // determine primary assembly
        CdsUtils utils = new CdsUtils(tdao, assembly.getKey());

        // map protein domain to protein
        Protein p = dao.getProteinByUniProtId(uniProtAccId);
        if( p==null ) {
            log.warn("*** ERROR: unexpected: no protein for "+uniProtAccId);
            return results;
        }

        // map protein to ncbi protein acc ids
        List<XdbId> ncbiProtAccIds = dao.getXdbIdsByRgdId(XdbId.XDB_KEY_GENEBANKPROT, p.getRgdId());
        // get NCBI transcript given ncbi protein acc ids
        for( XdbId xid: ncbiProtAccIds ) {
            List<Transcript> transcripts = tdao.getTranscriptsByProteinAccId(xid.getAccId());
            for (Transcript tr : transcripts) {
                for (MapData md : tr.getGenomicPositions()) {
                    if( md.getMapKey()!=assembly.getKey() ) {
                        continue;
                    }
                    List<CodingFeature> cfs = utils.buildCfList(md);

                    List<MapData> mds;
                    if( md.getStrand().equals("-") ) {
                        mds = handleNegativeStrand(cfs, nucStartPos, lenToGo, pd.geInRgd.getRgdId(), uniProtAccId);
                    } else {
                        mds = handlePositiveStrand(cfs, nucStartPos, lenToGo, pd.geInRgd.getRgdId(), uniProtAccId);
                    }

                    // merge all cds chunks of protein domain position into one range
                    MapData mdRange = null;
                    for( MapData m: mds ) {
                        if( mdRange==null ) {
                            mdRange = m.clone();
                        } else {
                            if( m.getStartPos() < mdRange.getStartPos() ) {
                                mdRange.setStartPos( m.getStartPos() );
                            }
                            if( m.getStopPos() > mdRange.getStopPos() ) {
                                mdRange.setStopPos( m.getStopPos() );
                            }
                        }
                    }
                    if( mdRange!=null ) {
                        results.add(mdRange);
                    }
                }
            }
        }
        return results;
    }

    List<MapData> handlePositiveStrand(List<CodingFeature> cfs, int nucStartPos, int lenToGo, int domainRgdId, String uniProtAccId) {

        List<MapData> results = new ArrayList<>();

        // iterate over CDS features
        for (CodingFeature cf : cfs) {
            if (cf.getFeatureType() == TranscriptFeature.FeatureType.CDS) {
                // we found a CDS
                int cdsLen = cf.getStopPos() - cf.getStartPos() + 1;
                if (nucStartPos >= cdsLen) {
                    // nucStartPos is outside of this CDS
                    nucStartPos -= cdsLen;
                    continue;
                }

                // nucStartPos is within this cds
                //
                MapData mdDomain = new MapData();
                mdDomain.setStartPos(cf.getStartPos() + nucStartPos);
                mdDomain.setChromosome(cf.getChromosome());
                mdDomain.setMapKey(cf.getMapKey());
                mdDomain.setSrcPipeline(getSrcPipeline());
                mdDomain.setStrand(cf.getStrand());
                mdDomain.setRgdId(domainRgdId);
                results.add(mdDomain);
                mdDomain.setNotes(uniProtAccId+" part "+results.size());

                // is nucStopPos entirely within this CDS
                int nucStopPos = nucStartPos + lenToGo;
                if (nucStopPos < cdsLen) {
                    mdDomain.setStopPos(cf.getStartPos() + nucStopPos);
                    return results;
                } else {
                    // add this part of protein
                    mdDomain.setStopPos(cf.getStopPos());

                    int domainPartLen = mdDomain.getStopPos() - mdDomain.getStartPos() + 1;
                    nucStartPos = 0;
                    lenToGo -= domainPartLen;
                }
            }
        }
        logStrandProblem.info("plus strand problem "+uniProtAccId+" "+lenToGo);
        counters.increment("STRAND_PROBLEMS (details in strand_problem.log)");
        return results;
    }

    List<MapData> handleNegativeStrand(List<CodingFeature> cfs, int nucStartPos, int lenToGo, int domainRgdId, String uniProtAccId) {

        List<MapData> results = new ArrayList<>();

        // iterate over CDS features
        for( int i=cfs.size()-1; i>=0; i-- ) {
            CodingFeature cf = cfs.get(i);
            if (cf.getFeatureType() != TranscriptFeature.FeatureType.CDS ) {
                continue;
            }

            // we found a CDS
            int cdsLen = cf.getStopPos() - cf.getStartPos() + 1;
            if (nucStartPos >= cdsLen) {
                // nucStartPos is outside of this CDS
                nucStartPos -= cdsLen;
                continue;
            }

            MapData mdDomain = new MapData();
            mdDomain.setChromosome(cf.getChromosome());
            mdDomain.setMapKey(cf.getMapKey());
            mdDomain.setSrcPipeline(getSrcPipeline());
            mdDomain.setStrand(cf.getStrand());
            mdDomain.setRgdId(domainRgdId);
            results.add(mdDomain);
            mdDomain.setNotes(uniProtAccId+" part "+results.size());

            // nucStartPos is within this cds
            // is nucStopPos entirely within this CDS
            int stopPos = cf.getStopPos() - nucStartPos + 1;
            int startPos = stopPos - lenToGo + 1;
            if (startPos >= cf.getStartPos()) {
                // whole part is within this CDS (i.e. part of this CDS
                mdDomain.setStartPos(startPos);
                mdDomain.setStopPos(stopPos);
                return results;
            } else {
                // add this part of protein
                mdDomain.setStartPos(cf.getStartPos());
                mdDomain.setStopPos(stopPos);

                int domainPartLen = mdDomain.getStopPos() - mdDomain.getStartPos() + 1;
                nucStartPos = 0;
                lenToGo -= domainPartLen;
            }
        }
        logStrandProblem.info("minus strand problem "+uniProtAccId+" "+lenToGo);
        counters.increment("STRAND_PROBLEMS (details in strand_problem.log)");
        return results;
    }

    void load(List<ProteinDomain> domains) throws Exception {

        // merge protein domains
        for (ProteinDomain pd: domains) {

            if (pd.geInRgd == null) {
                GenomicElement ge = dao.insertDomainName(pd.getDomainName());
                // if ge has been inserted, its objectStatus property will be null
                if (ge.getObjectStatus() == null) {
                    counters.increment("DOMAINS INSERTED");
                } else {
                    counters.increment("DOMAINS UP-TO-DATE");
                }
                pd.geInRgd = ge;
            } else {
                counters.increment("DOMAINS UP-TO-DATE");
            }

            List<ProteinDomain> list = domainsMap.get(pd.geInRgd.getRgdId());
            if (list == null) {
                list = new ArrayList<>();
                domainsMap.put(pd.geInRgd.getRgdId(), list);
            }
            list.add(pd);
        }

        handleProteinToDomainAssociations(domains);
        handleDomainLoci();
    }

    void handleProteinToDomainAssociations(List<ProteinDomain> domains) throws Exception {

        log.debug("start associations");

        // create incoming associations as a set to filter out duplicates
        Collection<Association> assocIncoming = new HashSet<>();
        for (ProteinDomain domain : domains) {
            if (domain.geInRgd == null) {
                continue;
            }

            Protein protein = dao.getProteinByUniProtId(domain.uniprotAcc);

            Association assoc = new Association();
            assoc.setSrcPipeline("UNIPROTKB");
            assoc.setAssocType("protein_to_domain");
            assoc.setMasterRgdId(protein.getRgdId());
            assoc.setDetailRgdId(domain.geInRgd.getRgdId());
            assocIncoming.add(assoc);
        }

        log.debug("incoming assocs created");

        // get existing associations in RGD
        List<Association> assocInRgd = dao.getAssociationsByType("protein_to_domain", assembly.getSpeciesTypeKey());

        log.debug("in-rgd assocs loaded");

        Collection<Association> assocMatched = CollectionUtils.intersection(assocIncoming, assocInRgd);
        Collection<Association> assocToBeInserted = CollectionUtils.subtract(assocIncoming, assocInRgd);
        Collection<Association> assocToBeDeleted = CollectionUtils.subtract(assocInRgd, assocIncoming);

        log.debug("assocs qc-ed");

        dao.insertAssociations(assocToBeInserted);
        dao.deleteAssociations(assocToBeDeleted);

        if( !assocMatched.isEmpty() ) {
            counters.add("ASSOC protein_to_domain MATCHED ", assocMatched.size());
        }
        if( !assocToBeInserted.isEmpty() ) {
            counters.add("ASSOC protein_to_domain INSERTED", assocToBeInserted.size());
        }
        if( !assocToBeDeleted.isEmpty() ) {
            counters.add("ASSOC protein_to_domain DELETED ", assocToBeDeleted.size());
        }
    }


    void handleDomainLoci() throws Exception {

        log.debug(" loading protein domains ...");

        Set<MapDataEx> proteinDomainLociIncoming = new HashSet<>();

        List<Integer> proteinDomainRgdIds = new ArrayList<>(domainsMap.keySet());
        for( int proteinDomainRgdId: proteinDomainRgdIds ) {
            log.debug("processing PD " + proteinDomainRgdId);

            List<MapData> domainLoci = new ArrayList<>();
            for( ProteinDomain pd : domainsMap.get(proteinDomainRgdId) ) {
                if (pd.loci != null) {
                    for (MapData md : pd.loci) {
                        addDomainLoci(domainLoci, md);
                    }
                }
            }

            for( MapData md: domainLoci ) {
                proteinDomainLociIncoming.add(new MapDataEx(md));
            }
        }

        List<MapData> inRgd = dao.getProteinDomainLociForMapAndSource(assembly.getKey(), getSrcPipeline());
        List<MapData> duplicateLoci = new ArrayList<>();
        Set<MapDataEx> proteinDomainLociInRgd = new HashSet<>();
        for( MapData md: inRgd ) {
            if( !proteinDomainLociInRgd.add(new MapDataEx(md)) ) {
                duplicateLoci.add(md);
            }
        }
        log.debug("in-rgd loci loaded");
        counters.add("PROTEIN DOMAIN LOCI INCOMING ", proteinDomainLociIncoming.size());
        counters.add("PROTEIN DOMAIN LOCI IN RGD ", proteinDomainLociInRgd.size());

        if( !duplicateLoci.isEmpty() ) {
            dao.deleteMapData(duplicateLoci);

            counters.add("PROTEIN DOMAIN LOCI DELETED DUPLICATE", duplicateLoci.size());
        }


        Collection<MapDataEx> lociMatched = CollectionUtils.intersection(proteinDomainLociIncoming, proteinDomainLociInRgd);
        Collection<MapDataEx> lociToBeInserted = CollectionUtils.subtract(proteinDomainLociIncoming, proteinDomainLociInRgd);
        Collection<MapDataEx> lociToBeDeleted = CollectionUtils.subtract(proteinDomainLociInRgd, proteinDomainLociIncoming);

        log.debug("loci qc-ed");

        dao.insertMapData(new ArrayList<MapData>(lociToBeInserted));
        dao.deleteMapData(new ArrayList<MapData>(lociToBeDeleted));

        if( !lociMatched.isEmpty() ) {
            counters.add("PROTEIN DOMAIN LOCI MATCHED ", lociMatched.size());
        }
        if( !lociToBeInserted.isEmpty() ) {
            counters.add("PROTEIN DOMAIN LOCI INSERTED", lociToBeInserted.size());
        }
        if( !lociToBeDeleted.isEmpty() ) {
            counters.add("PROTEIN DOMAIN LOCI DELETED ", lociToBeDeleted.size());
        }
    }

    void addDomainLoci(List<MapData> loci, MapData md) {

        // see for duplicate loci, with same chr, strand, start and stop
        for( MapData mdLoci: loci ) {
            if( mdLoci.getChromosome().equals(md.getChromosome())
                    && mdLoci.getStrand().equals(md.getStrand())
                    && mdLoci.getStartPos().equals(md.getStartPos())
                    && mdLoci.getStopPos().equals(md.getStopPos()) ) {

                log.debug("-- protein-domain: merging duplicate loci for PROTEIN_DOMAIN_RGD:"+md.getRgdId());
                if( Utils.isStringEmpty(mdLoci.getNotes()) ) {
                    mdLoci.setNotes(md.getNotes());
                } else if( !Utils.isStringEmpty(md.getNotes()) ) {
                    if( !mdLoci.getNotes().contains(md.getNotes()) ) {
                        // merge notes alphabetically
                        String[] tokens = mdLoci.getNotes().split("; ");
                        Set<String> tokenSet = new TreeSet<>(Arrays.asList(tokens));
                        tokenSet.add(md.getNotes());
                        String newNotes = Utils.concatenate(tokenSet, "; ");
                        mdLoci.setNotes(newNotes);
                    }
                }
                return;
            }
        }

        loci.add(md);
    }

    public void setProcessedMapKeys(List<Integer> processedMapKeys) {
        this.processedMapKeys = processedMapKeys;
    }

    public List<Integer> getProcessedMapKeys() {
        return processedMapKeys;
    }

    public void setSrcPipeline(String srcPipeline) {
        this.srcPipeline = srcPipeline;
    }

    public String getSrcPipeline() {
        return srcPipeline;
    }

    /// this class is same as MapData, but it has refined equals() and hashCode() methods suitable for sorting of genomic positions on multiple assemblies
    class MapDataEx extends MapData {

        public MapDataEx(MapData md) {
            setKey(md.getKey());
            setMapKey(md.getMapKey());
            setRgdId(md.getRgdId());
            setChromosome(md.getChromosome());
            setStartPos(md.getStartPos());
            setStopPos(md.getStopPos());
            setStrand(md.getStrand());
            setSrcPipeline(md.getSrcPipeline());
            setNotes(md.getNotes());
        }

        public boolean equals(Object o) {
            MapData md = (MapData) o;
            return super.equals(o)
                && Utils.stringsAreEqual(getChromosome(), md.getChromosome())
                && Utils.intsAreEqual(getStartPos(), md.getStartPos())
                && Utils.intsAreEqual(getStopPos(), md.getStopPos())
                && Utils.stringsAreEqual(getNotes(), md.getNotes())
                && Utils.stringsAreEqual(getSrcPipeline(), md.getSrcPipeline());
        }

        public int hashCode() {
            return super.hashCode()
                ^ Utils.defaultString(getChromosome()).hashCode()
                ^ (getStartPos()==null ? 0 : getStartPos())
                ^ (getStopPos()==null ? 0 : getStopPos())
                ^ Utils.defaultString(getNotes()).hashCode()
                ^ Utils.defaultString(getSrcPipeline()).hashCode();
        }
    }
}
