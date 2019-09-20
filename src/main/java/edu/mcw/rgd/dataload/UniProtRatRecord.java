package edu.mcw.rgd.dataload;

import java.util.*;

/**
 * @author mtutaj
 * @since Apr 6, 2010
 * represents a single record as read from from source file
 */
public class UniProtRatRecord {

    String srcPipeline;

    /// list of all external db links read from source file
    Map<String, List<String>> xdbMap = new HashMap<>();

    /// record nr, as read from the source file; first record has recno 1
    int recNo;

    /// flag to mark the record data has been dumped to database pipeline logging system
    boolean dumpedToDatabase;

    /// primary accession id
    String uniProtAccId;

    // /matching result: map of rgd ids to list of xdb ids
    List<Integer> matchingRgdIds = new ArrayList<>();

    // recommended protein name, or if not available, submitted protein name
    String proteinName;

    String geneName; // as found in 'GN' field

    String proteinSequence;

    public UniProtRatRecord(String srcPipeline) {
        setSrcPipeline(srcPipeline);
    }

    public void addEntry(String xdbName, String accId, String link1, String link2) {

        // lazy-create of the list for given xdbName
        List<String> xdbList = xdbMap.get(xdbName);
        if( xdbList==null ) {
            xdbList = new ArrayList<>();
            xdbMap.put(xdbName, xdbList);
        }

        // always add ACC_ID
        xdbList.add(accId);

        // combine link1 and link2 into one valid link text
        String vlink1 = (link1!=null && link1.length()>=3) ? link1 : null;
        String vlink2 = (link2!=null && link2.length()>=3) ? link2 : null;
        if( vlink1!=null && vlink2!=null ) {
            vlink1 += "; "+vlink2;
        }
        // if LINK is the same as ACC_ID, no need to create LINK field
        if( vlink1!=null && vlink1.equals(accId) ) {
            vlink1 = null;
        }
        xdbList.add(vlink1);
    }

    public void addEnsemblEntry(String ensemblID, boolean isHuman) {
        // it could be one of this: ENSRNOT00000016981; ENSRNOP00000016981; ENSRNOG00000010945 for rat
        //   or  ENSSSCT00000017610; ENSSSCP00000017132; ENSSSCG00000016174. for pig
        // so for non-human genes the mask is  len=18
        //   genes:       ENS???G???????????
        //   proteins:    ENS???P???????????
        //   transcripts: ENS???T???????????
        //and for human:
        //   genes:       ENSG??????????? len=15
        //   proteins:    ENSP???????????
        //   transcripts: ENST???????????
        if( ensemblID!=null && ensemblID.startsWith("ENS") ) {

            if( ensemblID.length()<15 ) {
                System.out.println("WARN: unexpected Ensembl Acc: "+ensemblID);
                return;
            }
            char idType = isHuman ? ensemblID.charAt(3) : ensemblID.charAt(6);
            if( idType=='T' ) {
                addEntry("EnsemblT", ensemblID, ensemblID, null);
            } else if( idType=='P' ) {
                addEntry("EnsemblP", ensemblID, ensemblID, null);
            } else if( idType=='G' ) {
                addEntry("EnsemblG", ensemblID, ensemblID, null);
            }
        }
    }

    // dump all data as one string
    public String makeRecAsString() {

        // to prevent outputting multiple instance of data dump to database,
        // we are checking UniProtRatRecord flag to see if we already put the data into db
        if( isDumpedToDatabase() )
            return null;

        // create new record dump
        StringBuilder buf = new StringBuilder(8000);
        Map<String, List<String>> xdbMap = getXdbMap();
        for( Map.Entry<String, List<String>> entry: xdbMap.entrySet() ) {
            String xdbName = entry.getKey();
            buf.append("<b>").append(xdbName).append(":</b> ");  // xdb name made bold
            List<String> xdbList = entry.getValue();
            // on list, process by pairs, first goes the xid, then link-info
            for( int i=0; i<xdbList.size(); i+=2 ) {
                String accId = xdbList.get(i);
                if( accId==null )
                    continue;
                buf.append(accId);
                String linkText = xdbList.get(i+1);
                if( linkText!=null && linkText.length()>3 )
                    buf.append('#').append(linkText);
                buf.append(", ");
            }
            buf.append("<br/>\n");
        }

        // mark as dumped
        setDumpedToDatabase(true);

        // return the dump as string
        return buf.toString();
    }

    public void parseSecondaryAccessionIds(String line) {
        for( String id: line.split("[; ]+") ) {
            if( !id.isEmpty() ) {
                this.addEntry("UniProtSecondary", id, id, null);
            }
        }
    }

    static List<String> _emptyList = new ArrayList<>();
    public List<String> getXdbInfo(String xdbName) {
        List<String> list = xdbMap.get(xdbName);
        return list==null ? _emptyList : list;
    }

    public String getSrcPipeline() {
        return srcPipeline;
    }

    public void setSrcPipeline(String srcPipeline) {
        this.srcPipeline = srcPipeline;
    }

    public Map<String, List<String>> getXdbMap() {
        return xdbMap;
    }

    public int getRecNo() {
        return recNo;
    }

    public void setRecNo(int recNo) {
        this.recNo = recNo;
    }

    public boolean isDumpedToDatabase() {
        return dumpedToDatabase;
    }

    public void setDumpedToDatabase(boolean dumpedToDatabase) {
        this.dumpedToDatabase = dumpedToDatabase;
    }

    public String getUniProtAccId() {
        return uniProtAccId;
    }

    public void setUniProtAccId(String uniProtAccId) {
        this.uniProtAccId = uniProtAccId;
    }

    public String getProteinName() {
        return proteinName;
    }

    public void setProteinName(String proteinName) {
        this.proteinName = proteinName;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }
}
