package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.GenomicElement;
import edu.mcw.rgd.datamodel.MapData;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by mtutaj on 1/31/2018.
 */
public class ProteinDomain {

    public String uniprotAcc;
    public int aaStartPos;
    public int aaStopPos;
    private String domainName;

    public GenomicElement geInRgd;
    public List<MapData> loci;


    void cleanupDomainName() {

        // remove comments from the domain name
        //  EGF-like; calcium-binding       --> EGF-like
        //  Sushi 1; atypical; lacks a Cys  --> Sushi 1
        // TODO: consider incorpating comments into notes of the locus in MAPS_DATA table
        int commentPos = domainName.indexOf("; ");
        if( commentPos>0 ) {
            domainName = domainName.substring(0, commentPos);
        }
    }

    static Pattern domainNameCounterPattern = Pattern.compile(" \\d+$");

    public void qcDomainName() {

        cleanupDomainName();

        // many proteins have the same domain appearing multiple times
        // f.e. P58365 protein has Cadherin domain appearing 27 times!
        //   so the domain names are like this: 'Cadherin 1', 'Cadherin 2', ... 'Cadherin 27'
        // therefore we have to strip the last part and keep only the domain name f.e. 'Cadherin'
        Matcher m = domainNameCounterPattern.matcher(getDomainName());
        int counterPos = -1;
        while (m.find()) {
            counterPos = m.start();
        }
        if (counterPos > 0) {
            setDomainName(getDomainName().substring(0, counterPos));
        }
    }

    public String getDomainName() {
        return domainName;
    }

    public void setDomainName(String domainName) {
        this.domainName = domainName;
    }
}
