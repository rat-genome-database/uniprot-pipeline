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
    public String notes;
    public String domainSeq;

    public GenomicElement geInRgd;
    public List<MapData> loci;

    public void cleanupDomainName() {

        int pos;
        if( domainName.lastIndexOf('}')>0 ) {
            // is there a matching '{'
            pos = domainName.lastIndexOf('{');
            if( pos > 0 ) {
                domainName = domainName.substring(0, pos).trim();
            } else {
                pos = domainName.lastIndexOf(' ');
                domainName = domainName.substring(0, pos).trim();
            }
        }

        pos = domainName.lastIndexOf("ProRule:");
        if( pos>0 ) {
            domainName = domainName.substring(0, pos).trim();
        }
    }

    static Pattern domainNameCounterPattern = Pattern.compile(" \\d+$");

    public void qcDomainName() {
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
