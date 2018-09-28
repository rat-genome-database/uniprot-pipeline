package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.GenomicElement;
import edu.mcw.rgd.datamodel.MapData;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by mtutaj on 1/31/2018.
 */
public class ProteinDomain {

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

    public String getDomainName() {
        return domainName;
    }

    public void setDomainName(String domainName) {
        this.domainName = domainName;
    }
}
