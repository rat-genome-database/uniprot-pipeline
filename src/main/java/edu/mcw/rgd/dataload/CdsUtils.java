package edu.mcw.rgd.dataload;

/**
 * Created by mtutaj on 3/22/2018.
 */

import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.TranscriptFeature;
import edu.mcw.rgd.process.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by mtutaj on 9/20/2017.
 * utility class to generate CDS features based on EXONS and UTR3,UTR5 features
 */
public class CdsUtils {

    private TranscriptDAO dao;
    private Integer mapKey;

    public CdsUtils(TranscriptDAO dao, Integer mapKey) {
        this.dao = dao;
        this.mapKey = mapKey;
    }

    public List<CodingFeature> buildCfList(MapData trMd) throws Exception {

        String strand = trMd.getStrand();

        List<TranscriptFeature> transcriptFeats = getTranscriptFeatures(trMd);
        List<TranscriptFeature> exonList = new ArrayList<>(transcriptFeats.size());
        TranscriptFeature utr5 = null;
        TranscriptFeature utr3 = null;

        for(TranscriptFeature transf : transcriptFeats){
            if(transf.getFeatureType()== TranscriptFeature.FeatureType.UTR5){
                utr5 = transf;
            }
            if(transf.getFeatureType()== TranscriptFeature.FeatureType.EXON){
                exonList.add(transf);
            }
            if(transf.getFeatureType()== TranscriptFeature.FeatureType.UTR3){
                utr3 = transf;
            }
        }
        if( exonList.isEmpty() )
            return Collections.emptyList();


        List<CodingFeature> cfList = new ArrayList<>();

        boolean wasUtr = false;
        boolean wasCds = false;

        for( TranscriptFeature exon: exonList ){
            // copy exon
            CodingFeature cf = new CodingFeature(exon);
            cf.setFeatureType(TranscriptFeature.FeatureType.EXON);
            cfList.add(cf);

            if(strand.equals("+")){

                // handle UTR5
                if( !wasUtr ) {
                    // utr5 region was not processed yet entirely
                    if( utr5==null ) {
                        // there is no utr5 region, proceed to cds region
                        wasUtr = true;
                    } else {
                        // there is utr5 region -- see if utr5 covers entire exon
                        if( utr5.getStartPos()<=exon.getStartPos() && utr5.getStopPos()>=exon.getStopPos() ) {
                            // add utr5 object
                            cf = new CodingFeature(exon);
                            cf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                            cfList.add(cf);
                            continue;
                        } else {
                            // partial utr5, partial cds
                            if( utr5.getStopPos()>=exon.getStartPos() ) {
                                cf = new CodingFeature(exon);
                                cf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                                cf.setStopPos(utr5.getStopPos());
                                cfList.add(cf);
                            }
                            wasUtr = true;
                        }
                    }
                }

                // handle CDS
                if( !wasCds ) {

                    cf = new CodingFeature(exon);
                    cf.setFeatureType(TranscriptFeature.FeatureType.CDS);
                    cfList.add(cf);

                    // remove utr5 block
                    if( utr5!=null && utr5.getStopPos()>=exon.getStartPos() ) {
                        cf.setStartPos(utr5.getStopPos()+1);
                    }
                    // remove utr3 block
                    if( utr3!=null && utr3.getStartPos()<=exon.getStopPos() ) {
                        cf.setStopPos(utr3.getStartPos()-1);
                        wasCds = true;
                    }

                    if( cf.getStartPos() > cf.getStopPos() ) {
                        // empty CDS
                        cfList.remove(cfList.size()-1);
                    }
                }

                // handle UTR3
                if( utr3!=null ) {
                    // full utr3
                    if( utr3.getStartPos()<=exon.getStartPos() && utr3.getStopPos()>=exon.getStopPos() ) {
                        cf = new CodingFeature(exon);
                        cf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                        cfList.add(cf);
                    }
                    else // partial utr3
                        if( utr3.getStartPos()<=exon.getStopPos() && utr3.getStopPos()>=exon.getStopPos() ) {
                            cf = new CodingFeature(exon);
                            cf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                            cf.setStartPos(utr3.getStartPos());
                            cfList.add(cf);
                        }
                }
            }else if(strand.equals("-")){

                // handle UTR3
                if( !wasUtr ) {
                    // utr3 region was not processed yet entirely
                    if( utr3==null ) {
                        // there is no utr3 region, proceed to cds region
                        wasUtr = true;
                    } else {
                        // there is utr3 region -- see if utr5 covers entire exon
                        if( utr3.getStartPos()<=exon.getStartPos() && utr3.getStopPos()>=exon.getStopPos() ) {
                            // add utr5 object
                            cf = new CodingFeature(exon);
                            cf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                            cfList.add(cf);
                            continue;
                        } else {
                            // partial utr3, partial cds
                            if( utr3.getStopPos()>=exon.getStartPos() ) {
                                cf = new CodingFeature(exon);
                                cf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                                cf.setStopPos(utr3.getStopPos());
                                cfList.add(cf);
                            }
                            wasUtr = true;
                        }
                    }
                }

                // handle CDS
                if( !wasCds ) {

                    cf = new CodingFeature(exon);
                    cf.setFeatureType(TranscriptFeature.FeatureType.CDS);
                    cfList.add(cf);

                    // remove utr3 block
                    if( utr3!=null && utr3.getStopPos()>=exon.getStartPos() ) {
                        cf.setStartPos(utr3.getStopPos()+1);
                    }
                    // remove utr5 block
                    if( utr5!=null && utr5.getStartPos()<=exon.getStopPos() ) {
                        cf.setStopPos(utr5.getStartPos()-1);
                        if( cf.getStartPos()>cf.getStopPos() ) {
                            // empty CDS
                            cfList.remove(cfList.size()-1);
                        }
                        // empty CDS or not, we are again in utr area
                        wasCds = true;
                    }
                }

                // handle UTR5
                if( utr5!=null ) {
                    // full utr5
                    if( utr5.getStartPos()<=exon.getStartPos() && utr5.getStopPos()>=exon.getStopPos() ) {
                        cf = new CodingFeature(exon);
                        cf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                        cfList.add(cf);
                    }
                    else // partial utr5
                        if( utr5.getStartPos()<=exon.getStopPos() && utr5.getStopPos()>=exon.getStopPos() ) {
                            cf = new CodingFeature(exon);
                            cf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                            cf.setStartPos(utr5.getStartPos());
                            cfList.add(cf);
                        }
                }
            }
        }

        validateCfList(cfList);
        computeCodingPhase(cfList);
        return cfList;
    }

    void validateCfList(List<CodingFeature> cfList) {

        // QC: report any features that start pos <= stop pos
        for( CodingFeature cf: cfList ) {
            if( cf.getStartPos() > cf.getStopPos() ) {
                System.out.println("problem with coding feature trId="+cf.getTranscriptRgdId());
                break;
            }
        }
    }

    List<TranscriptFeature> getTranscriptFeatures(MapData trMd) throws Exception {

        String key = trMd.getRgdId()+"|"+mapKey;
        List<TranscriptFeature> tfs = _featureMap.get(key);
        if( tfs==null ) {
            tfs = dao.getFeatures(trMd.getRgdId(), mapKey);
            _featureMap.put(key, tfs);
        }

        // note: multiple loci could be returned
        // filter out those outside range
        Iterator<TranscriptFeature> it = tfs.iterator();
        while( it.hasNext() ) {
            TranscriptFeature tf = it.next();
            if( !transcriptPositionOverlapsGenePosition(trMd, tf) ) {
                it.remove();
            }
        }
        return tfs;
    }

    static ConcurrentHashMap<String, List<TranscriptFeature>> _featureMap = new ConcurrentHashMap<>();

    public boolean transcriptPositionOverlapsGenePosition(MapData md1, MapData md2) {
        // map keys must match
        if( !md1.getMapKey().equals(md2.getMapKey()) )
            return false;
        // chromosomes must match
        if( !Utils.stringsAreEqualIgnoreCase(md1.getChromosome(), md2.getChromosome()) )
            return false;
        // positions must overlap
        if( md1.getStopPos() < md2.getStartPos() )
            return false;
        if( md2.getStopPos() < md1.getStartPos() )
            return false;
        return true;
    }

    public void computeCodingPhase(List<CodingFeature> cfList) {
        List<CodingFeature> cdsList = new ArrayList<>(cfList.size());
        for( CodingFeature cf: cfList ) {
            if( cf.getFeatureType()==TranscriptFeature.FeatureType.CDS )
                cdsList.add(cf);
        }
        if( cdsList.isEmpty() )
            return;

        int phase = 0;

        if(cdsList.get(0).getStrand().equals("+")){
            cdsList.get(0).setCodingPhase(phase);

            for(int p=0; p<cdsList.size()-1; p++){
                //-----------adding phase information------
                phase = (3 - (cdsList.get(p).getStopPos()-cdsList.get(p).getStartPos()+1-phase)%3)%3;
                cdsList.get(p+1).setCodingPhase(phase);
            }
        }

        if(cdsList.get(0).getStrand().equals("-")){
            cdsList.get(cdsList.size()-1).setCodingPhase(phase);

            for(int j=(cdsList.size()-1); j>0; j--) {
                //-----------adding phase information------
                phase = (3 - (cdsList.get(j).getStopPos()-cdsList.get(j).getStartPos()+1-phase)%3)%3;
                cdsList.get(j-1).setCodingPhase(phase);
            }
        }
    }
}

class CodingFeature extends TranscriptFeature{
    int codingPhase=0;

    public int getCodingPhase() {
        return codingPhase;
    }

    public void setCodingPhase(int codingPhase) {
        this.codingPhase = codingPhase;
    }

    public String getCodingPhaseStr() {
        return codingPhase<0 ? "." : Integer.toString(codingPhase);
    }

    public CodingFeature(TranscriptFeature tf) {
        setTranscriptRgdId(tf.getTranscriptRgdId());
        setChromosome(tf.getChromosome());
        setStrand(tf.getStrand());
        setRgdId(tf.getRgdId());
        setMapKey(tf.getMapKey());
        setCodingPhase(-1);
        setStartPos(tf.getStartPos());
        setStopPos(tf.getStopPos());
    }
}
