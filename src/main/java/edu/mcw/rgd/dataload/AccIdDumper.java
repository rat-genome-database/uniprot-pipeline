package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.Utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 5/22/12
 * Time: 11:23 AM
 * <p>
 * dumps a list of primary to secondary accession ids;
 * tab separated format;
 * 1st column contains primary accession id
 * 2nd column contains a list of secondary acc ids separated by '|'
 */
public class AccIdDumper {

    protected BufferedWriter writer;

    /**
     * open the file for writing and write a number of comments in the header
     * @param fileName file name to be written (previous contents is overwritten)
     * @throws IOException
     */
    public void writeFileHeader(String fileName) throws IOException {

        writer = new BufferedWriter(new FileWriter(fileName.toLowerCase()));

        writer.write("#list of primary to secondary accession ids - UniProtKB, SPROT");
        writer.newLine();
        writer.write("#tab separated format");
        writer.newLine();
        writer.write("#1st column contains primary accession id");
        writer.newLine();
        writer.write("#2nd column contains a list of secondary acc ids separated by '|'");
        writer.newLine();
    }

    /**
     * write a new accession id  and secondary ids into file
     * @param primaryAccId primary accession id
     * @param secondaryAccIds list of secondary accession ids, separated
     * @throws IOException
     */
    public void writeData(String primaryAccId, String secondaryAccIds) throws IOException {

        String[] ids = secondaryAccIds.split("[; ]+");

        writer.write(primaryAccId);
        writer.write('\t');
        writer.write(Utils.concatenate(Arrays.asList(ids), "|"));
        writer.newLine();
    }

    /**
     * close the file
     * @throws IOException
     */
    public void close() throws IOException {

        writer.close();
    }
}
