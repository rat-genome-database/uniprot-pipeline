package edu.mcw.rgd.dataload;

import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.XdbId;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/*** some SwissProt acc ids are being assigned incorrectly by the current algorithm;
 * this code patches
 */
public class SwissProtFixer {

    private Map<String,Integer> ids; // map of Swiss-Prot ids to gene rgd ids

    public void run() throws Exception {

        // currently disable this module
        if( true ) {
            return;
        }

        UniProtDAO dao = new UniProtDAO();
        int swissProtIdsNotInRgd = 0;

        System.out.println("mappings from config file processed: "+ids.size());

        List<XdbId> idsForDelete = new ArrayList<>();

        for( Map.Entry<String,Integer> entry: ids.entrySet() ) {
            String swissProtId = entry.getKey();
            int geneRgdId = entry.getValue();

            List<XdbId> xdbIds = dao.getXdbIds(XdbId.XDB_KEY_UNIPROT, swissProtId, "UniProtKB/Swiss-Prot");
            if( xdbIds.isEmpty() ) {
                swissProtIdsNotInRgd++;
            } else {
                for( XdbId xdbId: xdbIds ) {
                    if( xdbId.getRgdId()!=geneRgdId ) {
                        idsForDelete.add(xdbId);
                    }
                }
            }
        }
        dao.deleteXdbIds(idsForDelete, RgdId.OBJECT_KEY_GENES);

        System.out.println("xdb ids deleted: "+idsForDelete.size());
        System.out.println("swiss prot ids not in RGD: "+swissProtIdsNotInRgd);
    }

    public Map<String, Integer> getIds() {
        return ids;
    }

    public void setIds(Map<String, Integer> ids) {
        this.ids = ids;
    }
}
