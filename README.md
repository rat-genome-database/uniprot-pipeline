# uniprot-pipeline
Imports external db ids, protein objects and sequences from UniProtKB.

PROTEIN LOADER

- load protein_to_gene associations (into RGD_ASSOCIATIONS table)
- load secondary uniprot xdb ids (into RGD_ACC_XDB table);
  secondary accession ids are imported as xdb ids with XDB_KEY=60
- load protein sequences (into RGD_SEQUENCES table);
  * the canonical protein sequences (available in the incoming data)
  are loaded with seq_type 'uniprot_seq'
  * the old protein sequences are created with seq_type 'old_uniprot_seq';
  they are created when the incoming sequence differs from the sequence in the database;
  sequences of type 'old_uniprot_seq' combined are considered the 'consensus sequence history'


FILE PARSER

- processing of AC lines:
  * primary uniprot acc id is the first acc id on the first AC line
  * when a record had multiple AC lines, the other AC lines must be ignored;
- processing os OS lines: (species)
  * logic must be aware that OS value can span multiple lines
- processing of Gene3D accession ids:
  * they should be loaded just as they are; however in the past, Gene3D acc ids were in the format 'G3DSA:XX.XX.XX.XX'

added processing of chinchilla, bonobo, dog and squirrel


fixed matching logic to exclude self-matching (matching by UniProtKB ids brought in
  by UniProtKB pipeline in the past)

deletion of UniProt id results in insertion of an 'old_protein_id' alias


MATCHING
1) matching algorithm change:
  - incoming data is matched to RGD data in the following order (as below);
  matching is terminated as soon there is some matching data
  - OLD logic:
    RGD, GeneId, HGNC, MGI, UniProt, RefSeq Nucl, Ensembl
  - NEW logic:
    GeneId, HGNC, MGI, RefSeq Nucl, UniProt, Ensembl
2) improved handling of multis:
 - where one UniProtKB id matches multiple genes in RGD
   previous logic was not properly handling of up to 0.1% of uniprot ids
3) discontinued ids_merged.log
  - ids_inserted and ids_deleted logs until June 22, 2015 contain a lot of duplication:
  due to code flaw some entries are unnecessarily inserted, and during next pipeline run they are deleted and so on;
  therefore as of June 22, 2015 these logs are discontinued;
  inserted.log and deleted.log are used instead