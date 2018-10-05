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
    they are created when the incoming sequence differs from the sequence in the database
