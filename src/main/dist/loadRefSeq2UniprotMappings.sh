# for PROTEIN objects, load corresponding RefSeq accession ids
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST="mtutaj@mcw.edu jthota@mcw.edu jrsmith@mcw.edu"
fi

$APP_HOME/_run.sh --loadRefSeq2UniProt $@  2>&1  > cron2.log

mailx -s "[$SERVER] UniProtKB: RefSeq2UniProt mappings loaded" $EMAIL_LIST< cron2.log
