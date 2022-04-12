# load canonical proteins for species defined in properties file
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu
fi

cd $APP_HOME

$APP_HOME/_run.sh --loadCanonicalProteins $@  > $APP_HOME/canonical_proteins.log

mailx -s "[$SERVER] UniProtKB canonical proteins done" $EMAIL_LIST < $APP_HOME/canonical_proteins.log