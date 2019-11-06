# script to run UniProtKB pipeline
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu
fi

cd $APP_HOME

# processed species and assemblies are read from properties/AppConfigure.xml
$APP_HOME/_run.sh -species "$species" --loadProteinDomains $@  2>&1

mailx -s "[$SERVER] UniProtKB protein domain done" $EMAIL_LIST < $APP_HOME/logs/domains.log
