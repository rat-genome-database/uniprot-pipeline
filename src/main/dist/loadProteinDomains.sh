# script to run UniProtKB pipeline
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu
fi

cd $APP_HOME

speciesList=( "rat" "mouse" "human" "dog" "bonobo" "squirrel" "chinchilla" "pig" )

for species in "${speciesList[@]}"; do
    $APP_HOME/_run.sh -species "$species" --loadProteinDomains $@  2>&1
done

mailx -s "[$SERVER] UniProtKB protein domain done" $EMAIL_LIST < $APP_HOME/logs/domains.log
