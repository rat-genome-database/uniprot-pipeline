# script to run UniProtKB pipeline
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu,jthota@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=rgd.developers@mcw.edu,rgd.pipelines@mcw.edu
fi

cd $APP_HOME

#initialize summary.log
echo "===" > $APP_HOME/logs/summary.log

speciesList=( "rat" "mouse" "human" "dog" "bonobo" "squirrel" "chinchilla" "pig" )

for species in "${speciesList[@]}"; do
    $APP_HOME/_run.sh -species "$species" $@  2>&1
    cat $APP_HOME/logs/main_summary.log >> $APP_HOME/logs/summary.log
done

mailx -s "[$SERVER] UniProtKB data loading is done" $EMAIL_LIST < $APP_HOME/logs/summary.log
