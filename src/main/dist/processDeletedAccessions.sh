# download files with deleted accessions from UniProtKB
# validate them against protein objects in RGD
# and withdraw the protein objects in RGD that became deleted at UniProtKB
#
APP_HOME=/home/rgddata/pipelines/UniProtPipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,jthota@mcw.edu,jrsmith@mcw.edu
fi

$APP_HOME/_run.sh --deletedAccessions $@  2>&1  > cron3.log

mailx -s "[$SERVER] UniProtKB: deleted accessions" $EMAIL_LIST< $APP_HOME/logs/deleted_acc_summary.log
