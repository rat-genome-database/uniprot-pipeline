#!/usr/bin/env bash
#
# wrapper script to call gradle-generated script
#
. /etc/profile
APPNAME=UniProtPipeline

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR
DB_OPTS="-Dspring.config=$APPDIR/../properties/default_db.xml"
LOG4J_OPTS="-Dlog4j.configuration=file://$APPDIR/properties/log4j.properties"
declare -x "UNI_PROT_PIPELINE_OPTS=$DB_OPTS $LOG4J_OPTS"
bin/$APPNAME "$@" 2>&1