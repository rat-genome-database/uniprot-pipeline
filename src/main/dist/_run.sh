#!/usr/bin/env bash
#
# wrapper script
#
. /etc/profile
APPNAME=UniProtPipeline

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

bin/$APPNAME -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar $APPNAME.jar "$@" 2>&1