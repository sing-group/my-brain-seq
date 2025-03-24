#!/bin/bash
set -o nounset
set -o errexit

if [ $# -ne 2 ]; then
	show_error "[ERROR]: This script requires two arguments (the path to the project dir and the new version)"
	exit 1
fi

PROJECT_DIR=$1
NEW_VERSION=$2
OLD_VERSION=$(cat ${PROJECT_DIR}/pipeline.xml | grep '<version>' | sed 's/<version>//; s/<\/version>//' | tr -d $'\t')

echo "Updating version from ${OLD_VERSION} to ${NEW_VERSION}"

sed -i "s#<version>.*</version>#<version>${NEW_VERSION}</version>#" ${PROJECT_DIR}/pipeline.xml
sed -i "s#^MYBRAIN_SEQ_VERSION=.*#MYBRAIN_SEQ_VERSION=\${MYBRAIN_SEQ_VERSION-${NEW_VERSION}}#" resources/init-working-dir/run.sh

year=$(date '+%Y')
echo "Updating the footer of the visual console with the new version and the current year"

mv ${PROJECT_DIR}/resources/visual_console/foot.txt ${PROJECT_DIR}/resources/visual_console/foot.txt.old
head -n -1 ${PROJECT_DIR}/resources/visual_console/foot.txt.old > ${PROJECT_DIR}/resources/visual_console/foot.txt
rm ${PROJECT_DIR}/resources/visual_console/foot.txt.old
echo "..............................................................................." > ${PROJECT_DIR}/resources/visual_console/foot.txt
echo "                                              © MyBrain-Seq v${NEW_VERSION}, 2020-${year}" >> ${PROJECT_DIR}/resources/visual_console/foot.txt

echo "You should also have a look at the following files as they contain the old version (${OLD_VERSION})":
rgrep -l --exclude-dir './image-files' --exclude '*.R' ${OLD_VERSION} ./*
