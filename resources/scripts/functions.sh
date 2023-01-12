
# test if file is locked, then cp the script to working-dir
function cp_and_lock {
# $1 : script to copy  # $2 : task name  # $3 : cp source directory
	(
	flock 200 || echo "[PIPELINE -- "${2}"]: ${1} is locked, cp omitted."
	# selects the target directory using the source directory as reference
	if [[ $3 == ${scriptsDir} ]];     then target_dir='/compi_scripts/';
	elif [[ $3 == ${databasesDir} ]]; then target_dir='/compi_databases/'; fi
	#Makes a copy of the scripts used in the analysis to working-dir
	cp -n ${3}/${1} ${workingDir}${target_dir}${1}
	) 200>/var/lock/${1}.lock
}
