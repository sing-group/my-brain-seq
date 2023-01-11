
# test if file is locked, then cp the script to working-dir
function cp_and_lock {
# $1 : script to copy  # $2 : task name
	(
	flock 200 || echo "[PIPELINE -- "${2}"]: ${1} is locked, cp omitted."
	#Makes a copy of the scripts used in the analysis to working-dir
	cp -n ${scriptsDir}/${1} ${workingDir}/compi_scripts/${1}
	) 200>/var/lock/${1}.lock
}