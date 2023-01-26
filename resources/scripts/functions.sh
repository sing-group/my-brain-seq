
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

# get output path for the files saved in 6_deseq2, and 6_edger directories
# to use this function assign it to a variable like this: result=$(get_output_dir X Y)
function get_output_dir {
	# $1 : software used (deseq/edger)  # $2 : ${comparison}
	local comparison_label=$(echo "$2" | cut -d'=' -f1 | tr -d \" | sed -e 's/ /\\ /g' | xargs)
	if [[ $1 == 'deseq' ]];     then local soft_dir="${dsqOut}";
	elif [[ $1 == 'edger' ]];   then local soft_dir="${edgOut}"; fi
	local path_output="${workingDir}/${outDir}/${soft_dir}/${comparison_label}/"
	echo "$path_output"
}
