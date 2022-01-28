#!/bin/bash

if [[ $# -gt 0 && "$1" == "--help" || "$1" == "-h" ]]; then
	/compi run -p /pipeline*.xml "$@"
	echo -e "\nAdditional CGA commands:"
	echo -e "\t- init_working_dir.sh /working_dir: initializes the working directory with a base parameters file."
	echo -e "\t- process_logs.sh /working_dir: processes the Compi logs directory to create a report of tasks with errors or warnings."
else
	if [ $# -gt 0 ] && [ "$1" == "resume" ]; then
		/compi resume -p /pipeline*.xml "${@:2}"
	else
		if [ $# -gt 0 ] && [ "$1" == "init_working_dir.sh" ]; then
			/scripts/init-working-dir.sh "${@:2}"
		else
			if [ $# -gt 0 ] && [ "$1" == "process_logs.sh" ]; then
				process_logs.sh "${@:2}"
			else
				/compi run -p /pipeline*.xml "$@"
			fi
		fi
	fi
fi
