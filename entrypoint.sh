#!/bin/bash

if [[ $# -gt 0 && "$1" == "--help" || "$1" == "-h" ]]; then
	/compi run -p /pipeline*.xml "$@"
	echo -e "\t- init_working_dir.sh /working_dir: initializes the working directory with a base parameters file."
	echo -e "\t- make_run-sh.sh /path/to/compi.parameters: build the runner for myBrain-Seq analysis."
	echo -e "\t- visual_console.sh: terminal user interface that allows to initialize the working directory, 
	build the runner and run the analysis interactively."

elif [ $# -gt 0 ] && [ "$1" == "resume" ]; then
	/compi resume -p /pipeline*.xml "${@:2}"

elif [ $# -gt 0 ] && [ "$1" == "init_working_dir.sh" ]; then
	/scripts/init_working_dir.sh "${@:2}"

elif [ $# -gt 0 ] && [ "$1" == "run" ]; then
	/init-working-dir/run.sh "${@:2}"

elif [ $# -gt 0 ] && [ "$1" == "visual_console.sh" ]; then
	/visual_console/visual_console.sh "${@:2}"

else
	/compi run -p /pipeline*.xml "$@"

fi
