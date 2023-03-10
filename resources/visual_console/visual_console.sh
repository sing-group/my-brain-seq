#!/bin/bash
################################################################################
#                    ____               _                _____                 #
#                   |  _ \             (_)              / ____|                #
#  _ __ ___   _   _ | |_) | _ __  __ _  _  _ __  ______| (___    ___   __ _    #
# | '_ ` _ \ | | | ||  _ < | '__|/ _` || || '_ \|______|\___ \  / _ \ / _` |   #
# | | | | | || |_| || |_) || |  | (_| || || | | |       ____) ||  __/| (_| |   #
# |_| |_| |_| \__, ||____/ |_|   \__,_||_||_| |_|      |_____/  \___| \__, |   #
#              __/ |                                           __        | |   #
#             |___/     _  __      ____ ___   ___   ___ ___   / /___     |_|©  #
#                      | |/ /  _  / __// _ \ / _ \ (_-</ _ \ / // -_)          #
#                      |___/  (_) \__/ \___//_//_//___/\___//_/ \__/           #
#                                                                              #
################################################################################
set -o nounset
set -o errexit

header='/visual_console/header.txt'
foot='/visual_console/foot.txt'
curr_wd=' Working directory: [ NOT SELECTED ]'

temp_file='/tmp/mbs_visual-console_temp_file.txt'
temp_file_exists='FALSE'

# If the temporary file don't exist, creates it
if [[ -f "${temp_file}" ]]; then
    working_dir=$(cut -d '=' -f 2 "${temp_file}")
    temp_file_exists='TRUE'
fi

# Colors
ESC=$(printf '\033') RESET="${ESC}[0m" BLACK="${ESC}[30m" RED="${ESC}[31m"
GREEN="${ESC}[32m" YELLOW="${ESC}[33m" BLUE="${ESC}[34m" MAGENTA="${ESC}[35m"
GREY="${ESC}[1;30m" WHITE="${ESC}[37m" DEFAULT="${ESC}[39m"

# Color Functions
function greenprint { printf "${GREEN}%s${RESET}\n" "$1"; }
function blueprint { printf "${BLUE}%s${RESET}\n" "$1"; }
function redprint { printf "${RED}%s${RESET}\n" "$1"; }
function yellowprint { printf "${YELLOW}%s${RESET}\n" "$1"; }
function magentaprint { printf "${MAGENTA}%s${RESET}\n" "$1"; }
function greyprint { printf "${GREY}%s${RESET}\n" "$1"; }
function fn_bye { clear; exit 0; }
function fn_fail { echo "Wrong option." exit 1; }
function enter_to_continue { echo; echo 'Press [ENTER] to continue'; read; }

function clear_view {
    clear
    cat "${header}"
    greyprint "${curr_wd}"
    #greenprint $(cat "${header}")
}

function update_wd {
    if [[ ! -z "${working_dir+x}" ]]; then
        curr_wd=" Working directory: ${working_dir}"
    fi
}

function print_foot { cat "${foot}"; }

function exit_options {
        echo "Proceed?"
        case `select_opt "Yes" "No" "Exit"` in
        0)
            exit_loop='y'
            ;;
        1)
            exit_loop='n'
            ;;
        2)
            exit_loop='exit'
            ;;
        esac
}

function wd_initialization_menu {
    clear_view
    echo
    echo -ne "$(greenprint 'WORKING DIRECTORY INITIALIZATION:')"
    echo
    echo "Please, paste the path where the working directory should be placed."
    echo
    echo "Working directory:"
    echo -ne "    > "
    read working_dir
    echo

    # display options
    exit_options
}

# TO MAKE THE SELECTIONS USING THE ARROW KEYS
function select_option {
    # little helpers for terminal print control and key input
    ESC=$( printf "\033")
    cursor_blink_on()  { printf "$ESC[?25h"; }
    cursor_blink_off() { printf "$ESC[?25l"; }
    cursor_to()        { printf "$ESC[$1;${2:-1}H"; }
    print_option()     { printf "   $1 "; }
    print_selected()   { printf "  $ESC[7m $1 $ESC[27m"; }
    get_cursor_row()   { IFS=';' read -sdR -p $'\E[6n' ROW COL; echo ${ROW#*[}; }
    key_input()        { read -s -n3 key 2>/dev/null >&2
                         if [[ $key = $ESC[A ]]; then echo up;    fi
                         if [[ $key = $ESC[B ]]; then echo down;  fi
                         if [[ $key = ""     ]]; then echo enter; fi; }

    # initially print empty new lines (scroll down if at bottom of screen)
    for opt; do printf "\n"; done

    # determine current screen position for overwriting the options
    local lastrow=`get_cursor_row`
    local startrow=$(($lastrow - $#))

    # ensure cursor and input echoing back on upon a ctrl+c during read -s
    trap "cursor_blink_on; stty echo; printf '\n'; exit" 2
    cursor_blink_off

    local selected=0
    while true; do
        # print options by overwriting the last lines
        local idx=0
        for opt; do
            cursor_to $(($startrow + $idx))
            if [ $idx -eq $selected ]; then
                print_selected "$opt"
            else
                print_option "$opt"
            fi
            ((idx++))
        done

        print_foot

        # user key control
        case `key_input` in
            enter) break;;
            up)    ((selected--));
                   if [ $selected -lt 0 ]; then selected=$(($# - 1)); fi;;
            down)  ((selected++));
                   if [ $selected -ge $# ]; then selected=0; fi;;
        esac
    done

    # cursor position back to normal
    cursor_to $lastrow
    printf "\n"
    cursor_blink_on

    return $selected
}

function select_opt {
    select_option "$@" 1>&2
    local result=$?
    echo $result
    return $result
}

# WORKING-DIR INITIALIZATION
function initialize {
    clear_view
    exit_loop='n'
    while [[ "${exit_loop}" == 'n' ]]
    do
        wd_initialization_menu

        # back to main menu
        if [[ "${exit_loop}" == 'exit' ]]
        then
            clear
            break
        fi
    done

    if [[ "${exit_loop}" != 'exit' ]]
    then
        # initializes the working-directory
        clear_view
        echo "
        Initializing the working directory..."
        sleep 1
        docker run --rm -v ${working_dir}:${working_dir} \
            -u "$(id -u)":"$(id -g)" \
            singgroup/my-brain-seq init_working_dir.sh ${working_dir} #FIXME

        # creates the temp file with the path to the new working_directory
        echo "WD=${working_dir}" > "${temp_file}"
        temp_file_exists='TRUE'
        update_wd

        enter_to_continue
    fi
}

# RUN THE ANALYSIS
function run_analysis {
    clear_view
    echo
    exit_loop='n'
    while [[ "${exit_loop}" == 'n' ]]
    do
        exit_loop='n'
        
        # if the working dir was set 
        if [[ "${temp_file_exists}" == 'TRUE' ]]
        then
            clear_view
            echo
            echo
            echo -ne "$(greenprint 'RUNNING MYBRAIN-SEQ ANALYSIS:')"
            echo
            echo "Run the analysis in the current working directory?"
            echo
            
            # display options
            exit_options
            
            # get compi parameters path or remove saved working dir
            if [[ "${exit_loop}" == 'y' ]]; then
                compi_parameters="${working_dir}/input/compi.parameters"
            elif [[ "${exit_loop}" == 'n' ]]; then
                temp_file_exists='FALSE'
            fi
        fi

        # if the working dir was not set
        if [[ "${temp_file_exists}" == 'FALSE' ]]
        then
            clear_view
            echo
            echo
            echo -ne "$(greenprint 'RUNNING MYBRAIN-SEQ ANALYSIS:')"
            echo
            echo "Please, paste the full path to the \"compi.parameters\" file."
            echo "e.g.: /home/user/analysis_1/input/compi.parameters"
            echo
            echo "Path to \"compi.parameters\" file:"
            echo -ne "    > "
            read compi_parameters

            # display options
            exit_options

            # get working dir from compi.parameters and save it
            if [[ "${exit_loop}" == 'y' ]]; then
                working_dir = "$(dirname $(dirname ${compi_parameters}))"
                echo "WD=${working_dir}" > "${temp_file}"
                temp_file_exists='TRUE'
                update_wd
            fi
        fi

        # back to main menu
        if [[ "${exit_loop}" == 'exit' ]]; then
            clear
            break
        fi
    done

    # run the analysis
    if [[ "${exit_loop}" != 'exit' ]]
    then
        clear_view
        echo "
        Running myBrain-Seq analysis..."
        sleep 1
        docker run --rm -it \
            -v '/var/run/docker.sock':'/var/run/docker.sock' \
            -v ${working_dir}:${working_dir} \
            -v ${compi_parameters}:${compi_parameters} \
            singgroup/my-brain-seq run ${compi_parameters}

        enter_to_continue        
    fi
}

update_wd

function mainmenu {
    clear_view
    echo
    echo -ne "$(magentaprint 'MAIN MENU')"
    echo
    case `select_opt "Initialize the working-directory" "Run myBrain-Seq analysis" "Exit"` in
    0)
        initialize
        mainmenu
        ;;
    1) 
        run_analysis
        mainmenu
        ;;
    2)
        fn_bye
        ;;
    esac
}

mainmenu

#command to run this script
#docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock -v /tmp:/tmp --entrypoint /visual_console/visual_console.sh singgroup/my-brain-seq