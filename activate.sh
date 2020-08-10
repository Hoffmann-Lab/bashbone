#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && {
	[[ "${BASH_SOURCE[0]}" == "$0" ]] && {
		echo ":ERROR: script needs to be sourced. do" >&2
		echo ":ERROR: source $(basename "$0")" >&2
		exit 1
	} || {
		[[ ! $OSTYPE =~ linux ]] && echo "unsupported operating system" || {
			if [[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]; then
				insdir_bashbone="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
				insdir_tools_bashbone="$(dirname "$insdir_bashbone")"
				unset OPTIND activate_conda_bashbone
				while getopts :i:c: arg; do
					case $arg in
						i) insdir_tools_bashbone="$OPTARG";;
						c) activate_conda_bashbone="$OPTARG";;
						:) echo ":ERROR: argument missing" >&2 ; return 1;;
					esac
				done
				IFS=$'\n'
				for f in "$insdir_bashbone/lib/"*.sh; do
					source "$f"
				done && {
					unset IFS
					configure::environment -i "$insdir_tools_bashbone" -b "$insdir_bashbone" -c ${activate_conda_bashbone:-false} && {
						[[ $activate_conda_bashbone ]] || {
							echo ":INFO: to activate conda environment do"
							echo ":INFO: source $(basename "${BASH_SOURCE[0]}") -c true"
						}
						INSDIR="$insdir_bashbone"
						bashbone() {
							_usage(){
								commander::print {COMMANDER[0]}<<- EOF
									Bashbone

									A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

									Usage:
									-h | help
									-l | list functions for users
									-e | list functions for experienced users
									-d | list functions for developers
									-a | list all, including non-standalone, functions
								EOF
								return 0
							}
							local OPTIND arg
							while getopts 'hleda' arg; do
								case $arg in
								h) _usage; return 0;;
								l) declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
								e) declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
								d) declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
								a) declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
								*) _usage; return 1;;
								esac
							done
							_usage
							return 0
						}
					} || {
						echo ":ERROR: bashbone environment setup failed! do" >&2
						echo ":ERROR: source $(basename "${BASH_SOURCE[0]}") -i <path/to/install/dir>" >&2
						echo ":ERROR: or to disable tools and conda activation do" >&2
						echo ":ERROR: source $(basename "${BASH_SOURCE[0]}") -c false" >&2
						return 1
					}
				} || {
					echo ":ERROR: bashbone activation failed" >&2
					echo ":ERROR: unexpected error in source code - please contact developer" >&2
					return 1
				}
			else
				echo ":ERROR: requieres bash version 4.4 or above" >&2
				return 1
			fi
		}
	}
} || {
	echo ":ERROR: loading library requieres bash" >&2
	return 1
}
