#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && {
	[[ "${BASH_SOURCE[0]}" == "${0}" ]] && {
		echo "script needs to be sourced"
		echo "do: source $(readlink -e $0)"
		exit 1
	} || {
		[[ ! $OSTYPE =~ linux ]] && echo "unsupported operating system" || {
			([[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]) && {
				activate_insdir_bashbone=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
				unset OPTIND
				while getopts :i: ARG; do
					case $ARG in
						i) activate_insdir_bashbone="$OPTARG";;
						:) echo "argument missing for option -i"; exit 1;;
					esac
				done
				IFS=$'\n'
				for f in $(readlink -e "$activate_insdir_bashbone")/latest/bashbone/lib/*.sh; do
					source "$f"
				done && {
					unset IFS
					configure::environment -i $activate_insdir_bashbone
				} || {
					unset IFS
					echo "install directory cannot be found"
					echo "please use parameter -i <path/to/install/dir>"
				}
				alias bashbone="declare -f | grep -P '::[^_].+\(\)' | grep -vF compile:: | sort -V | sed -r 's/\s+\(\)\s+$//'; kill -SIG $$;"
			} || echo "requieres bash version 4.4 or above"
		}
	}
} || echo "loading library requieres bash"
