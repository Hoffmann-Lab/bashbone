#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" != "bash" ]] && {
	echo "loading library requieres bash"
} || {
	[[ "${BASH_SOURCE[0]}" == "${0}" ]] && {
		echo "script needs to be sourced"
		echo "do: source $(readlink -e $0)"
		exit 1
	} || {
		[[ ! $OSTYPE =~ linux ]] && echo "unsupported operating system" || {
			([[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]) && {
				init_insdir_bashbone=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
				unset OPTIND
				while getopts i: ARG; do
					case $ARG in
						i) init_insdir_bashbone="$OPTARG";;
					esac
				done
				IFS=$'\n'
				for f in $(readlink -e "$init_insdir_bashbone")/latest/bashbone/lib/*.sh; do
					source "$f"
				done && {
					unset IFS
					configure::environment -i $init_insdir_bashbone
				} || {
					unset IFS
					echo "install directory cannot be found"
					echo "please use parameter -i <path/to/install/dir>"
				}
			} || echo "requieres bash version 4.4 or above"
		}
	}
}
