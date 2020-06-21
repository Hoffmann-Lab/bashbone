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
				activate_conda_bashbone=true
				error=false
				unset OPTIND
				while getopts :i:c: ARG; do
					case $ARG in
						i) activate_insdir_bashbone="$OPTARG";;
						c) activate_conda_bashbone="$OPTARG";;
						:) echo "argument missing"; error=true;;
					esac
				done
				$error || {
					IFS=$'\n'
					for f in $(readlink -e "$activate_insdir_bashbone")/latest/bashbone/lib/*.sh; do
						source "$f"
					done && {
						unset IFS
						configure::environment -i $activate_insdir_bashbone -c $activate_conda_bashbone
					} || {
						unset IFS
						echo "install directory cannot be found to activate conda"
						echo "please use parameter -i <path/to/install/dir>"
						echo "or"
						echo "disable conda activation by parameter -c false (not recommended)"
						echo "ATTENTION: disabeling conda activation heavily limitates bashbone functionality"
					}
					# add kill to not raise error in case someone adds a parmeter like -h or --help
					alias bashbone="declare -f | grep -P '::[^_].+\(\)' | grep -vF compile:: | sort -V | sed -r 's/\s+\(\)\s+$//'; kill -SIGINT $$"
				}
			} || echo "requieres bash version 4.4 or above"
		}
	}
} || echo "loading library requieres bash"
