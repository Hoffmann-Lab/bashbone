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
				insdir_bashbone=$(dirname $(readlink -e ${BASH_SOURCE[0]}))
				insdir_tools_bashbone=$(dirname $insdir_bashbone)
				activate_conda_bashbone=true
				error=false
				unset OPTIND
				while getopts :i:c: ARG; do
					case $ARG in
						i) insdir_tools_bashbone="$OPTARG";;
						c) activate_conda_bashbone="$OPTARG";;
						:) echo "argument missing"; error=true;;
					esac
				done
				$error || {
					IFS=$'\n'
					for f in "$insdir_bashbone/lib/"*.sh; do
						source "$f"
					done && {
						unset IFS
						configure::environment -i $insdir_tools_bashbone -b $insdir_bashbone -c $activate_conda_bashbone || {
							unset IFS
							echo ":ERROR: bashbone $version activation failed" >&2
							echo ":ERROR: tried $insdir_tools_bashbone but"
							echo ":ERROR: tools installation directory was not found to setup environment" >&2
							echo ":ERROR: please use or check parameter -i <path/to/install/dir>" >&2
							echo ":ERROR: OR" >&2
							echo ":ERROR: disable tools and conda activation by using parameter -c false" >&2
							echo ":ERROR: the latter will heavily limitate bashbone functionality" >&2
						}
						# add kill to not raise error in case someone adds a parmeter like -h or --help
						alias bashbone="declare -f | grep -P '::[^_].+\(\)' | grep -vF compile:: | sort -V | sed -r 's/\s+\(\)\s+$//'; kill -SIGINT $$"
					} || {
						echo ":ERROR: bashbone $version activation failed" >&2
						echo ":ERROR: unexpected error in source code - please contact developer" >&2
					}
				}
			} || echo "requieres bash version 4.4 or above"
		}
	}
} || echo "loading library requieres bash"
