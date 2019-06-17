#! /usr/bin/env bash
# (c) Konstantin Riege
shopt -s extglob
trap 'die' INT TERM
#trap 'kill -PIPE 0' EXIT # kills parental processes as well - shlog conflict
#trap 'kill -PIPE -- -$$' EXIT # kill all childs - works only if $$ is process group leader
#trap 'kill -PIPE $(jobs -p)' EXIT # same as above
trap 'kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT # parse pstree
# AVOID DOUBLE FORKS -> run(){cmd &}; run & -> i.e. cmd gets new process group and cannot be killed

die() {
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

usage() {
	cat <<- EOF
		DESCRIPTION
		setup the RIP peak caller pipeline RIPPCHEN

		VERSION
		$version

		SYNOPSIS
		$(basename $0) -i [all|upgrade]

		OPTIONS
		-i | --install [all|upgrade] : install/upgrade into exported directory: export RIPPCHEN=/path/to/install/dir
		-t | --threads [value]       : threads - predicted default: $threads
		-l | --log [path]            : log file - default: RIPPCHEN/install.log
		-v | --verbose               : enable verbose mode
		-h | --help                  : prints this message

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 0
}

checkopt (){
	local arg=false
	case $1 in
		-h | --help) usage;;
		-v | --verbose) verbose=1; return 0;;
		-t | --threads) arg=true; threads=$2;;
		-l | --log) arg=true; log=$2;;
		-i | --install) arg=true; mapfile -d ',' -t install <<< $2;;
		-*) echo ":ERROR: illegal option $1"; return 1;; 
		*) echo ":ERROR: illegal option $2"; return 1;;
	esac
	$arg && {
		[[ ! $2 ]] && echo ":ERROR: argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	}
}

parse(){
	[[ $# -eq 0 ]] && usage
	[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && echo "illegal option $1" >&2 && return 1
	for i in $(seq 1 $#); do
		if [[ ${!i} =~ ^- ]]; then
			j=$((i+1))
			checkopt "${!i}" "${!j}" || return 1
		else 
			((++i))
		fi
	done

	return 0
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
bash --version | head -1 | cut -d ' ' -f 4 | cut -d '.' -f 1-2 | awk '$0<4.4{exit 1}' || die "requieres bash version 4.4 or above"
[[ ! $RIPPCHEN ]] && die "can not find installation. please do: export RIPPCHEN=/path/to/install/dir"
src=$(cd $(dirname $0) && echo $PWD)
threads=$(cat /proc/cpuinfo | grep -cF processor)
memory=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.8/1024)}')
insdir=$RIPPCHEN
verbose=''
log=''

source $src/lib/version.sh
parse "$@" || die "parameterization issue"

[[ ! $install ]] && die "mandatory parameter -i missing"
[[ $insdir =~ ^$src ]] && die "unsupported installation path $insdir. please do: export RIPPCHEN=/path/to/install/dir outside of source directory $src"
[[ ! $log ]] && log=$insdir/install.log

mkdir -p $(dirname $log) || die "can not create log file directory $(dirname $log)"
mkdir -p $insdir || die "can not create installation directory $insdir"

echo ":INFO: installation started. please be patient. check logs via: tail -f $log"
rm -f $log
touch $log
source $src/lib/compile.sh
if [[ $verbose ]]; then
	for i in ${install[@]}; do # do not quote!! mapfile appends newline to last element
		compile::$i 2>&1 | tee -a $log
		[[ ${PIPESTATUS[0]} -gt 0 ]] && die
	done
else
	source $src/lib/progress.sh
	progress::bar &
	progress::log1 $log &
	for i in ${install[@]}; do # do not quote!! mapfile appends newline to last element
		compile::$i &>> $log
		[[ $? -gt 0 ]] && die
	done
fi

echo ":INFO: success"
exit 0
