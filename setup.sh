#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
shopt -s extglob

die() {
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
bash --version | head -1 | cut -d ' ' -f 4 | cut -d '.' -f 1-2 | awk '$0<4.4{exit 1}' || die "requieres bash version 4.4 or above"

for f in "$(readlink -e $(dirname $0))"/lib/*.sh; do
	source $f
done

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0

options::parse "$@" || die "parameterization issue"

[[ ! $INSTALL ]] && die "mandatory parameter -i missing"
[[ ! $INSDIR ]] && die "mandatory parameter -d missing"

for f in ${TOSOURCE[@]}; do
	source $f || die
done

mkdir -p $INSDIR || die "cannot access $INSDIR"
INSDIR=$(readlink -e $INSDIR)
[[ ! $LOG ]] && LOG=$INSDIR/install.log
touch $LOG && rm -f $LOG || die "cannot access $LOG"

commander::print "installation started. please be patient." | tee -a $LOG
progress::log -v $VERBOSITY -o $LOG

for i in ${INSTALL[@]}; do # do not quote!! mapfile appends newline to last element
	compile::$i -i $INSDIR -t $THREADS > >(tee -a $LOG) 2> >(tee -a $LOG >&2) || die 
done

commander::print "success" | tee -a $LOG
exit 0
