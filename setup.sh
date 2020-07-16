#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; echo -e "\r "; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT

die() {
	echo -ne "\e[0;31m"
	echo "\r:ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" == "bash" ]] || die "loading library requieres bash"
([[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]) || die "requieres bash version 4.4 or above"

for f in "$(readlink -e $(dirname $0))"/lib/*.sh; do
	source $f || die "unexpected error in source code - please contact developer"
done

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0

options::parse "$@" || die "parameterization issue"

[[ ! $INSTALL ]] && die "mandatory parameter -i missing"
[[ ! $INSDIR ]] && die "mandatory parameter -d missing"

for f in "${TOSOURCE[@]}"; do
	source $f || die "unexpected error in source code - please contact developer"
done

mkdir -p $INSDIR || die "cannot access $INSDIR"
INSDIR=$(readlink -e $INSDIR)
[[ ! $LOG ]] && LOG=$INSDIR/install.log

echo > $LOG || die "cannot access $LOG"
progress::log -v $VERBOSITY -o $LOG
commander::print "installation started. please be patient." >> $LOG

for i in "${INSTALL[@]}"; do
	compile::$i -i $INSDIR -t $THREADS >> $LOG 2> >(tee -a $LOG >&2) || die 
done

commander::print "success" >> $LOG
exit 0
