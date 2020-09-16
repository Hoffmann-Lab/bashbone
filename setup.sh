#! /usr/bin/env bash
# (c) Konstantin Riege
set -o pipefail

# defines INSDIR
source $(dirname $(readlink -e $0))/activate.sh -c false || exit 1

unset ERROR
trap 'configure::exit -p $$' EXIT
trap 'ERROR="killed"' INT TERM
trap 'configure::err -x $? -s "$0" -l $LINENO -e "$ERROR" -c "$BASH_COMMAND"; exit $?' ERR

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0
ERROR="parameterization issue"
options::parse "$@"

ERROR="mandatory parameter -i missing"
[[ $INSTALL ]]

ERROR="cannot access $INSDIR"
mkdir -p $INSDIR
INSDIR=$(readlink -e $INSDIR)

[[ $LOG ]] || LOG=$INSDIR/install.log
ERROR="cannot access $LOG"
printf '' > $LOG

unset ERROR
progress::log -v $VERBOSITY -o $LOG
commander::printinfo "installation started. please be patient." >> $LOG

for i in "${INSTALL[@]}"; do
	ERROR="compilation of $i failed"
	compile::$i -i $INSDIR -t $THREADS 2> >(tee -ai $LOG >&2) >> $LOG
done

unset ERROR
commander::printinfo "success" >> $LOG

exit 0
