#! /usr/bin/env bash
# (c) Konstantin Riege

source $(dirname $(readlink -e $0))/activate.sh -c false || exit 1

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0
BASHBONE_ERROR="parameterization issue"

options::parse "$@"

BASHBONE_ERROR="mandatory parameter -i missing"
[[ $INSTALL ]]
BASHBONE_ERROR="mandatory parameter -d missing"
[[ $INSDIR ]]

BASHBONE_ERROR="cannot access $INSDIR"
mkdir -p $INSDIR
INSDIR="$(readlink -e "$INSDIR")"

[[ $LOG ]] || LOG="$INSDIR/install.log"
BASHBONE_ERROR="cannot access $LOG"
progress::log -v $VERBOSITY -o "$LOG"

unset BASHBONE_ERROR
commander::printinfo "installation started. please be patient." >> "$LOG"
for i in "${INSTALL[@]}"; do
	BASHBONE_ERROR="compilation of $i failed"
	progress::observe -v $VERBOSITY -o "$LOG" -f compile::$i -i "$INSDIR" -t $THREADS
done

unset BASHBONE_ERROR
commander::printinfo "success" >> "$LOG"

exit 0
