#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(readlink -e "$0")")/activate.sh" -c false -a "$@" || exit 1

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0
BASHBONE_ERROR="parameterization issue"

options::parse "$@"

BASHBONE_ERROR="mandatory parameter -i missing"
[[ $INSTALL ]]
BASHBONE_ERROR="mandatory parameter -d missing"
[[ $INSDIR ]]

BASHBONE_ERROR="cannot access $INSDIR"
mkdir -p "$INSDIR"
INSDIR="$(readlink -e "$INSDIR")"

[[ $LOG ]] || LOG="$INSDIR/install.log"
BASHBONE_ERROR="cannot access $LOG"

commander::printinfo "installation started. please be patient." | tee -i "$LOG"
source "$(dirname "$(readlink -e "$0")")/activate.sh" -c false -i "$INSDIR"
for i in "${INSTALL[@]}"; do
	BASHBONE_ERROR="compilation of $i failed"
	progress::log -v $VERBOSITY -o "$LOG" -f compile::$i -i "$INSDIR" -t $THREADS -g ${USECONFIG:-false}
done

unset BASHBONE_ERROR
commander::printinfo "success" | tee -ia "$LOG"

exit 0
