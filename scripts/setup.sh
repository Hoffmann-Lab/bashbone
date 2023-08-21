#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(dirname "$(readlink -e "$0")")")/activate.sh" -c false -r true -x cleanup -a "$@" || exit 1

cleanup() {
	[[ -e "$LOG" ]] && {
		echo "date: $(date)" | tee -ia "$LOG"
		[[ $1 -eq 0 ]] && echo "success" | tee -ia "$LOG" || echo "failed" | tee -ia "$LOG"
	}
}

CMD="$(basename "$0") $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
VERBOSITY=0
BASHBONE_ERROR="parameterization issue"

options::parse "$@"

BASHBONE_ERROR="mandatory parameter -i missing"
[[ $INSTALL ]]
BASHBONE_ERROR="mandatory parameter -d missing"
[[ $INSDIR ]]

BASHBONE_ERROR="cannot access $INSDIR"
mkdir -p "$INSDIR"
BASHBONE_TOOLSDIR="$(readlink -e "$INSDIR")"

[[ $LOG ]] || LOG="$BASHBONE_TOOLSDIR/install.log"
BASHBONE_ERROR="cannot access $LOG"

commander::printinfo "bashbone $BASHBONE_VERSION installation started with command: $CMD" | tee -i "$LOG"
commander::printinfo "date: $(date)" | tee -ia "$LOG"

for i in "${INSTALL[@]}"; do
	BASHBONE_ERROR="compilation of $i failed"
	progress::log -v $VERBOSITY -o "$LOG" -f compile::$i -i "$BASHBONE_TOOLSDIR" -t $THREADS -g ${USECONFIG:-false}
done
unset BASHBONE_ERROR

exit 0
