#! /usr/bin/env bash
# (c) Konstantin Riege
trap '
	sleep 1
	pids=($(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -KILL "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null
	echo -e "\r "
' EXIT
trap 'die "killed by sigint or sigterm"' INT TERM

die() {
	echo -ne "\e[0;31m"
	echo -e "\r:ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" == "bash" ]] || die "loading library requieres bash"
([[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]) || die "requieres bash version 4.4 or above"

shopt -s extglob
for f in "$(dirname $(readlink -e $0))"/lib/*.sh; do
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

# solve duplicate entries by reverting order of writeout: function >> $LOG 2> >(tee -a $LOG >&2) to function 2> >(tee -a $LOG >&2) >> $LOG
# slove sigpipe due to proken pipe by sigint by protecting tee's process with either trap '' INT or use tee -i to ignore termination signals
# script's trap is called in scope of function - i.e. using its file descriptors. to print trap messages as usual to fd 1 and 2 instead,
# either duplicate the main shell's original stdout and stderr for future use via exec 3>&1 4>&2 and do die(){ echo error } 1>&3 2>&4 or
# duplicate stderr and stdout of the function to alternative file descriptors. thus trap messages echoed in scope of the function goes to all file descriptors
exec 3>> $LOG
exec 4> >(tee -ai $LOG >&2)
for i in "${INSTALL[@]}"; do
	compile::$i -i $INSDIR -t $THREADS 2>&4 >&3 || die
done
exec 3>&-
exec 4>&-

commander::print "success" >> $LOG
exit 0
