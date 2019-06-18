#! /usr/bin/env bash
# (c) Konstantin Riege
shopt -s extglob
trap 'die' INT TERM
#trap 'kill -PIPE 0' EXIT # kills parental processes as well - shlog conflict
#trap 'kill -PIPE -- -$$' EXIT # kill all childs - works only if $$ is process group leader
trap 'kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT # parse pstree
# AVOID DOUBLE FORKS -> run(){cmd &}; run & -> i.e. cmd gets new process group and cannot be killed

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
bash --version | head -1 | cut -d ' ' -f 4 | cut -d '.' -f 1-2 | awk '$0<4.4{exit 1}' || die "requieres bash version 4.4 or above"

for f in "$(readlink -e $(dirname $0))"/lib/*.sh; do
	source $f
done

# overload fuctions

die() {
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

compile::src() {
	commander::print "installing bashbone"
	{	mkdir -p "$INSDIR/bashbone-$VERSION" && \
		cp -r "$(readlink -e $(dirname $0))"/* "$INSDIR/bashbone-$VERSION" && \
		mkdir -p "$INSDIR/latest" && \
		ln -sfn "$INSDIR/bashbone-$VERSION" "$INSDIR/latest/bashbone"
	} || return 1
	return 0
}

compile::all(){
	compile::src || return 1
	return 0
}

THREADS=$(cat /proc/cpuinfo | grep -cF processor)
VERBOSITY=0

options::parse "$@" || die "parameterization issue"

[[ ! $INSTALL ]] && die "mandatory parameter -i missing"
[[ ! $INSDIR ]] && die "mandatory parameter -d missing"
mkdir -p $INSDIR || die "cannot access $INSDIR"
INSDIR=$(readlink -e $INSDIR)
[[ ! $LOG ]] && LOG=$INSDIR/install.log
echo '' > $LOG || die "cannot access $LOG"

commander::print "installation started. please be patient."
progress::log -v $VERBOSITY -o $LOG

for i in ${INSTALL[@]}; do # do not quote!! mapfile appends newline to last element
	compile::$i -i $INSDIR &>> $LOG || die 
done

commander::print "success"
exit 0
