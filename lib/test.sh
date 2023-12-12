#! /usr/bin/env bash
# (c) Konstantin Riege

function test::makecmd(){
	declare -a cmds
	for i in {1..20}; do
		commander::makecmd -m -a cmds -s ";" -c {COMMANDER[0]}<<- CMD
			sleep $((RANDOM%10+1))
			echo foo$i
			sleep 2
			echo bar$i
		CMD
	done
	commander::printcmd -a cmds
	declare -a runcmds
	commander::makecmd -a runcmds -v cmds -c <<-'CMD'
		commander::runcmd -v -i 10 -a cmds
	CMD
	commander::runcmd -v -i 1 -a runcmds
}

function test::test(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-h | help
		EOF
		return 1
	}

	local OPTIND arg tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_coexpression _idfiles_coexpression
	while getopts 'h' arg; do
		case $arg in
			h)	_usage || true; return 0;;
			*)	_usage;;
		esac
	done
	commander::printinfo "running test::test"
	echo "echo test::test cleanup" >> "$BASHBONE_CLEANUP"
	declare -a cmd=('sleep 2')
	commander::runcmd -v -b -i 1 -a cmd
}

function test::start(){
	echo test::start under $BASHPID >&2
	cat <(test::inner; echo mut not show test::start >&2)
	# test::inner
	echo mut not show test::start >&2
}

function test::inner(){
	echo test::inner under $BASHPID >&2
	# cat <(test::error; echo mut not show test::inner >&2)
	test::error
	echo mut not show test::inner >&2
}

function test::error(){
	echo test::error under $BASHPID >&2
	# cat sdfdsfdf
	# (cat sdfdsfdf; echo mut not show test::error >&2)
	x=$(cat sdfdsfdf; echo mut not show test::error >&2)
	# cat <(cat sdfdsfdf; echo mut not show test::error >&2)
	# { cat sdfdsfdf; echo mut not show test::error >&2; } &
	# wait $!
	# dont capture tests! # if [[ $(cat sdfdsfdf) == foo ]]; then echo mut not show test::error >&2; fi
	# for i in $(cat sdfdsfdf; exit 1); do sleep 10; echo mut not show test::error loop >&2; done
	sleep 10
	echo mut not show test::error >&2
}
