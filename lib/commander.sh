#!/usr/bin/env bash
# (c) Konstantin Riege

declare -a COMMANDER

commander::print(){
	[[ $* ]] && echo ":INFO: $*"
	local fd
	declare -a mapdata
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t mapdata
		printf '%s\n' "${mapdata[@]}"
		# while read -u $fd -r tmp; do
		# 	echo ":INFO: $tmp"
		# done
	done
	COMMANDER=()
	return 0
}

commander::printinfo(){
	[[ $* ]] && echo ":INFO: $*"
	local fd
	declare -a mapdata
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t mapdata
		printf ':INFO: %s\n' "${mapdata[@]}"
	done
	COMMANDER=()
	return 0
}

commander::warn(){
	[[ $* ]] && echo ":WARNING: $*"
	local fd
	declare -a mapdata
	for fd in "${COMMANDER[@]}"; do
		mapdata
		mapfile -u $fd -t mapdata
		printf ':WARNING: %s\n' "${mapdata[@]}"
	done
	COMMANDER=()
	return 0
}

commander::printerr(){
	[[ $* ]] && echo ":ERROR: $*" 1>&2
	local fd
	declare -a mapdata
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t mapdata
		printf ':ERROR: %s\n' "${mapdata[@]}" 1>&2
	done
	COMMANDER=()
	return 0
}

commander::makecmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-a <cmds>      | array of
			-s <seperator> | string for
			-o <outfile>   | stdout redirection to
			-c <cmd|fd3..> | ALWAYS LAST OPTION
			                 command line string(s) and or
			                 file descriptor(s) starting from 3
			example:
			$funcname -a cmds -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			    perl -sl - -y=$x <<< '
			        print "$x";
			        print "\$y";
			    '
			CMD
			    awk '{print $0}'
			CMD
		EOF
		return 0
	}

	local OPTIND arg mandatory sep='|' suffix='' fd tmp
	declare -n _cmds_makecmd # be very careful with circular name reference
	declare -a mapdata cmd_makecmd # be very careful with references name space
	while getopts ':a:o:s:c' arg; do
		case $arg in
			a)	mandatory=1; _cmds_makecmd=$OPTARG;;
			s)	sep=$(echo -e "${OPTARG:- }");; # echo -e to make e.g. '\t' possible
			o)	suffix=' > '"$OPTARG";;
			c)	[[ ! $mandatory ]] && { _usage; return 1; }
				shift $((OPTIND-1)) # remove '-a <cmd>' '-s <char> '-o <file>' '-c' from $@
				for fd in "${COMMANDER[@]}"; do
					mapfile -u $fd -t mapdata
					cmd_makecmd+=("${mapdata[*]}") # * instead of @ to concatenate a non-splitted sentence instead of single non-splitted words
					exec {fd}>&- # just to be safe
					# exec {fd}>&- to close fd in principle not necessary since heredoc is read only and handles closure after last EOF
				done
				COMMANDER=()
				tmp="${cmd_makecmd[*]/#/$sep }" # concatenate CMD* with seperator
				_cmds_makecmd+=("$* ${tmp/#$sep /}$suffix")
				return 0;;
			*)	_usage; return 1;;
		esac
	done

	_usage; return 1
}

commander::printcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-a <cmds> | array of
		EOF
		return 0
	}

	local OPTIND arg
	declare -n _cmds_printcmd # be very careful with circular name reference
	while getopts 'a:' arg; do
		case $arg in
			a)	_cmds_printcmd=$OPTARG
				[[ "${#_cmds_printcmd[@]}" -gt 0 ]] && printf ':CMD: %s\n' "${_cmds_printcmd[@]}"
				return 0;;
			*)	_usage; return 1;;
		esac
	done

	_usage; return 1
}

commander::runcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v           | verbose on
			-b           | benchmark on
			-c <env>     | run with conda
			-t <threads> | number of
			-a <cmds>    | ALWAYS LAST OPTION
			               array of
			example:
			$funcname -v -b -t 1 -a cmd
		EOF
		return 0
	}

	local OPTIND arg threads=1 verbose=false benchmark=false conda
	declare -n _cmds_runcmd # be very careful with circular name reference
	while getopts 'vbt:c:a:' arg; do
		case $arg in
			t)	threads=$OPTARG;;
			v)	verbose=true;;
			b)	benchmark=true;;
			c)	conda=$OPTARG;;
			a)	_cmds_runcmd=$OPTARG
				[[ $_cmds_runcmd ]] || return 0
				$verbose && {
					commander::printinfo "running commands of array ${!_cmds_runcmd}"
					commander::printcmd -a _cmds_runcmd
				}
				local i sh tmpdir
				# better write to file to avoid xargs argument too long error due to -I {}
				# old: printf '%s\0' "${_cmds_runcmd[@]}" | xargs -0 -P $threads -I {} bash -c {}
				# use a subshell with its own trap, destructed if subshell terminates
				$benchmark && {
					(	trap 'rm -rf "$tmpdir"' EXIT
						trap 'exit $?' ERR INT TERM
						set -e -o pipefail
						tmpdir=$(mktemp -d -p /dev/shm jobs.XXXXXXXXXX)
						for i in "${!_cmds_runcmd[@]}"; do
							sh="$(mktemp -p "$tmpdir" job.XXXXXXXXXX.sh)"
							echo "#!/usr/bin/env bash" > "$sh"
							echo "set -e -o pipefail" >> "$sh"
							[[ $conda ]] && echo "source $CONDA_PREFIX/bin/activate $conda" >> "$sh"
							printf '%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
							echo "$sh"
						done | command time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -P $threads -I {} bash {}
					)
				} || {
					(	trap 'rm -rf "$tmpdir"' EXIT
						trap 'exit $?' ERR INT TERM
						set -e -o pipefail
						tmpdir=$(mktemp -d -p /dev/shm jobs.XXXXXXXXXX)
						for i in "${!_cmds_runcmd[@]}"; do
							sh="$(mktemp -p "$tmpdir" job.XXXXXXXXXX.sh)"
							echo "#!/usr/bin/env bash" > "$sh"
							echo "set -e -o pipefail" >> "$sh"
							[[ $conda ]] && echo "source $CONDA_PREFIX/bin/activate $conda" >> "$sh"
							printf '%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
							echo "$sh"
						done | xargs -P $threads -I {} bash {}
					)
				}
				return $?;;
			*)	_usage; return 1;;
		esac
	done

	_usage; return 1
}

commander::qsubcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v           | verbose on
			-b           | benchmark on
			-c <env>     | run with conda
			-h <hosts>   | sge digestable list of
			-l <logfile> | sge nodes shared path to
			-t <threads> | number of
			-p <env>     | name of parallel
			-a <cmds>    | ALWAYS LAST OPTION
			               array of
			example:
			$funcname -v -h "!bcl102&!bcl103" -p envname -t 4 -a cmd
			$funcname -v -h "!bcl102&!bcl103" -a cmd
		EOF
		return 0
	}

	local OPTIND arg threads=1 verbose=false benchmark=false conda penv log hosts
	declare -n _cmds_qsubcmd # be very careful with circular name reference
	declare -a mapdata
	while getopts 'vbt:l:h:p:c:a:' arg; do
		case $arg in
			v)	verbose=true;;
			b)	benchmark=true;;
			c)	conda=$OPTARG;;
			t)	threads=$OPTARG;;
			l)	log="$OPTARG";;
			h)	hosts="-l h=$OPTARG";;
			p)	penv="-pe $OPTARG";;
			a)	_cmds_qsubcmd=$OPTARG
				[[ $_cmds_qsubcmd ]] || return 0
				$verbose && {
					commander::printinfo "running commands of array ${!_cmds_qsubcmd}"
					commander::printcmd -a _cmds_qsubcmd
				}
				[[ $penv ]] && penv+=" $threads"
				local i sh e tmpdir ex
				local jobname
				(	trap '[[ $jobname ]] && qdel "$jobname.*"; rm -rf "$tmpdir"; rm -f "$ex"' EXIT
					trap 'exit $?' ERR INT TERM
					set -e -o pipefail
					tmpdir=$(mktemp -d -p /dev/shm jobs.XXXXXXXXXX)
					jobname="$(basename "$tmpdir")"
					jobname="X${jobname#*.}" # ensure first character to be a letter
					[[ $log ]] && ex="$(dirname "$log")"/$jobname.exitcodes || ex="$tmpdir/exitcodes"
					for i in "${!_cmds_qsubcmd[@]}"; do
						sh="$(mktemp -p "$tmpdir" job.XXXXXXXXXX.$i.sh)"
						[[ ! $log ]] && log="${sh%.*}.out"

						echo "#!/usr/bin/env bash" > "$sh"
						echo "set -e -o pipefail" >> "$sh"
						[[ $conda ]] && echo "source $CONDA_PREFIX/bin/activate $conda" >> "$sh"
						printf '%s\n' "${_cmds_qsubcmd[$i]}" >> "$sh"
						echo "echo \$? >> '$ex'" >> "$sh"

						qsub $penv $hosts -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -e "$log" -o "$log" -N $jobname.$i "$sh" > /dev/null
					done
					$benchmark && {
						command time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" \
						qsub $penv $hosts -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -b y -sync y -e /dev/null -o /dev/null -hold_jid "$jobname.*" -N $jobname.wait true > /dev/null
						e=$?
					} || {
						qsub $penv $hosts -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -b y -sync y -e /dev/null -o /dev/null -hold_jid "$jobname.*" -N $jobname.wait true > /dev/null
						e=$?
					}
					unset jobname # do this for qdel trap handling, since check for $? -gt 0 may call qdel just because cmd failed
					[[ $log && -e "$ex" ]] && {
						mapfile -t mapdata < "$ex"
						exit $((${mapdata[@]/%/+}0))
					} || {
						exit $e
					}
				)
				return $?;;
			*)	_usage;	return 1;;
		esac
	done

	_usage; return 1
}

commander::_test(){
	local x=${1:-'hello world'} threads=${2:-1} i=2 s

	while read -u $((++i)) -r s 2> /dev/null; do
		echo $s
	done 3<<< 'foo' 4<<- EOF 5< <(echo baz)
		bar
	EOF

	declare -a cmd
	commander::makecmd -a cmd -c perl -sl - -x="'$x'" \<\<\< \''print "$x"'\' \| awk \''{print $0}'\'

	commander::makecmd -a cmd -c perl -sl - -x="'$x'" {COMMANDER[0]}<<- 'CMD'
		<<< '
			print "$x";
		' | awk '{print $0}'
	CMD

	commander::makecmd -a cmd -c perl -l - {COMMANDER[0]}<<- CMD
		<<< '
			print "$x";
		' | awk '{print \$0}'
	CMD

	commander::makecmd -a cmd -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		perl -l - <<< '
			print "$x";
		'
	CMD
		awk '{print $0}'
	CMD

	commander::makecmd -a cmd -s ' ' -c perl -s - -x="'$x'" {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		<<< '
			print "$x";
	CMD
			print " $x\n";
		'
	CMD

	commander::makecmd -a cmd -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -sE '
			say "$x";
		'
	CMD
		-- -x="$x"
	CMD

	commander::runcmd -v -b -t $threads -a cmd || commander::printerr "failed"

	commander::warn ${FUNCNAME[0]} EXPECTED OUTPUT to stderr
	commander::printerr {COMMANDER[0]}<<-OUT
		foo
		bar
		baz
		:INFO: running commands of array cmd
		:CMD: perl -sl - -x='hello world' <<< 'print "$x"' | awk '{print $0}'
		:CMD: perl -sl - -x='hello world' <<< ' print "$x"; ' | awk '{print $0}'
		:CMD: perl -l - <<< ' print "hello world"; ' | awk '{print $0}'
		:CMD: perl -l - <<< ' print "hello world"; ' | awk '{print $0}'
		:CMD: perl -s - -x='hello world' <<< ' print "$x";  print " hello world\n"; '
		:CMD: perl -sE ' say "$x"; '   -- -x="hello world"
		hello world
		hello world
		hello world
		hello world
		hello world hello world
		hello world
		:BENCHMARK: runtime 0:00.03 [hours:]minutes:seconds
		:BENCHMARK: memory 4328 Kbytes
	OUT
}


: <<-INFO
### idea simplified
makecmd(){
	declare -n arref=$1
	shift
	arref+=("$*\0")
}
runcmd(){
	threads=$1
	declare -n arref=$2
	echo -ne "${arref[@]}" | xargs -0 -P $threads -I {} bash -c {}
}
cmd=()
for i in {1..3}; do
	makecmd cmd echo $i \| perl -lane \''
		print "\"$_\"";
	'\'
	makecmd cmd << EOF
		echo $i | perl -lane '
		print "\$_";
	'
	EOF
done
runcmd 1 cmd

### file descriptors
exec 3> tmp
echo 'echo baz' >&3
exec 3>&-
exec 3< tmp
cat <&3 tmp

# here document implicitly calls 'exec 0< /dev/fd/0; echo .. >&0; exec 0<&- /dev/fd/0'
cat <<-CMD
	echo echo0
CMD
# equal to (0 = STDIN)
cat 0<<-CMD /dev/fd/0
	echo echo0
CMD
# using different descriptors !=1 (1 = STDOUT, 2 = STDERR)
cat 2<<-CMD /dev/fd/2
	echo echo2
CMD
cat 3<<-CMD /dev/fd/3
	echo echo3
CMD

# allow or deny all kinds of parameter expansion
cat 0<<-CMD 2<<-'CMD' 3<<-'CMD' /dev/fd/0 /dev/fd/2 /dev/fd/3
	echo $x |
CMD
	perl -lane '
		print "$_"
	' |
CMD
	awk '{
		print $0
	}'
CMD

INFO
