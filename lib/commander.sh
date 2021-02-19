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
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-a <cmds>      | array of
			-s <seperator> | string for
			-o <outfile>   | stdout redirection to
			-c <cmd|fd3..> | ALWAYS LAST OPTION
			                 command line string(s) and or
			                 file descriptor(s) starting from 3
			example:
			${FUNCNAME[1]} -a cmds -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			    perl -sl - -y=$x <<< '
			        print "$x";
			        print "\$y";
			    '
			CMD
			    awk '{print $0}'
			CMD
		EOF
		return 1
	}

	local OPTIND arg mandatory sep='|' suffix='' fd tmp
	declare -n _cmds_makecmd # be very careful with circular name reference
	declare -a mapdata cmd_makecmd # be very careful with references name space
	while getopts ':a:o:s:c' arg; do
		case $arg in
			a)	mandatory=1; _cmds_makecmd=$OPTARG;;
			s)	sep=$(echo -e "$OPTARG");; # echo -e to make e.g. '\t' possible
			o)	suffix=' > '"$OPTARG";;
			c)	[[ ! $mandatory ]] && _usage
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
			*)	_usage;;
		esac
	done

	_usage
}

commander::printcmd(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-a <cmds> | array of
		EOF
		return 1
	}

	local OPTIND arg
	declare -n _cmds_printcmd # be very careful with circular name reference
	while getopts 'a:' arg; do
		case $arg in
			a)	_cmds_printcmd=$OPTARG
				[[ "${#_cmds_printcmd[@]}" -gt 0 ]] && printf ':CMD: %s\n' "${_cmds_printcmd[@]}"
				return 0;;
			*)	_usage;;
		esac
	done

	_usage
}

commander::runcmd(){
	local tmpdir
	_cleanup::commander::runcmd(){
		rm -rf "$tmpdir"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v             | verbose on
			-b             | benchmark on
			-c <env>       | run with conda
			-t <instances> | number of parallel
			-a <cmds>      | ALWAYS LAST OPTION
			                 array of
			example:
			${FUNCNAME[1]} -v -b -t 1 -a cmd
		EOF
		return 1
	}

	local OPTIND arg threads=1 verbose=false benchmark=false cenv
	declare -n _cmds_runcmd # be very careful with circular name reference
	while getopts 'vbt:c:a:' arg; do
		case $arg in
			t)	threads=$OPTARG;;
			v)	verbose=true;;
			b)	benchmark=true;;
			c)	cenv=$OPTARG;;
			a)	_cmds_runcmd=$OPTARG
				[[ $_cmds_runcmd ]] || return 0
				$verbose && {
					commander::printinfo "running commands of array ${!_cmds_runcmd}"
					commander::printcmd -a _cmds_runcmd
				}
				local i sh
				tmpdir=$(mktemp -d -p /dev/shm jobs.XXXXXXXXXX)
				# better write to file to avoid xargs argument too long error due to -I {}
				# old: printf '%s\0' "${_cmds_runcmd[@]}" | xargs -0 -P $threads -I {} bash -c {}
				# use a subshell with its own trap, destructed if subshell terminates
				if $benchmark; then
					for i in "${!_cmds_runcmd[@]}"; do
						sh="$(mktemp -p "$tmpdir" job.XXXXXXXXXX.sh)"
						echo "#!/usr/bin/env bash" > "$sh"
						#echo "set -o pipefail" >> "$sh"
						#echo "trap 'e=\$?; if [[ \$e -ne 141 ]]; then exit \$e; fi' ERR" >> "$sh" # do not use set -e, to ignore sigpipe caused by e.g. samtools view | head
						#[[ $cenv ]] && echo "source $CONDA_PREFIX/bin/activate $cenv" >> "$sh"
						if [[ $cenv ]]; then
							echo "source '$BASHBONE_DIR/activate.sh' -c true" >> "$sh"
							echo "conda activate $cenv" >> "$sh"
						else
							echo "source '$BASHBONE_DIR/activate.sh' -c false" >> "$sh"
						fi
						printf '%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
						echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
						echo "$sh"
					done | $(command -v time) -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -P $threads -I {} bash {}
					# time may be a shell keyword, so use full path
					# due to set -E based error tracing, command time may leads to *** longjmp causes uninitialized stack frame ***: bash terminated
					# workaround: use which or command -v
				else
					for i in "${!_cmds_runcmd[@]}"; do
						sh="$(mktemp -p "$tmpdir" job.XXXXXXXXXX.sh)"
						echo "#!/usr/bin/env bash" > "$sh"
						# echo "set -o pipefail" >> "$sh"
						# echo "trap 'e=\$?; if [[ \$e -ne 141 ]]; then exit \$e; fi' ERR" >> "$sh"
						# [[ $cenv ]] && echo "source $CONDA_PREFIX/bin/activate $cenv" >> "$sh"
						if [[ $cenv ]]; then
							echo "source '$BASHBONE_DIR/activate.sh' -c true" >> "$sh"
							echo "conda activate $cenv" >> "$sh"
						else
							echo "source '$BASHBONE_DIR/activate.sh' -c false" >> "$sh"
						fi
						printf '%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
						echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
						echo "$sh"
					done | xargs -P $threads -I {} bash {}
				fi
				return 0
			;;
			*)	_usage;;
		esac
	done

	_usage
}

commander::qsubcmd(){
	local jobname
	_cleanup::commander::qsubcmd(){
		[[ "$dowait" == "y" ]] && qdel -f "$jobname" &> /dev/null || true
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v             | verbose on
			-b             | benchmark on
			-w             | do wait for jobs and receive a non-null exit code if a single job fails
			-n <name>      | prefix of logs and jobs to wait for - should be unique
			-c <env>       | run with conda
			-l <complex>   | sge digestable list of consumables as key="value" pairs (see qconf -sc or -mc)
			-o <path>      | shared among machines for scripts, logs and exit codes
			-r             | override existing logs
			-i <instances> | number of parallel instances
			-q <queue>     | name of sge queue
			-p <env>       | name of parallel sge environment
			-t <threads>   | to be allocated per instance in parallel environment
			-a <cmds>      | ALWAYS LAST OPTION
			                 array of
			example:
			${FUNCNAME[1]} -v -l hostname="!bcl102&!bcl103" -l mem_free="50G" -c segemehl -p threads -t 4 -w -o ~/logs -a cmd
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 instances verbose=false benchmark=false dowait="n" override=false cenv penv queue logdir complex params
	declare -n _cmds_qsubcmd # be very careful with circular name reference
	declare -a mapdata complexes
	while getopts 'vbwrt:i:o:l:p:q:c:n:a:' arg; do
		case $arg in
			v)	verbose=true;;
			b)	benchmark=true;;
			w)	dowait="y";;
			r)	override=true;;
			c)	cenv=$OPTARG;;
			t)	threads=$OPTARG;;
			i)	instances=$OPTARG;;
			o)	((++mandatory)); logdir="$OPTARG"; mkdir -p "$logdir";;
			l)	complexes+=("-l $OPTARG");;
			p)	((++mandatory)); penv="-pe $OPTARG";;
			q)	((++mandatory)); queue="-q $OPTARG";;
			n)	jobname="$OPTARG";;
			a)	_cmds_qsubcmd=$OPTARG
				[[ $mandatory -lt 2 ]] && _usage
				[[ $_cmds_qsubcmd ]] || return 0

				$verbose && {
					commander::printinfo "running commands of array ${!_cmds_qsubcmd}"
					commander::printcmd -a _cmds_qsubcmd
				}

				[[ $penv ]] && params="$penv $threads" || params="$queue"
				[[ $jobname ]] || jobname=$(mktemp -u -p "$logdir" XXXXXXXXXX)
				local ex="$logdir/exitcodes.$jobname"
				local log="$logdir/job.$jobname.\$TASK_ID.log" # not SGE_TASK_ID

				local i id sh
				for i in "${!_cmds_qsubcmd[@]}"; do
					id=$((i+1))
					if $override; then
						rm -f "$logdir/job.$jobname.$id.log" "$ex"
					fi
					sh="$logdir/job.$jobname.$id.sh"
					cat <<- EOF > "$sh"
						#!/usr/bin/env bash
						exit::$jobname.$id(){
							echo "$jobname.$id exited with exit code \$1" >> "$ex"
						}
					EOF
					if [[ $cenv ]]; then
						echo "source '$BASHBONE_DIR/activate.sh' -c true -x exit::$jobname.$id" >> "$sh"
						echo "conda activate $cenv" >> "$sh"
					else
						echo "source '$BASHBONE_DIR/activate.sh' -c false -x exit::$jobname.$id" >> "$sh"
					fi
					printf '%s\n' "${_cmds_qsubcmd[$i]}" >> "$sh"
					chmod 755 "$sh"
				done

				[[ $instances ]] || instances=$id
				if $benchmark; then
					TIMEFORMAT=':BENCHMARK: runtime %3lR [hours][minutes]seconds' # different to /usr/bin/time, bash builtin time can handle: time echo "sleep 2" | bash
					time echo "$logdir/job.$jobname.\$SGE_TASK_ID.sh" | qsub -sync $dowait $params ${complexes[@]} -t 1-$id -tc $instances -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -o "$log" -j y -N $jobname |& sed -u -E '/exited/!d; s/Job [0-9]+\.(.+)\./job.'$jobname'.\1/;t;s/Job [0-9]+ (.+)\./job.'$jobname'.1 \1/'
					# in case of job exit code > 0, leads to *** longjmp causes uninitialized stack frame ***: bash terminated
				else
					echo "$logdir/job.$jobname.\$SGE_TASK_ID.sh" | qsub -sync $dowait $params ${complexes[@]} -t 1-$id -tc $instances -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -o "$log" -j y -N $jobname |& sed -u -E '/exited/!d; s/Job [0-9]+\.(.+)\./job.'$jobname'.\1/;t;s/Job [0-9]+ (.+)\./job.'$jobname'.1 \1/'
				fi
				return $?
			;;
			*)	_usage;;
		esac
	done

	_usage
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
