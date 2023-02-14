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

commander::makecmd2(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-a <cmds>      | array of
			-v <variable>  | env variable to pass to command
			-o <outfile>   | stdout redirection to
			-c <cmd|fd3..> | ALWAYS LAST OPTION
			                 command line string(s) and or
			                 file descriptor(s) starting from 3
			example 1:
			${FUNCNAME[1]} -a cmds -c perl -le \''print "foo"'\'

			example 2:
			x=1
			${FUNCNAME[1]} -v x -a cmds -c <<-'CMD'
			    perl -sl - -x=\$x <<< '
			        print "\$x";
			        print "\$y";
			    ' | awk '{print \$0}'
			CMD
		EOF
		return 1
	}

	local OPTIND arg mandatory vars suffix
	declare -a mapdata
	declare -n _cmds_makecmd
	while getopts ':v:a:o:c' arg; do
		case $arg in
			v)	declare -n _var_makecmd=$OPTARG
				vars+="$OPTARG=$(printf '%q;' "$_var_makecmd") " # to use multi-line variable assignemnts for job shell use '%q\n'
			;;
			a)	mandatory=1; _cmds_makecmd=$OPTARG;;
			o)	suffix=" > '$OPTARG'";;
			c)	[[ ! $mandatory ]] && _usage
				[[ -t 0 ]] || mapfile -t mapdata # to keep multi-line cmds use read -r -d '' cmd || true
				shift $((OPTIND-1))
				_cmds_makecmd+=("$vars${mapdata[*]}$*$suffix")
				return 0
			;;
			*) _usage;;
		esac
	done

	_usage
}

commander::makecmd(){
	_cleanup::commander::makecmd(){
		COMMANDER=()
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-a <cmds>      | array of
			-s <seperator> | string for (default: |)
			-o <outfile>   | stdout redirection to
			-c <cmd|fd3..> | ALWAYS LAST OPTION
			                 command line string(s) and or
			                 file descriptor(s) starting from 3
			example 1:
			${FUNCNAME[1]} -a cmds -c perl -le \''print "foo"'\'

			example 2:
			${FUNCNAME[1]} -a cmds -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			    perl -sl - -y=\$x <<< '
			        print "\$x";
			        print "\\\$y";
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
			o)	suffix=" > '$OPTARG'";;
			c)	[[ ! $mandatory ]] && _usage
				shift $((OPTIND-1)) # remove '-a <cmd>' '-s <char> '-o <file>' '-c' from $@
				for fd in "${COMMANDER[@]}"; do
					mapfile -u $fd -t mapdata
					cmd_makecmd+=("${mapdata[*]}") # * instead of @ to concatenate a non-splitted sentence instead of single non-splitted words
					exec {fd}>&- # just to be safe
					# exec {fd}>&- to close fd in principle not necessary since heredoc is read only and handles closure after last EOF
				done
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
	local tmpdir pgid logdir
	_cleanup::commander::runcmd(){
		# solution1: setsid xargs needs kill to pgid on SIG
		# pgrep -g $pgid  OR  ps -o pid= -g $pgid # | ps -o ppid= -g $pgid column -t | grep -Fxc $pgid  OR  pgrep -P $pgid
		#[[ $pgid ]] && { env kill -INT -- -$pgid; sleep 1; env kill -TERM -- -$pgid; wait $pgid $(ps -o pid= --ppid $pgid); } &> /dev/null || true
		[[ $logdir ]] || rm -rf "$tmpdir"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v             | verbose on
			-b             | benchmark on
			-c <env>       | run with conda
			-i <instances> | number of parallel
			-t <instances> | obsolete synonym for -i
			-s <idx[:idx]> | execute only jobs from cmds array starting from given index or range (default: 1)
			-n <name>      | optional. prefix of logs and jobs - should be unique
			-o <path>      | optional. for scripts, logs and exit codes
			-r             | optional. override existing logs
			-a <cmds>      | array of
			example:
			${FUNCNAME[1]} -v -b -i 2 -a cmd
			example2:
			${FUNCNAME[1]} -i 2 -n current -o ~/jobs -r -a cmd
		EOF
		return 1
	}

	local OPTIND arg mandatory instances=1 verbose=false benchmark=false cenv startid=1 stopid override=false jobname="current"
	declare -n _cmds_runcmd # be very careful with circular name reference
	while getopts 'vbt:i:c:a:s:o:n:r' arg; do
		case $arg in
			t)	instances=$OPTARG;; # obsolete, for compatibility
			i)	instances=$OPTARG;;
			v)	verbose=true;;
			b)	benchmark=true;;
			c)	cenv=$OPTARG;;
			a)	mandatory=1; _cmds_runcmd=$OPTARG;;
			s)	IFS=":" read -r startid stopid <<< "$OPTARG";;
			r)	override=true;;
			o)	logdir="$OPTARG"; mkdir -p "$logdir";;
			n)	jobname="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ ! $mandatory ]] && _usage
	[[ $_cmds_runcmd ]] || return 0
	declare -f "$BASHBONE_HOOKCMD" &> /dev/null && $BASHBONE_HOOKCMD _cmds_runcmd

	$verbose && commander::printinfo "running commands of array ${!_cmds_runcmd}"

	local i id sh ex log
	[[ $logdir ]] && tmpdir="$logdir" || tmpdir=$(mktemp -d -p /tmp jobs.XXXXXXXXXX)
	[[ $jobname ]] || jobname="$(basename "$(mktemp -u -p "$tmpdir" XXXXXXXXXX)")"
	echo $instances > "$tmpdir/instances.$jobname"
	ex="$tmpdir/exitcodes.$jobname"
	touch "$ex" # for runstat
	[[ $stopid ]] || stopid=${#_cmds_runcmd[@]}

	# better write to file to avoid xargs argument too long error due to -I {}
	# old: printf '%s\0' "${_cmds_runcmd[@]}" | xargs -0 -P $instances -I {} bash -c {}
	# upon error return 255 to prevent xargs to load further jobs
	# use exit function on PPID to kill all sibling processes executed by xargs
	declare -a scripts
	# necessary for runstat
	cat <<- EOF > "$tmpdir/info.$jobname"
		$startid
		$stopid
		$(date "+%a-%b-%d-%T-%Y")
		$instances
	EOF

	((startid--))
	((stopid--))
	for i in $(seq $startid $stopid); do
		id=$((i+1));
		sh="$tmpdir/job.$jobname.$id.sh"
		log="$tmpdir/job.$jobname.$id.log"
		$override && rm -f "$log" "$ex"

		echo '#!/usr/bin/env bash' > "$sh"
		# if [[ $logdir ]]; then
			#solution1: pgid=

			echo "exit::$jobname.$id(){" >> "$sh"
			echo "    echo \"$jobname.$id (\$((\$(ps -o ppid= -p \$\$ 2> /dev/null)))) exited with exit code \$1\" >> '$ex'" >> "$sh"
			echo '}' >> "$sh"
			echo "export BASHBONE_NOSETSID=true" >> "$sh"

			if [[ $cenv ]]; then
				echo "source '$BASHBONE_DIR/activate.sh' -c true -x exit::$jobname.$id -i '$BASHBONE_TOOLSDIR'" >> "$sh"
				echo "conda activate --no-stack $cenv" >> "$sh"
			else
				echo "source '$BASHBONE_DIR/activate.sh' -c false -x exit::$jobname.$id -i '$BASHBONE_TOOLSDIR'" >> "$sh"
			fi
		# else
		# 	if [[ $cenv ]]; then
		# 		echo "source '$BASHBONE_DIR/activate.sh' -c true -i '$BASHBONE_TOOLSDIR'" >> "$sh"
		# 		echo "conda activate --no-stack $cenv" >> "$sh"
		# 	else
		# 		echo "source '$BASHBONE_DIR/activate.sh' -c false -i '$BASHBONE_TOOLSDIR'" >> "$sh"
		# 	fi
		# fi

		# $verbose && echo 'tail -3 "$0" | head -1 | paste -d " " <(echo ":CMD:") -' >> "$sh"
		# echo '{\n%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
		# echo "} 2> >(tee -ia '$log' >&2) > >(tee -ia '$log') | cat" >> "$sh"

		$verbose && echo 'tail -2 "$0" | head -1 | paste -d " " <(echo ":CMD:") -' >> "$sh"
		echo "exec 1> >(tee -ai '$log')" >> "$sh"
		echo "exec 2> >(tee -ai '$log' >&2)" >> "$sh"
		echo "${_cmds_runcmd[$i]}" >> "$sh"
		echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
		# $verbose && echo 'tail -3 "$0" | head -1 | paste -d " " <(echo ":CMD:") -' >> "$sh"
		#printf '{\n%s\n} &\nwait $!\nexit 0\n' "${_cmds_runcmd[$i]}" >> "$sh"
		# run asynchronous and use wait to get rid of terminated messages. but this will print this for loop as terminated job command at wait below
		# thus, use INT signal, but attention: kill -INT is ignored in asynchronous commands with disabled job control due to a wierd POSIX requirement: set +m; bash -c 'trap "echo INT" INT; trap -p' & wait $!
		# workaround via env: set +m; env --default-signal=SIGINT,SIGQUIT bash -c 'trap "echo INT" INT; trap -p' & wait $!
		scripts+=("$(realpath -se "$sh")") # necessary for runstat
	done

	echo -e 'will cite' | parallel --citation &> /dev/null || true
	if $benchmark; then
		# solution1: setsid xargs, see also this _cleanup and configure::exit
		#printf '%q\n' "${scripts[@]}" | setsid --wait env --default-signal=INT time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -P $instances -I {} bash {} &
		#pgid=$!
		#wait $pgid

		# wait: capture e.g. sigint wich kills wait but setid job still running via kill at return trap
		# time: is also a bash shell keyword which works differently and does not record memory usage
		# due to set -E based error tracing, command time may leads to *** longjmp causes uninitialized stack frame ***: bash terminated
		# workaround: use full path, which, env or $(command -v time) <- prefer env to use env bash too

		# solution2:
		printf '%q\n' "${scripts[@]}" | env time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" parallel --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer -P "$tmpdir/instances.$jobname" -I {} bash {}
	else
		# printf '%q\n' "${scripts[@]}" | setsid --wait env --default-signal=INT xargs -P $instances -I {} bash {} &
		# pgid=$!
		# wait $pgid
		printf '%q\n' "${scripts[@]}" | parallel --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer -P "$tmpdir/instances.$jobname" -I {} bash {}
	fi

	return 0
}

commander::runalter_xargs(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of parallel
			-p <id|name>   | xargs process id or job name
		EOF
		return 1
	}

	local OPTIND arg mandatory pid instances
	while getopts 'p:i:' arg; do
		case $arg in
			p)	((++mandatory)); pid=$OPTARG;;
			i)	((++mandatory)); instances=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	BASHBONE_ERROR="no such job name: $pid"
	[[ $pid =~ ^[0-9]+$ ]] || pid=$(($(ps -o pgid= -p $(pgrep -o -f "job.$pid.*.sh"))))
	BASHBONE_ERROR="not a valid job id: $pid"
	[[ "$(ps -o comm= -p $pid)" == xargs ]]

	local i=$(pgrep -c -P $pid)
	i=$((i-instances))
	[[ $i -eq 0 ]] && return 0

	if [[ $i -gt 0 ]]; then
		# decrement
		while [[ $((i--)) -gt 0 ]]; do
			kill -USR2 $pid
			sleep 0.1
		done
	else
		# increment
		while [[ $((i++)) -lt 0 ]]; do
			kill -USR1 $pid
			sleep 0.1
		done
	fi

	return 0
}

commander::runalter(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of parallel
			-p <id|name>   | gnu parallel process id or job name
		EOF
		return 1
	}

	local OPTIND arg mandatory pid instances
	while getopts 'p:i:' arg; do
		case $arg in
			p)	((++mandatory)); pid=$OPTARG;;
			i)	((++mandatory)); instances=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	BASHBONE_ERROR="no such job name: $pid"

	[[ $pid =~ ^[0-9]+$ ]] || pid=$(($(ps -o ppid= -p $(pgrep -o -f "job.$pid.*.sh"))))
	BASHBONE_ERROR="not a valid job id: $pid"
	[[ "$(ps -o cmd= -p $pid)" =~ perl[[:space:]]+[^[:space:]]+parallel[[:space:]] ]]
	echo $instances > "$(ps -o cmd= --ppid $pid | head -1 | sed -E 's/^bash\s+//' | xargs -I {} bash -c 'echo "$(dirname {})/instances.$(basename {} | cut -d "." -f 2)"')"

	return 0
}

commander::runstat(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-p <id|name>   | xargs process id or job name
			-u <name>      | user
		EOF
		return 1
	}

	local OPTIND arg pid name user=".+"
	while getopts 'p:u:' arg; do
		case $arg in
			p)	pid="$OPTARG";;
			u)	user="$OPTARG";;
			*)	_usage;;
		esac
	done

	declare -a mapdata pgids
	mapfile -t mapdata < <(ps -o lstart,pid,ppid,pgid,rss,etime,user,comm,args -U $(users | tr -s ' ' ',') | sed -E '1!{s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/\1-\2-\3-\4-\5/}')

	if [[ $pid =~ ^[0-9]+$ ]]; then
		name=".+"
	elif [[ $pid ]]; then
		name="$pid"
		pid="\d+"
	else
		name=".+"
		pid="\d+"
	fi

	{	echo "${mapdata[0]}" | perl -lane '$F[4]="MEMORY"; print join" ",("NAME",$F[0],"JOBID",@F[4..6],"STATE","TASKID")'
		printf "%s\n" "${mapdata[@]}" | grep -E -f <(printf "%s\n" "${mapdata[@]}" | grep -E 'perl\s\S+/parallel\s' | awk '{print "^\\S+\\s+[0-9]+\\s+"$2"\\s+.+job\\..+\\.[0-9]\\.sh$"}') | perl -M'List::Util qw(max)' -slane 'if(1==1){$dirs{$F[2]}=$_=~s/.+bash (\/.*)job\.$n\.\d+\.sh$/$1/r; $jobs{$F[2]}=0; $names{$F[2]}=(split/\./,$F[-1])[-3]; $users{$F[2]}=$F[6]; $F[-1]=~s/.+\.(\d+)\.sh$/$1/; $jobs{$F[2]}=max($jobs{$F[2]},$F[-1]); print join" ",($names{$F[2]},$F[0],$F[1],@F[4..6],"r",$F[-1])}; END{for my $j (keys %jobs){ open F, "<$dirs{$j}/info.$names{$j}"; chomp(@l=<F>); close F; $njobs=$l[1]; $time=$l[2]; @l=($jobs{$j}); open F, "<$dirs{$j}/exitcodes.$names{$j}"; while(<F>){$_=~/^\S+\.(\d+)\s/; push @l,$1}; close F; $nrun=max(@l); print join" ",($names{$j},$time,$j,"* *",$users{$j},"q",++$nrun."-".$njobs) if $nrun < $njobs; } }' -- -n="$name" -p="$pid" -u="$user" | sort -k4,4V -k11,11n
	} | column -t || true

	return 0
}

commander::qsubcmd(){
	local jobid pid
	_cleanup::commander::qsubcmd(){
		[[ $jobid ]] && qdel $jobid &> /dev/null || true
		[[ $pid ]] && { env kill -PIPE $pid && wait $pid; } &> /dev/null || true
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v             | verbose on
			-w             | do wait for jobs and receive a non-null exit code if a single job fails
			-b             | benchmark on if do wait
			-n <name>      | prefix of logs and jobs to wait for - should be unique
			-c <env>       | run with conda
			-l <complex>   | sge digestable list of consumables as key="value" pairs (see qconf -sc or qconf -mc)
			-o <path>      | shared among machines for scripts, logs and exit codes
			-r             | override existing logs
			-i <instances> | number of parallel instances (in- or decrease afterwards via qalter -tc <instances> <name>)
			-q <queue>     | name of sge queue
			-p <env>       | name of parallel sge environment
			-t <threads>   | to be allocated per instance in parallel environment
			-s <idx[:idx]> | submit only jobs from cmds array starting from given index or range (default: 1)
			-d <jobid|name>| depends on and start after
			-a <cmds>      | array of
			example:
			${FUNCNAME[1]} -v -l hostname="!bcl102&!bcl103" -l mem_free="50G" -c base -p threads -t 4 -i 2 -w -o ~/logs -a cmd
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 instances verbose=false benchmark=false dowait="n" override=false cenv penv queue logdir complex params startid=1 stopid depends
	declare -n _cmds_qsubcmd # be very careful with circular name reference
	declare -a mapdata complexes logs
	while getopts 'vbwrt:i:o:l:p:q:c:n:a:s:d:' arg; do
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
			a)	((++mandatory)); _cmds_qsubcmd=$OPTARG;;
			s)	IFS=":" read -r startid stopid <<< "$OPTARG";;
			d)	depends="-hold_jid $OPTARG";; # do not quote optarg!
			*)	_usage;;
		esac
	done

	[[ $mandatory -lt 3 ]] && _usage
	[[ $_cmds_qsubcmd ]] || return 0
	declare -f "$BASHBONE_HOOKCMD" &> /dev/null && $BASHBONE_HOOKCMD _cmds_runcmd

	$verbose && commander::printinfo "running commands of array ${!_cmds_qsubcmd}"

	[[ $stopid ]] || stopid=${#_cmds_qsubcmd[@]}
	[[ $instances ]] || instances=$((stopid-startid+1))

	[[ $penv ]] && params="$penv $threads" || params="$queue"
	[[ $depends ]] && params+=" $depends"

	[[ $jobname ]] || jobname="$(basename "$(mktemp -u -p "$logdir" XXXXXXXXXX)")"
	local ex="$logdir/exitcodes.$jobname"
	local log="$logdir/job.$jobname.\$TASK_ID.log" # not SGE_TASK_ID

	export BASHBONE_SGEPATH="$PATH"
	local i id sh
	for i in "${!_cmds_qsubcmd[@]}"; do
		id=$((i+1))
		if [[ $id -ge $startid && $id -le $stopid ]]; then
			$override && rm -f "$logdir/job.$jobname.$id.log" "$ex"
			touch "$logdir/job.$jobname.$id.log" # for tail -f
			logs+=("$logdir/job.$jobname.$id.log")
		fi
		sh="$logdir/job.$jobname.$id.sh"
		echo '#!/usr/bin/env bash' > "$sh"
		echo "exit::$jobname.$id(){" >> "$sh"
		echo "    echo \"$jobname.$id (\$JOB_ID) exited with exit code \$1\" >> '$ex'" >> "$sh"
		[[ "$dowait" == "y" ]] && echo '    [[ $1 -gt 0 ]] && qdel $JOB_ID &> /dev/null || true' >> "$sh"
		echo '}' >> "$sh"
		echo 'PATH="${BASHBONE_SGEPATH:-$PATH}"' >> "$sh"
		echo "BASHBONE_SETSID=false" >> "$sh"

		if [[ $cenv ]]; then
			echo "source '$BASHBONE_DIR/activate.sh' -c true -x exit::$jobname.$id -i '$BASHBONE_TOOLSDIR'" >> "$sh"
			echo "conda activate --no-stack $cenv" >> "$sh"
		else
			echo "source '$BASHBONE_DIR/activate.sh' -c false -x exit::$jobname.$id -i '$BASHBONE_TOOLSDIR'" >> "$sh"
		fi
		$verbose && echo 'tail -2 "$0" | head -1 | paste -d " " <(echo ":CMD:") -' >> "$sh"
		printf '%s\n' "${_cmds_qsubcmd[$i]}" >> "$sh"
		echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
		chmod 755 "$sh"
	done

	# compared to /usr/bin/time, bash builtin time can handle: time echo "sleep 2" | bash
	# cons:
	# - benchmarks only qsub not job itself
	# - no memory consumption logged
	# - in case of job exit code > 0, leads to *** longjmp causes uninitialized stack frame ***: bash terminated
	# TIMEFORMAT=':BENCHMARK: runtime %3lR [hours][minutes]seconds'
	# time echo "$logdir/job.$jobname.\$SGE_TASK_ID.sh" | qsub -sync $dowait $params ${complexes[@]} -t 1-$id -tc $instances -S "$(/usr/bin/env bash -c 'which bash')" -V -cwd -o "$log" -j y -N $jobname |& sed -u -E '/exited/!d; s/Job [0-9]+\.(.+)\./job.'$jobname'.\1/;t;s/Job [0-9]+ (.+)\./job.'$jobname'.1 \1/'
	# attention: -S /bin/bash cannot be -S "/bin/bash --noprofile" and thus sources bash_profile and bashrc which in worst case modifies PATH so that conda may not serve its binaries first

	if [[ "$dowait" == "y" ]]; then
		local l
		tail -q -f "${logs[@]}" & pid=$!
		while read -r l; do
			[[ $jobid ]] || jobid=$(cut -d '.' -f 1 <<< $l)
			# requires 1>&2 : echo "$l" | sed -E '/exited/!d; s/Job ([0-9]+)\.(.+)\./\1 job.'$jobname'.\2/;t;s/Job ([0-9]+) (.+)\./\1 job.'$jobname'.1 \2/'
		done < <(echo "$logdir/job.$jobname.\$SGE_TASK_ID.sh" | BASH_EXECUTION_STRING="shournal" qsub -terse -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$(env bash -c 'which bash')" -V -cwd -o "$log" -j y -N $jobname 2> /dev/null || true)

		# use command/env qstat in case someone like me makes use of an alias :)
		# wait until accounting record is written to epilog after jobs post-processing metrics collection
		while env qstat -j $jobid &> /dev/null; do
			sleep 1
		done
		touch "$ex" #nfs requirend to make file visible in current shell on current node
		ex=$(awk '{print $NF}' "$ex" | sort -rg | head -1)
		# ex=$(qacct -j $jobid | awk '/^exit_status/{if($NF>x){x=$NF}}END{print x}')
		if $benchmark; then
			qacct -j $jobid | perl -M'List::Util qw(min max)' -lanE 'if($F[0] eq "start_time"){$d=join(" ",@F[1..$#F]); $d=`date -d "$d" +%s`; $sta=min($d,$sta?$sta:$d)} if($F[0] eq "end_time"){$d=join(" ",@F[1..$#F]); $d=`date -d "$d" +%s`; $sto=max($d,$sto?$sto:$d)} END{$s=$sto-$sta; if($s>3600){$d=3600}else{$d=60}; $hm=sprintf("%.0d",$s/$d); $hm=0 unless $hm; $ms=sprintf("%05.2f",($s/$d-$hm)*60); say ":BENCHMARK: runtime $hm:$ms [hours:]minutes:seconds"}'
			printf ":BENCHMARK: memory %s Kbytes\n" $(qacct -j $jobid | sed -nE 's/^ru_maxrss\s+([0-9.]+)\s*$/\1/p' | sort -gr | head -1)
		fi
		return $ex
	else
		echo "$logdir/job.$jobname.\$SGE_TASK_ID.sh" | BASH_EXECUTION_STRING="shournal" qsub -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$(env bash -c 'which bash')" -V -cwd -o "$log" -j y -N $jobname
		return 0
	fi
}

commander::qalter(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of parallel
			-l <complex>   | sge digestable list of consumables as key="value" pairs (see qconf -sc or qconf -mc)
			-p <id|name>   | xargs process id or job name
		EOF
		return 1
	}

	local OPTIND arg mandatory pid instances
	declare -a complexes
	while getopts 'p:i:l:' arg; do
		case $arg in
			p)	((++mandatory)); pid=$OPTARG;;
			i)	((++mandatory)); instances="-tc $OPTARG";;
			l)	complexes+=("-l $OPTARG");;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	qalter ${complexes[@]} $instances $pid

	return 0
}

commander::qstat(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-p <id|name>   | job id or job name
			-u <name>      | user
		EOF
		return 1
	}

	local OPTIND arg pid name user=".+"
	while getopts 'p:u:' arg; do
		case $arg in
			p)	pid="$OPTARG";;
			u)	user="$OPTARG";;
			*)	_usage;;
		esac
	done

	if [[ $pid =~ ^[0-9]+$ ]]; then
		name=".+"
	elif [[ $pid ]]; then
		name="$pid"
		pid="\d+"
	else
		name=".+"
		pid="\d+"
	fi

	{	echo "JOBID PRIOR NAME USER STATE STARTED QUEUE SLOTS TASKID"
		env qstat -xml | tr '\n' ' ' | sed 's/<job_list[^>]*>/\n/g;s/<[^>]*>//g' | tr -s ' ' ' ' | perl -slane 'next unless $F[0]=~/^$p$/ || $F[2]=~/^$n$/ || $F[3]=~/^$u$/; unless($F[8]){$F[8]=$F[7]; $F[7]=$F[6]; $F[6]="*"} print join" ",@F; ' -- -p="$pid" -n="$name" -u="$user"
	} | column -t

	return 0
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

	commander::runcmd -v -b -i $threads -a cmd || commander::printerr "failed"

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
