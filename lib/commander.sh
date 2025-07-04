#!/usr/bin/env bash
# (c) Konstantin Riege

declare -a COMMANDER

function commander::print(){
	if [[ $COMMANDER ]]; then
		local fd
		declare -a mapdata
		for fd in "${COMMANDER[@]}"; do
			mapfile -u $fd -t mapdata
			printf '%s\n' "${mapdata[@]}"
		done
		COMMANDER=()
	else
		if [[ $1 ]]; then
			echo "$*"
		else
			cat
		fi
	fi

	return 0
}

function commander::printinfo(){
	if [[ $COMMANDER ]]; then
		local fd
		declare -a mapdata
		for fd in "${COMMANDER[@]}"; do
			mapfile -u $fd -t mapdata
			printf ':INFO: %s\n' "${mapdata[@]}"
		done
		COMMANDER=()
	else
		if [[ $1 ]]; then
			echo ":INFO: $*"
		else
			awk '{print ":INFO: "$0}'
		fi
	fi

	return 0
}

function commander::warn(){
	if [[ $COMMANDER ]]; then
		local fd
		declare -a mapdata
		for fd in "${COMMANDER[@]}"; do
			mapfile -u $fd -t mapdata
			printf ':WARNING: %s\n' "${mapdata[@]}"
		done
		COMMANDER=()
	else
		if [[ $1 ]]; then
			echo ":WARNING: $*"
		else
			awk '{print ":WARNING: "$0}'
		fi
	fi

	return 0
}

function commander::printerr(){
	if [[ $COMMANDER ]]; then
		local fd
		declare -a mapdata
		for fd in "${COMMANDER[@]}"; do
			mapfile -u $fd -t mapdata
			printf ':ERROR: %s\n' "${mapdata[@]}" >&2
		done
		COMMANDER=()
	else
		if [[ $1 ]]; then
			echo ":ERROR: $*" >&2
		else
			awk '{print ":ERROR: "$0}' >&2
		fi
	fi

	return 0
}

function commander::make(){
	commander::makecmd "$@"
}

function commander::makecmd(){
	function _usage(){
		commander::print <<- 'EOF'
			commander::makecmd

			easily compose and store command in an array for later execution

			-a <cmds>          | mandatory. array to append the command onto
			-v <variable>      | optional. env variable(s) necessary for command. can be used multiple times
			-s <separator>     | optional. character to separate command chunks, when file descriptor array COMMANDER is used. default: '\n'
			-m                 | force multi line support in case of activated bashbone legacy mode, which inerts line breaks
			-o <outfile>       | optional. redirect stdout to output file
			-O <outfile>       | optional. append stdout to output file
			-c <cmd|{fd[0]}..> | ALWAYS LAST OPTION. command line string(s) and or here-documents or file descriptor array COMMANDER

			example 1:
			commander::makecmd -a cmds -c perl -le \''print "foo"'\'

			example 2:
			x=1
			commander::makecmd -a cmds -v x -c echo '$x'

			example 3:
			x=1
			commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			    perl -sl - -y=$x <<< '
			        print "$x";
			        print "\$y";
			    '
			CMD
			    awk '{print $0}'
			CMD

			example 4:
			x=1
			commander::makecmd -a cmds -v x -c <<-'CMD'
			    perl -sl - -y=$x <<< '
			        print "$x";
			        print "$y";
			    ' | awk '{print $0}'
			CMD
		EOF
		return 1
	}

	local OPTIND arg mandatory sep suffix vars='' multiline
	# declare -a vars
	declare -n _cmds_makecmd # be very careful with circular name reference
	while getopts 'v:a:o:O:s:mc' arg; do
		case $arg in
			v)	vars+="$(declare -p $OPTARG)"$'\n';;
			a)	((++mandatory)); _cmds_makecmd=$OPTARG;;
			s)	sep=$(echo -e "$OPTARG");;
			m)	multiline=true;;
			o)	suffix=" > '$OPTARG'";;
			O)	suffix=" >> '$OPTARG'";;
			c)	((++mandatory)); shift $((OPTIND-1)); break;; # do not use getopts c: which requires shift OPTIND-2 and to implement :) case with shift OPTIND-1 if OPTARG=="c" and return 1 if OPTARG!="c"
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage
	[[ ! $multiline ]] && ${BASHBONE_LEGACY:-false} && multiline=false
	[[ $sep ]] || { $multiline && sep=$'\n' || sep=''; }

	local fd tmp
	declare -a mapdata cmd_makecmd # be very careful with references name space

	if [[ $COMMANDER ]]; then
		for fd in "${COMMANDER[@]}"; do
			if $multiline; then
				cmd_makecmd+=("$(< /dev/fd/$fd)")
			else
				mapfile -u $fd -t mapdata
				cmd_makecmd+=("${mapdata[*]}") # * instead of @ to concatenate a non-splitted sentence instead of single non-splitted words
			fi
			exec {fd}>&- # just to be safe. closing fd not necessary since heredoc is read only and handles closure after last EOF
		done
		tmp="${cmd_makecmd[*]/#/$sep }" # concatenate CMD* with separator. do not use $(echo -e ${tmp/#$sep /}) here
		tmp="${tmp/#$sep /}"
		_cmds_makecmd+=("$vars$tmp$suffix")
	else
		if [[ $1 ]]; then
			_cmds_makecmd+=("$vars$*$suffix")
		else
			if $multiline; then
				tmp="$(< /dev/stdin)"
			else
				mapfile -t mapdata /dev/stdin
				tmp="${mapdata[*]}"
			fi
			_cmds_makecmd+=("$vars$tmp$suffix")
		fi
	fi
	COMMANDER=()

	return 0
}

function commander::printcmd(){
	function _usage(){
		commander::print <<- EOF
			commander::printcmd

			prints the commands stored in an array

			-a <cmds> | mandatory. array of commands to be printed
		EOF
		return 1
	}

	local OPTIND arg mandatory
	declare -n _cmds_printcmd # be very careful with circular name reference
	while getopts 'a:' arg; do
		case $arg in
			a)	((++mandatory)); _cmds_printcmd=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 1 ]] && _usage

	local i
	for i in $(seq 0 $((${#_cmds_printcmd[@]}-1))); do
		echo ":CMD$((i+1)): start"
		printf '%s\n' "${_cmds_printcmd[$i]}"
		echo ":CMD$((i+1)): end"
	done

	return 0
}

function commander::run(){
	commander::runcmd "$@"
}

function commander::runcmd(){
	function _usage(){
		commander::print <<- EOF
			commander::runcmd

			orchestrates the execution of commands stored within an array by utilizing GNU parallel

			-v             | optional. toggle on verbose mode which will print the commands prior to their execution
			-b             | optional. toggle on benchmarking run time and peak memory consumption
			-c <env>       | optional. activate a bashbone conda environment prior to command execution
			-i <instances> | optional. number of parallel instances. see also commander::runalter. default: 1
			-t <instances> | optional. obsolete synonym for -i
			-s <idx[:idx]> | optional. restrict the execution of commands by array index start or range. default: 1
			-n <name>      | optional. prefix for logs and job scripts
			-o <path>      | optional. output directory for logs, job scripts and exit codes
			-r             | optional. override existing logs. default: append
			-a <cmds>      | mandatory. array of commands. see also commander::makecmd

			example 1:
			cmds=("echo job1" "echo job2" "echo job3")
			cmds+=("sleep 2; echo job4")
			commander::runcmd -v -b -i 2 -a cmds

			example 2:
			commander::runcmd -i 2 -n current -o ~/jobs -r -a cmds
		EOF
		return 1
	}

	local OPTIND arg mandatory instances=1 verbose=false benchmark=false cenv startid=1 stopid override=false jobname="current" logdir
	declare -n _cmds_runcmd # be very careful with circular name reference
	while getopts 'vbt:i:c:a:s:o:n:r' arg; do
		case $arg in
			t)	instances=$OPTARG;; # obsolete, for compatibility
			i)	instances=$OPTARG;;
			v)	verbose=true;;
			b)	benchmark=true;;
			c)	cenv=$OPTARG;;
			a)	((++mandatory)); _cmds_runcmd=$OPTARG;;
			s)	IFS=":" read -r startid stopid <<< "$OPTARG";;
			r)	override=true;;
			o)	logdir="$OPTARG"; mkdir -p "$logdir"; logdir="$(realpath -s "$logdir")";;
			n)	jobname="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 1 ]] && _usage
	[[ $_cmds_runcmd ]] || return 0
	"${BASHBONE_HOOKCMD:-:}" _cmds_runcmd

	$verbose && commander::printinfo "running commands of array ${!_cmds_runcmd}"

	local i id sh ex log tmpdir
	[[ $logdir ]] && tmpdir="$logdir" || tmpdir="$(mktemp -d -p "${TMPDIR:-/tmp}" jobs.XXXXXXXXXX)"
	[[ $jobname ]] || jobname="$(command mktemp -u XXXXXXXXXX)"
	echo $instances > "$tmpdir/instances.$jobname"
	ex="$tmpdir/exitcodes.$jobname"
	touch "$ex" # for runstat
	[[ $stopid ]] || stopid=${#_cmds_runcmd[@]}

	# attention: parallel sources bashrc which in worst case modifies PATH so that bashbone may not serve its binaries first
	# solution: use a fake shell via PARALLEL_SHELL to prevent sourcing bashrc
	# .. and do not exec bash "$@" which will source bashrc
	cat <<- 'EOF' > "$tmpdir/shell.$jobname"
		#!/usr/bin/env bash
		bash "$@"
	EOF
	chmod 755 "$tmpdir/shell.$jobname"

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
		id=$((i+1))
		sh="$tmpdir/job.$jobname.$id.sh"
		log="$tmpdir/job.$jobname.$id.log"
		$override && rm -f "$log" "$ex"
		echo '#!/usr/bin/env bash' > "$sh"
		echo "_on_exit(){ :; }" >> "$sh"
		echo "exit::$jobname.$id(){" >> "$sh"
		echo "    echo \"$jobname.$id (\$((\$(ps -o ppid= -p \$\$ 2> /dev/null)))) exited on \$HOSTNAME with exit code \$1\" >> '$ex'" >> "$sh"
		echo "    _on_exit \$1" >> "$sh"
		echo '}' >> "$sh"
		if [[ $cenv ]]; then
			echo "source '$BASHBONE_DIR/activate.sh' -l ${BASHBONE_LEGACY:-true} -s '$BASHBONE_EXTENSIONDIR' -c true -r false -x exit::$jobname.$id -i $BASHBONE_TOOLSDIR" >> "$sh"
			echo "conda activate --no-stack $cenv" >> "$sh"
		else
			echo "source '$BASHBONE_DIR/activate.sh' -l ${BASHBONE_LEGACY:-true} -s '$BASHBONE_EXTENSIONDIR' -c false -r false -x exit::$jobname.$id -i $BASHBONE_TOOLSDIR" >> "$sh"
		fi
		$verbose && echo "echo :BASE64$id: $(printf '%s\n' "${_cmds_runcmd[$i]}" | base64 -w 0)" >> "$sh"
		# redirect sterr to original stdout FD first to not duplicate stream without >&2
		echo "exec 2> >(trap '' INT TERM; exec stdbuf -o L tee -a '$log' | sed -u 's/^/:STDERR:/')" >> "$sh"
		echo "exec 1> >(trap '' INT TERM; exec stdbuf -o L tee -a '$log')" >> "$sh"

		printf '%s\n' "${_cmds_runcmd[$i]}" >> "$sh"
		echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
		scripts+=("$(realpath -se "$sh")") # necessary for runstat
	done

	if $benchmark; then
		# solution1: setsid xargs, and send termination sequence to its PGID upon receiving INT via ctr+c e.g. by using return trap
		# needs also termination sequence on PGID in exit function to kill all sibling processes upon error
		# printf '%q\n' "${scripts[@]}" | setsid --wait env --default-signal=INT time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -P $instances -I {} bash {} &
		# pgid=$!
		# wait $pgid

		# solution2: gnu parallel
		# drawback: in contrast to xargs, where USR signal to increase instace count has an instant effect, parallel reads instances file only after completing one job
		# -> implement somehow respawning dummy job?
		printf '%q\n' "${scripts[@]}" | env PARALLEL_SHELL="$tmpdir/shell.$jobname" time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" parallel --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer -P "$tmpdir/instances.$jobname" -I {} bash {}
		# time: is also a bash shell keyword which works differently and does not record memory usage
		# due to set -E based error tracing, command time may leads to *** longjmp causes uninitialized stack frame ***: bash terminated
		# workaround: use full path, which, env or command -v <- prefer env to use env bash too
	else
		printf '%q\n' "${scripts[@]}" | PARALLEL_SHELL="$tmpdir/shell.$jobname" parallel --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer -P "$tmpdir/instances.$jobname" -I {} bash {}
	fi | stdbuf -o L awk  '/^:BASE64[0-9]*:/{sub(":BASE64",":CMD",$1); print $1" start"; while(("echo "$2" | base64 -d" | getline l)){print l} print $1" end"; next} /^:STDERR:/{sub(/^:STDERR:/, "", $0); print > "/dev/stderr"; next} {print}'

	return 0
}

function commander::runalter(){
	function _usage(){
		commander::print <<- EOF
			commander::runalter

			increment or decrement the number of instances i.e. parallel running jobs of a xargs or GNU parallel process

			-i <instances> | mandatory. number of new parallel instances
			-p <id|name>   | mandatory. GNU parallel process id or job name. see also commander::runstat
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 2 ]] && _usage

	BASHBONE_ERROR="no such job name: $pid"
	[[ $pid =~ ^[0-9]+$ ]] || pid=$(($(ps -o ppid= -p $(pgrep -o -f "job.$pid.*.sh"))))

	BASHBONE_ERROR="not a valid job id: $pid"
	if [[ "$(ps -o comm= -p $pid)" == xargs ]]; then
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
	elif [[ "$(ps -o cmd= -p $pid)" =~ perl[[:space:]]+[^[:space:]]+parallel[[:space:]] ]]; then
		echo $instances > "$(ps -o cmd= --ppid $pid | head -1 | sed -E 's/^bash\s+//' | xargs -I {} bash -c 'echo "$(dirname "$1")/instances.$(basename "$1" | cut -d "." -f 2)"' bash {})"
	else
		false
	fi

	return 0
}

function commander::runstat(){
	function _usage(){
		commander::print <<- EOF
			commander::runstat

			inspect currently running GNU parallel processes

			-p <id|name>   | optional. GNU parallel process id or job name
			-u <name>      | optional. user name
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

function commander::qsub(){
	commander::qsubcmd "$@"
}

function commander::qsubcmd(){
	function _usage(){
		commander::print <<- EOF
			commander::qsubcmd

			orchestrates the execution of commands stored within an array by submitting them as tasks of an array job into the specified queue, parallel environment repectively, of a Univa/Sun Grid Engine/Son of Grid Engine (SGE)

			-v             | optional. toggle on verbose mode which will print the commands prior to their execution
			-b             | optional. given -w option, toggles on benchmarking run time and peak memory consumption
			-w             | optional. wait for the completion of all job tasks and receive a non-null exit code if a single one fails
			-n <name>      | optional. prefix for logs and job scripts
			-c <env>       | optional. activate a bashbone conda environment prior to command execution
			-l <complex>   | optional. list of consumables as key="value" pairs. can be used multiple times. (see 'qconf -sc' or 'qconf -mc')
			-o <path>      | mandatory. output directory for logs, job scripts and exit codes. needs to be accessible by all machines
			-r             | optional. override existing logs. default: append
			-i <instances> | optional. number of parallel instances/tasks. see also commander::qalter. default: all
			-q <queue>     | mandatory unless -p option is used. name of the SGE queue, the job will be submitted to
			-p <env>       | mandatory unless -q option is used. name of the SGE parallel environment, the job will be submitted to
			-t <threads>   | optional. given the -p option, the number of threads equatable to number of sockets to be allocated per instance. default: 1
			-m <memory>    | optional. the hard limit of usable memory per job task in MB. applied if h_vmem consumable can be set and distributed across sockets according to -t option.
			-s <idx[:idx]> | optional. restrict the execution of commands by array index start or range. default: 1
			-d <jobid|name>| optional. the job which is going to be submitted depends on results of a currently running job and will be put on halt. see also commander::qstat
			-a <cmds>      | mandatory. array of commands. see also commander::makecmd

			example 1:
			cmds=("echo job1" "echo job2" "echo job3")
			cmds+=("sleep 2; echo job4")
			commander::qsubcmd -v -q all.q -i 2 -w -b -o ~/logs -a cmds

			example 2:
			commander::qsubcmd -v -l hostname="!server1&!server2" -l mem_free="150G" -p all.pe -i 2 -r -n first -o ~/logs -a cmds1
			commander::qsubcmd -v -l hostname="server1|server2" -p all.pe -t 4 -m 200000 -r -n second -r -o ~/logs -a cmds2 -d first
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 instances verbose=false benchmark=false dowait="n" override=false cenv penv q queue logdir complex params startid=1 stopid depends memory
	declare -n _cmds_qsubcmd # be very careful with circular name reference
	declare -a mapdata complexes logs
	while getopts 'vbwrt:m:i:o:l:p:q:c:n:a:s:d:' arg; do
		case $arg in
			v)	verbose=true;;
			b)	benchmark=true;;
			w)	dowait="y";;
			r)	override=true;;
			c)	cenv=$OPTARG;;
			t)	threads=$OPTARG;;
			m)	memory=$OPTARG;;
			i)	instances=$OPTARG;;
			o)	((++mandatory)); logdir="$OPTARG"; mkdir -p "$logdir"; logdir="$(realpath -s "$logdir")";;
			l)	complexes+=("-l $OPTARG");;
			p)	((++mandatory)); q=$OPTARG; penv="-pe $q";;
			q)	((++mandatory)); q=$OPTARG; queue="-q $q";;
			n)	jobname="$OPTARG";;
			a)	((++mandatory)); _cmds_qsubcmd=$OPTARG;;
			s)	IFS=":" read -r startid stopid <<< "$OPTARG";;
			d)	depends="-hold_jid $OPTARG";; # do not quote optarg!
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage
	[[ $_cmds_qsubcmd ]] || return 0
	"${BASHBONE_HOOKCMD:-:}" _cmds_runcmd

	local conf="$(mktemp -p "${TMPDIR:-/tmp}" $q.XXXXXXXXXX.conf)"
	qconf -sq $q > "$conf" 2> /dev/null && {
		if [[ $(awk '/^terminate_method/{print $NF}' "$conf") != SIGINT ]]; then
			commander::warn "SGE termination method cannot be set to SIGINT. Cleanup upon error or job deletion not possible!"
			#sed -Ei 's/(terminate_method\s+).*/\1SIGINT/' "$conf"
			#qconf -Mq "$conf"
		fi
	} || commander::warn "SGE configuration cannot be read from current server $HOSTNAME. Is SGE termination method set to SIGINT?"

	$verbose && commander::printinfo "running commands of array ${!_cmds_qsubcmd}"

	[[ $stopid ]] || stopid=${#_cmds_qsubcmd[@]}
	[[ $instances ]] || instances=$((stopid-startid+1))

	[[ $penv ]] && params="$penv $threads" || params="$queue"
	[[ $depends ]] && params+=" $depends"

	if [[ $memory ]] && qconf -sc | grep -Fw h_vmem | awk 'tolower($6)~/^no/{exit 1}'; then
		[[ $penv ]] && params+=" -l h_vmem=$((memory/threads))M" || params+=" -l h_vmem=${memory}M"
	fi

	[[ $jobname ]] || jobname="$(command mktemp -u XXXXXXXXXX)"
	local ex="$logdir/exitcodes.$jobname"
	local log="$logdir/job.$jobname.\$TASK_ID.log" # not SGE_TASK_ID

	# attention: -S /bin/bash cannot be -S "/bin/bash --noprofile" and thus sources bash_profile and bashrc which in worst case modifies PATH so that bashbone may not serve its binaries first
	# solution: use a fake shell to prevent sourcing bashrc
	# or better source rc to allow user configs like shorunal and don't submit current environment
	# cat <<- 'EOF' > "$logdir/shell.$jobname"
	# 	#!/usr/bin/env bash
	# 	exec bash "$@"
	# EOF
	# chmod 755 "$logdir/shell.$jobname"

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
		echo "_on_exit(){ :; }" >> "$sh"
		echo "exit::$jobname.$id(){" >> "$sh"
		echo "    echo \"$jobname.$id (\$JOB_ID) exited on \$HOSTNAME with exit code \$1\" >> '$ex'" >> "$sh"
		echo "    _on_exit \$1" >> "$sh"
		[[ "$dowait" == "y" ]] && echo '    [[ $1 -gt 0 ]] && qdel $JOB_ID &> /dev/null || true' >> "$sh"
		echo '}' >> "$sh"
		# re-execution via setsid required!
		# 1 so that any SGE version that sends kill to PID or PGID lets bashbone kill its BASHBONE_PGID and thereby performs cleanup via exit trap before getting cut from terminal/pty
		# 2 shournal can be executed explicitly, so that no fork runs in backround that will be killed otherwise and thus leaves unlogged write events
		if [[ $cenv ]]; then
			echo "source '$BASHBONE_DIR/activate.sh' -l ${BASHBONE_LEGACY:-true} -s '$BASHBONE_EXTENSIONDIR' -c true -r true -x exit::$jobname.$id -i $BASHBONE_TOOLSDIR" >> "$sh"
			echo "conda activate --no-stack $cenv" >> "$sh"
		else
			echo "source '$BASHBONE_DIR/activate.sh' -l ${BASHBONE_LEGACY:-true} -s '$BASHBONE_EXTENSIONDIR' -c false -r true -x exit::$jobname.$id -i $BASHBONE_TOOLSDIR" >> "$sh"
		fi
		$verbose && echo "echo :BASE64$id: $(printf '%s\n' "${_cmds_qsubcmd[$i]}" | base64 -w 0)" >> "$sh"
		printf '%s\n' "${_cmds_qsubcmd[$i]}" >> "$sh"
		echo "exit 0" >> "$sh" # in case last command threw sigpipe, exit 0
		chmod 755 "$sh"
	done

	if [[ "$dowait" == "y" ]]; then
		local l jobid
		{ stdbuf -o L tail -q -f "${logs[@]}" | stdbuf -o L awk  '/^:BASE64[0-9]*:/{sub(":BASE64",":CMD",$1); print $1" start"; while(("echo "$2" | base64 -d" | getline l)){print l} print $1" end"; next} /^:STDERR:/{sub(/^:STDERR:/, "", $0); print > "/dev/stderr"; next} {print}'; } &
		echo "{ env kill -TERM $!; wait $!; } &> /dev/null" >> "$BASHBONE_CLEANUP"
		while read -r l; do
			if [[ ! $jobid ]]; then
				jobid=$(cut -d '.' -f 1 <<< $l)
				echo "qdel $jobid &> /dev/null" >> "$BASHBONE_CLEANUP"
			fi
		done < <(echo "\"$logdir/job.$jobname.\$SGE_TASK_ID.sh\"" | qsub -terse -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$(which bash)" -cwd -o "$log" -j y -N $jobname 2> /dev/null || true)
		# done < <(echo "\"$logdir/job.$jobname.\$SGE_TASK_ID.sh\"" | qsub -terse -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$logdir/shell.$jobname" -V -cwd -o "$log" -j y -N $jobname 2> /dev/null || true)

		# use command/env qstat in case someone like me makes use of an alias :)
		# wait until accounting record is written to epilog after jobs post-processing metrics collection
		# while command qstat -j $jobid &> /dev/null not enought
		qconf -sh &> /dev/null && {
			while ! qacct -j $jobid &> /dev/null; do
				sleep 1
			done
		} || {
			commander::warn "SGE accountings cannot be read from current server $HOSTNAME"
			benchmark=false
			while command qstat -j $jobid &> /dev/null; do
				sleep 1
			done
		}

		sleep 1
		touch "$ex" #nfs requirend to make file visible in current shell on current node
		ex=$(awk '{print $NF}' "$ex" | sort -rg | head -1)
		# ex=$(qacct -j $jobid | grep -F exit_status | awk '{print $NF}' "$ex" | sort -rg | head -1)

		if $benchmark; then
			qacct -j $jobid | perl -M'List::Util qw(min max)' -lanE 'if($F[0] eq "start_time"){$d=join(" ",@F[1..$#F]); $d=`date -d "$d" +%s`; $sta=min($d,$sta?$sta:$d)} if($F[0] eq "end_time"){$d=join(" ",@F[1..$#F]); $d=`date -d "$d" +%s`; $sto=max($d,$sto?$sto:$d)} END{$s=$sto-$sta; if($s>3600){$d=3600}else{$d=60}; $hm=sprintf("%.0d",$s/$d); $hm=0 unless $hm; $ms=sprintf("%05.2f",($s/$d-$hm)*60); say ":BENCHMARK: runtime $hm:$ms [hours:]minutes:seconds"}'
			printf ":BENCHMARK: memory %s\n" $(qacct -j $jobid | sed -nE 's/^ru_maxrss\s+([0-9.]+)\s*(.*)$/\1 \2/p' | sort -gr | head -1)
		fi
		return $ex
	else
		echo "\"$logdir/job.$jobname.\$SGE_TASK_ID.sh\"" | qsub -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$(which bash)" -cwd -o "$log" -j y -N $jobname
		# echo "\"$logdir/job.$jobname.\$SGE_TASK_ID.sh\"" | qsub -sync $dowait $params ${complexes[@]} -t $startid-$stopid -tc $instances -S "$logdir/shell.$jobname" -V -cwd -o "$log" -j y -N $jobname
		return 0
	fi
}

function commander::qalter(){
	function _usage(){
		commander::print <<- EOF
			commander::qalter

			increment or decrement the number of parallel running jobs within a task of a Univa/Sun Grid Engine/Son of Grid Engine (SGE) and/or change other consumables

			-i <instances> | mandatory. number of new parallel instances
			-l <complex>   | optional. list of consumables as key="value" pairs. can be used multiple times. (see 'qconf -sc' or 'qconf -mc')
			-p <id|name>   | mandatory. SGE task id or name. see also commander::qstat
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

function commander::qstat(){
	function _usage(){
		commander::print <<- EOF
			commander::qstat

			inspect currently running tasks and jobs of a Univa/Sun Grid Engine/Son of Grid Engine (SGE)

			-p <id|name>   | optional. job id or job name
			-u <name>      | optional. user name
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
		command qstat -u "*" -xml | tr '\n' ' ' | sed 's/<job_list[^>]*>/\n/g;s/<[^>]*>//g' | tr -s ' ' ' ' | perl -slane 'next unless $F[0]=~/^$p$/ && $F[2]=~/^$n$/ && $F[3]=~/^$u$/; unless($F[8]){$F[8]=$F[7]; $F[7]=$F[6]; $F[6]="*"} print join" ",@F; ' -- -p="$pid" -n="$name" -u="$user"
	} | column -t

	return 0
}

function commander::_test(){
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
