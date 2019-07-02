#!/usr/bin/env bash
# (c) Konstantin Riege

declare -a COMMANDER

commander::print(){
	[[ $* ]] && echo ":INFO: $*"
	local fd
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t
		printf ':INFO: %s\n' "${MAPFILE[@]}"
		# while read -u $fd -r tmp; do
		# 	echo ":INFO: $tmp"
		# done
	done
	COMMANDER=()
	return 0
}

commander::warn(){
	[[ $* ]] && echo ":WARNING: $*"
	local fd
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t
		printf ':WARNING: %s\n' "${MAPFILE[@]}"
	done
	COMMANDER=()

	return 0
}

commander::printerr(){
	echo -ne "\e[0;31m"
	[[ $* ]] && echo ":ERROR: $*" >&2
	local fd
	for fd in "${COMMANDER[@]}"; do
		mapfile -u $fd -t
		printf ':ERROR: %s\n' "${MAPFILE[@]}" >&2
	done
	COMMANDER=()
	echo -ne "\e[m"

	return 0
}

commander::makecmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
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
	declare -a cmd_makecmd # be very careful with references name space
	while getopts ':a:o:s:c' arg; do
		case $arg in
			a)	mendatory=1; _cmds_makecmd=$OPTARG;;
			s)	sep=$(echo -e "${OPTARG:- }");; # echo -e to make e.g. '\t' possible
			o)	suffix=' > '"$OPTARG";;
			c)	[[ ! $mendatory ]] && { _usage; return 1; }
				shift $((OPTIND-1)) # remove '-a <cmd> -s <char> -o <file> -c' from $*
				for fd in "${COMMANDER[@]}"; do
					mapfile -u $fd -t
					cmd_makecmd+=("${MAPFILE[*]}") # * instead of @ to concatenate
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
		commander::printerr {COMMANDER[0]}<<- EOF
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
				printf ':CMD: %s\n' "${_cmds_printcmd[@]}"
				return 0;;
			*)	_usage; return 1;;
		esac
	done

	_usage; return 1
}

commander::runcmd(){
	trap 'rm -rf /dev/shm/commander::runcmd::$$' RETURN INT TERM

	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v           | verbose on
			-b           | benchmark on
			-t <threads> | number of
			-a <cmds>    | ALWAYS LAST OPTION
			               array of
			example:
			$funcname -v -b -t 1 -a cmd
		EOF
		return 0
	}

	local OPTIND arg threads=1 verbose=false benchmark
	declare -n _cmds_runcmd # be very careful with circular name reference
	while getopts 'vbt:a:' arg; do
		case $arg in
			t)	threads=$OPTARG;;
			v)	verbose=true;;
			b)	benchmark=true;;
			a)	_cmds_runcmd=$OPTARG
				[[ $_cmds_runcmd ]] || return 0
				$verbose && {
					echo ":INFO: running commands of array ${!_cmds_runcmd}"
					commander::printcmd -a _cmds_runcmd
				}
				# better write to file to avoid xargs argument too long error
				local shdir="/dev/shm/commander::runcmd::$$"
				mkdir -p $shdir
				local i md5sh
				if [[ $benchmark ]]; then
					# printf '%s\0' "${_cmds_runcmd[@]}" | command time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -0 -P $threads -I {} bash -c {}
					for i in "${!_cmds_runcmd[@]}"; do
						md5sh=$(printf '%s\n' "$i${_cmds_runcmd[$i]}" | md5sum | cut -d ' ' -f 1)
						printf '%s\n' "${_cmds_runcmd[$i]}" > "$shdir/$md5sh"
						echo "$shdir/$md5sh"
					done | command time -f ":BENCHMARK: runtime %E [hours:]minutes:seconds\n:BENCHMARK: memory %M Kbytes" xargs -P $threads -I {} bash {}
				else
					# printf '%s\0' "${_cmds_runcmd[@]}" | xargs -0 -P $threads -I {} bash -c {}
					for i in "${!_cmds_runcmd[@]}"; do
						md5sh=$(printf '%s\n' "$i${_cmds_runcmd[$i]}" | md5sum | cut -d ' ' -f 1)
						printf '%s\n' "${_cmds_runcmd[$i]}" > "$shdir/$md5sh"
						echo "$shdir/$md5sh"
					done | xargs -P $threads -I {} bash {}
				fi
				return $((${PIPESTATUS[@]/%/+}0));;
			*)	_usage;	return 1;;
		esac
	done

	_usage; return 1
}

commander::_test(){
	local x=${1:-'hello world'} threads=${2:-1} i=3 s

	while read -u $((i++)) -r s 2> /dev/null; do
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

	commander::makecmd -a cmd -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
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

	commander::runcmd -v -b -t $threads -a cmd || commander::printerr "failed"

	commander::printerr ${FUNCNAME[0]} EXPECTED OUTPUT
	commander::printerr {COMMANDER[0]}<<-OUT
		foo
		bar
		baz
		:INFO: running commands of array cmd
		:CMD: perl -sl - -x='hello world' <<< 'print "$x"' | awk '{print $0}' 
		:CMD: perl -sl - -x='hello world' <<< ' print "$x"; ' | awk '{print $0}'
		:CMD: perl -l - <<< ' print "hello world"; ' | awk '{print $0}'
		:CMD: perl -l - <<< ' print "hello world"; ' |awk '{print $0}'
		:CMD: perl -s - -x='hello world' <<< ' print "$x";  print " hello world\n"; '
		hello world
		hello world
		hello world
		hello world
		hello world hello world
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
