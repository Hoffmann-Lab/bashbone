#! /usr/bin/env bash
# (c) Konstantin Riege

genome::mkdict() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, indexing respectively
			-t <threads>  | number of
			-i <genome>   | path to
			-p <tmpdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads genome tmpdir skip=false skipmd5=false
	while getopts 'S:s:5:t:i:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); genome="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	commander::printinfo "creating genome dictionary"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and genome dictionary creation respectively"
	} || {
		commander::printinfo "checking md5 sums"

		local instances ithreads jmem jgct jcgct
		read -r instances ithreads jmem jgct jcgct < <(configure::jvm -T $threads)
		declare -a cmd1 cmd2

		local dict="$(mktemp -u -p "$tmpdir" cleanup.XXXXXXXXXX.dict)"
		commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
			picard
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				CreateSequenceDictionary
				R="$genome"
				O="$dict"
				VERBOSITY=WARNING
		CMD

		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			grep -Eo 'SN:\S+' "$dict" | cut -d ':' -f 2- > "$genome.list"
		CMD
			mv "$dict" "${genome%.*}.dict"
		CMD

		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
			samtools faidx "$genome"
		CMD

		$skip && {
			commander::printcmd -a cmd1
			commander::printcmd -a cmd2
		} || {
			{	commander::runcmd -c picard -v -b -t $threads -a cmd1 || return 1
				local md5dict thismd5genome thismd5dict
				md5dict=$(md5sum "$dict" | cut -d ' ' -f 1)
				thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
				[[ -s "${genome%.*}.dict" ]] && thismd5dict=$(md5sum "${genome%.*}.dict" | cut -d ' ' -f 1)
				if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5dict" || "$thismd5dict" != "$md5dict" ]]; then
					commander::runcmd -v -b -t $threads -a cmd2 || return 1
				fi
			} || {
				rm -f "$dict"
				commander::printerr "$funcname failed"
				return 1
			}
		}
	}

	rm -f "$dict"
	return 0
}
