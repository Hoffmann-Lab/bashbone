#! /usr/bin/env bash
# (c) Konstantin Riege

visualize::venn() {
	declare -a tdirs
	_cleanup::visualize::venn(){
		rm -rf "${tdir[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-l <idfiles>  | array of bed or id files
			-n <names>    | array of
			-p <tmpdir>   | path to
			-o <outfile>  | path to prefix
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads
	declare -n _lists_venn _names_venn
	while getopts 'S:s:t:l:b:n:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			l)	_lists_venn=$OPTARG;;
			b)	_beds_venn=$OPTARG;;
			n)	_names_venn=$OPTARG;;
			o)	((++mandatory)); outfile="$OPTARG"; mkdir -p "$(dirname "$outfile")";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "plotting venn diagrams"

	declare -a cmd1
	local params
	[[ $_names_venn ]] && params="--names $(printf '%s,' "${_names_venn[@]}" | sed 's/,$//')"
	[[ $(head -1 ${_lists_venn[0]} | awk '{print NF}') -eq 1 ]] && params+=" --type list"

	tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.intervene)")
	commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		intervene venn
		$params
		-i $(printf '"%s" ' "${_lists_venn[@]}")
		-o "${tdirs[-1]}"
		--figtype pdf
	CMD
		mv "${tdirs[-1]}/Intervene_venn.pdf" "$outfile.pdf"
	CMD

	tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.intervene)")
	commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		intervene venn
		$params
		-i $(printf '"%s" ' "${_lists_venn[@]}")
		-o "${tdirs[-1]}"
		--figtype png
	CMD
		mv "${tdirs[-1]}/Intervene_venn.png" "$outfile.png"
	CMD

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c intervene -v -b -i $threads -a cmd1
	fi

	return 0
}
