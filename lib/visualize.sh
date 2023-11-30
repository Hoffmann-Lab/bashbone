#! /usr/bin/env bash
# (c) Konstantin Riege

function visualize::venn(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-l <idfiles>  | array of bed or id files
			-n <names>    | array of
			-o <outfile>  | path to prefix
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads tmpdir="${TMPDIR:-/tmp}"
	declare -n _lists_venn _names_venn
	while getopts 'S:s:t:l:b:n:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			l)	_lists_venn=$OPTARG;;
			b)	_beds_venn=$OPTARG;;
			n)	_names_venn=$OPTARG;;
			o)	((++mandatory)); outfile="$OPTARG"; mkdir -p "$(dirname "$outfile")";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "plotting venn diagrams"

	declare -a tdirs cmd1
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
