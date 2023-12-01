#! /usr/bin/env bash
# (c) Konstantin Riege

### arriba 1.2                   ### arriba 2                   ### star-fusion (fusioninspector) ### fusioncatcher
# 1 #gene1                       1  #gene1                      1  #FusionName                    1  Gene_1_symbol(5end_fusion_partner)
# 2 gene2                        2  gene2                       2  JunctionReadCount              2  Gene_2_symbol(3end_fusion_partner)
# 3 strand1(gene/fusion)         3  strand1(gene/fusion)        3  SpanningFragCount              3  Fusion_description
# 4 strand2(gene/fusion)         4  strand2(gene/fusion)        4  est_J                          4  Counts_of_common_mapping_reads
# 5 breakpoint1                  5  breakpoint1                 5  est_S                          5  Spanning_pairs
# 6 breakpoint2                  6  breakpoint2                 6  LeftGene                       6  Spanning_unique_reads
# 7 site1                        7  site1                       7  LeftLocalBreakpoint            7  Longest_anchor_found
# 8 site2                        8  site2                       8  LeftBreakpoint                 8  Fusion_finding_method
# 9 type                         9  type                        9  RightGene                      9  Fusion_point_for_gene_1(5end_fusion_partner)
# 10 direction1                  10 split_reads1                10 RightLocalBreakpoint           10 Fusion_point_for_gene_2(3end_fusion_partner)
# 11 direction2                  11 split_reads2                11 RightBreakpoint                11 Gene_1_id(5end_fusion_partner)
# 12 split_reads1                12 discordant_mates            12 SpliceType                     12 Gene_2_id(3end_fusion_partner)
# 13 split_reads2                13 coverage1                   13 LargeAnchorSupport             13 Exon_1_id(5end_fusion_partner)
# 14 discordant_mates            14 coverage2                   14 NumCounterFusionLeft           14 Exon_2_id(3end_fusion_partner)
# 15 coverage1                   15 confidence                  15 NumCounterFusionRight          15 Fusion_sequence
# 16 coverage2                   16 reading_frame               16 FAR_left                       16 Predicted_effect
# 17 confidence                  17 tags                        17 FAR_right
# 18 closest_genomic_breakpoint1 18 retained_protein_domains    18 LeftBreakDinuc
# 19 closest_genomic_breakpoint2 19 closest_genomic_breakpoint1 19 LeftBreakEntropy
# 20 filters                     20 closest_genomic_breakpoint2 20 RightBreakDinuc
# 21 fusion_transcript           21 gene_id1                    21 RightBreakEntropy
# 22 reading_frame               22 gene_id2                    22 FFPM
# 23 peptide_sequence            23 transcript_id1              23 annots
# 24 read_identifiers            24 transcript_id2              24 CDS_LEFT_ID
#                                25 direction1                  25 CDS_LEFT_RANGE
#                                26 direction2                  26 CDS_RIGHT_ID
#                                27 filters                     27 CDS_RIGHT_RANGE
#                                28 fusion_transcript           28 PROT_FUSION_TYPE
#                                29 peptide_sequence            29 FUSION_MODEL
#                                30 read_identifiers            30 FUSION_CDS
#                                                               31 FUSION_TRANSL
#                                                               32 PFAM_LEFT
#                                                               33 PFAM_RIGHT

# -> difference supporting read counts due to multimapped vs unique mapped (may change with arriba 2 able to use multimapped junctions)

# star-fusion to arriba: 1 | sed 's/--/\t/' (~1,2), 2 (~12+13), 3 (~14), 8 | cut -d ':' -f 1,2 (~5), 11 | cut -d ':' -f 1,2 (~6), 28 (~22)
# fusioncatcher to arriba: 1 (~1), 2 (~2), 6 (~12+13), 5 (~14), 9 | cut -d ':' -f 1,2 (~5), 10 | cut -d ':' -f 1,2 (~6), 16 (~22)
# -> 6 (different to JR + SR)

# f="$fusionsdir/fusioncatcher/$(ls "$fusionsdir/fusioncatcher" | awk -v s=$s 's~$1')/final-list_candidate-fusion-genes.txt"
# full: tail -n +2 "$f" | awk -F '\t' -v OFS='\t' -v s=$s -v t=$t '{gsub(":(+|-)$","",$9); gsub(":(+|-)$","",$10); if($6+$5>=10){print s,t,$1,$2,$6,$5,"chr"$9,"chr"$10,$16}}'
# filtered: tail -n +2 "$f" | awk -F '\t' -v OFS='\t' -v s=$s -v t=$t '{gsub(":(+|-)$","",$9); gsub(":(+|-)$","",$10); print s,t,$1,$2,$6,$5,"chr"$9,"chr"$10,$16}'

function fusions::starfusion(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-5 <skip>         | true/false md5sums, indexing respectively
			-t <threads>      | number of
			-g <genome>       | path to
			-o <outdir>       | path to
			-1 <fastq1>       | array of
			-2 <fastq2>       | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome outdir
	declare -n _fq1_starfusion _fq2_starfusion
	while getopts 'S:s:5:t:g:o:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/starfusion"; mkdir -p "$outdir";;
			1)	((++mandatory)); _fq1_starfusion=$OPTARG;;
			2)	_fq2_starfusion=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	commander::printinfo "detecting gene fusions star-fusion"

	# as of v1.8 and still true for v1.9, STAR 2.7.2b is requiered by plug-n-play CTAT_genome_lib, unless index is re-created
	# STAR-Fusion can use multimapped fusions and uses the following STAR parameters
	# -> incompatible with arriba (--chimOutType WithinBAM SoftClip), so no kickstart possible
	# -> incompatibe with fusioncatcher, since this pipeline comes with an own genome
	# --outSAMattrRGline ID:GRPundef
	# --alignInsertionFlush Right
	# --alignSplicedMateMapLmin 30: STAR default is 0.66 as fraction of read and not a fixed length
	# --alignSJDBoverhangMin 10: default in alignment::star
	# --alignMatesGapMax/alignIntronMax 100000: to avoid readthru fusions within 100kdue to polymerase misses STOP (see insertsize)
	# -> apply arriba blacklist filter

	# outputs: # abridged -> shortened version lacking the voluminous read identities
	# FusionInspector-validate/
	# finspector.FusionInspector.fusions.tsv # the 'in silico validated' fusion predictions from star-fusion.fusion_predictions.tsv
	# -> finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect
	# files for IGV:
	# finspector.(fa|bed|gtf)  # reference transcript structure and annotations in BED or GTF format
	# finspector.consolidated.bam # reads aligned to the fusion contigs
	# finspector.junction_reads.bam # junction / split-reads supporting fusions
	# finspector.spanning_reads.bam # fusion spanning fragment evidence

	# star array is implicitly defined globally in alignment::star function
	# to be able to use a temporary star array here, declare it explicitly i.e. locally
	declare -a mapper star
	declare -A strandness
	local params=""
	params+=" --chimSegmentMin 12"
	params+=" --chimJunctionOverhangMin 8"
	params+=" --chimOutJunctionFormat 1"
	params+=" --alignSJstitchMismatchNmax 5 -1 5 5"
	params+=" --chimMultimapScoreRange 3"
	params+=" --chimScoreJunctionNonGTAG -4"
	params+=" --chimMultimapNmax 20"
	params+=" --chimNonchimScoreDropMin 10"
	params+=" --peOverlapNbasesMin 12"
	params+=" --peOverlapMMp 0.1"
	params+=" --alignSplicedMateMapLminOverLmate 0"
	params+=" --alignSplicedMateMapLmin 30"

	# decouple star-fusion from star to check for genome index versions and to be able to use a recent version like arriba does
	alignment::star \
		-S false \
		-s $skip \
		-5 $skipmd5 \
		-1 _fq1_starfusion \
		-2 _fq2_starfusion \
		-i 100000 \
		-o "$outdir" \
		-t $threads \
		-g "$genome" \
		-x $genome.star.idx \
		-r mapper \
		-P "$params"

	local i j odir b e
	declare -a cmd1 cmd2
	for i in "${!_fq1_starfusion[@]}"; do
		helper::basename -f "${_fq1_starfusion[$i]}" -o b -e e
		j="$outdir/star/$b.Chimeric.out.junction"
		odir="$outdir/$b"
		mkdir -p $odir

		[[ -s "$j" ]] && params="-J '$j'"
		if [[ ${_fq2_starfusion[$i]} ]]; then
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "$odir"
			CMD
				STAR-Fusion
				$params
				--genome_lib_dir "$(dirname "$genome")"
				--CPU $threads
				--left_fq "${_fq1_starfusion[$i]}"
				--right_fq "${_fq2_starfusion[$i]}"
				--output_dir "$odir"
				--FusionInspector validate
				--examine_coding_effect
			CMD
		else
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "$odir"
			CMD
				STAR-Fusion
				$params
				--genome_lib_dir "$(dirname "$genome")"
				--CPU $threads
				--left_fq "${_fq1_starfusion[$i]}"
				--output_dir "$odir"
				--FusionInspector validate
				--examine_coding_effect
			CMD
		fi

		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			tail -n +2 "$odir/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect"
		CMD
			perl -F'\t' -lane '
				BEGIN{
					print "gene1\tgene2\tjunction_count\tspanning_pairs\tbreakpoint1\tbreakpoint2\teffect"
				}
				$F[27]="NA" if $F[27] eq ".";
				@g=split/--/,$F[0]; $F[7]=~s/:[+-]$//g; $F[10]=~s/:[+-]$//g;
				for my $g1 (split/,/,$g[0]){ $g1=~s/\(\d+\)$//;
					for my $g2 (split/,/,$g[1]){ $g2=~s/\(\d+\)$//;
						next if $g1 eq $g2;
						print $g1 lt $g2 ? join("\t",($g1,$g2,$F[1],$F[2],$F[7],$F[10],$F[27])) : join("\t",($g2,$g1,$F[1],$F[2],$F[10],$F[7],$F[27]));
					}
				}
			'
		CMD
			tee >(awk -F '\t' 'NR==1 || \$3+\$4>=10' > "$outdir/$b.tsv") > "$outdir/$b.full.tsv" | cat
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c starfusion -v -b -i 1 -a cmd1
	fi

	return 0
}

function fusions::arriba(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>      | true/false return
			-s <softskip>      | true/false only print commands
			-5 <skip>          | true/false md5sums, indexing respectively
			-t <threads>       | number of
			-g <genome>        | path to
			-v <genomeversion> | hg19/hg38/mm10
			-a <gtf>           | path to
			-o <outdir>        | path to
			-f <size>          | assumed mean fragment
			-1 <fastq1>        | array of
			-2 <fastq2>        | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome genomeversion gtf outdir tmpdir="${TMPDIR:-/tmp}" fragmentsize
	declare -n _fq1_arriba _fq2_arriba
	while getopts 'S:s:5:t:g:v:a:o:f:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			v)	genomeversion="$OPTARG";;
			a)	((++mandatory)); gtf="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/arriba"; mkdir -p "$outdir";;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			1)	((++mandatory)); _fq1_arriba=$OPTARG;;
			2)	_fq2_arriba=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage
	commander::printinfo "detecting gene fusions arriba"

	# Arriba STAR kickoff needs WithinBAM SoftClip or SeparateSAMold instead of Junctions (STAR-Fusion)
	# returning both in one STAR command impossible: related question (Junctions WithinBAM) @ https://github.com/alexdobin/STAR/issues/685
	# --chimMultimapNmax and --peOverlap* options incompatible with WithinBAM or SeparateSAMold
	# Most STAR-based fusion detection tools only consider chimeric alignments as evidence for gene fusions and are blind to focal deletions,
	# hence. As a workaround, these tools recommend reducing the value of the parameter --alignIntronMax. But this impairs the quality of alignment,
	# because it reduces the scope that STAR searches to find a spliced alignment.
	# -> alignment::star sets alignIntronMax = alignMatesGapMax = 200000 or user value by -i -> set to star default (0)

	# for arriba v2.x : star > 2.7.6 is required to allow for multimapped fusions
	# declare -a cmdchk=("conda list -f arriba | tail -1 | awk '{print \$2}' | cut -d '.' -f 1")
	declare -a cmdchk=("arriba -h | grep -m 1 Version | awk '{print \$2}' | cut -d '.' -f 1")
	local version=$(commander::runcmd -c arriba -a cmdchk) params
	if [[ $version -lt 2 ]]; then
		params+=" --outFilterMultimapNmax 1"
		params+=" --outFilterMismatchNmax 3"
		params+=" --chimOutJunctionFormat 1"
		params+=" --chimScoreMin 1"
		params+=" --chimScoreSeparation 1"
		params+=" --chimSegmentMin 10"
		params+=" --chimJunctionOverhangMin 10"
		params+=" --chimNonchimScoreDropMin 10"
		params+=" --alignSJstitchMismatchNmax 5 -1 5 5"
		params+=" --chimScoreDropMax 30"
		params+=" --chimOutType WithinBAM SoftClip"
		params+=" --chimScoreJunctionNonGTAG 0"
		params+=" --chimSegmentReadGapMax 3"
	else
		params+=" --outFilterMultimapNmax 50"
		params+=" --peOverlapNbasesMin 10"
		params+=" --alignSplicedMateMapLminOverLmate 0.5"
		params+=" --alignSJstitchMismatchNmax 5 -1 5 5"
		params+=" --chimSegmentMin 10"
		params+=" --chimOutType WithinBAM HardClip"
		params+=" --chimJunctionOverhangMin 10"
		params+=" --chimScoreDropMax 30"
		params+=" --chimScoreJunctionNonGTAG 0"
		params+=" --chimScoreSeparation 1"
		params+=" --chimSegmentReadGapMax 3"
		params+=" --chimMultimapNmax 50"
	fi

	# star array is implicitly defined globally in alignment::star function
	# to be able to use a temporary star array here, declare it explicitly i.e. locally
	declare -a mapper star
	declare -A strandness

	alignment::star \
		-S false \
		-s $skip \
		-5 $skipmd5 \
		-1 _fq1_arriba \
		-2 _fq2_arriba \
		-i 0 \
		-o "$outdir" \
		-t $threads \
		-g "$genome" \
		-x $genome.star.idx \
		-r mapper \
		-P "$params"

	alignment::postprocess \
		-S false \
		-s $skip \
		-j sort \
		-t $threads \
		-o "$outdir" \
		-r mapper

	alignment::postprocess \
		-S false \
		-s $skip \
		-j index \
		-t $threads \
		-o "$outdir" \
		-r mapper

	alignment::inferstrandness \
		-S false \
		-s $skip \
		-t $threads \
		-r mapper \
		-x strandness \
		-g "$gtf"

	cmdchk=('ls "$CONDA_PREFIX/var/lib/arriba/blacklist_"*'${genomeversion}'_*.gz 2> /dev/null || mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.arriba')
	local params="-b '$(commander::runcmd -c arriba -a cmdchk)'"
		# for arriba 2.x -T/-P/-I are now default and -T -T/-P -P/-I -I is now -X (ie report sequences and read ides in the discarded fusions file too)
	[[ $version -lt 2 ]] && params+=" -T -P -I"

	local m f odir b
	declare -a cmd1 cmd2
	for f in "${star[@]}"; do
		b=$(basename "$f" .sorted.bam)
		odir="$outdir/$b"
		mkdir -p "$odir"
		commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
			arriba
				-a "$genome"
				-g "$gtf"
				-x "$f"
				-o "$odir/$b.tsv"
				-O "$odir/$b.discarded.tsv"
				-s $(case ${strandness["$f"]} in 0) echo "no";; 1) echo "yes";; 2) echo "reverse";; *) echo "?";; esac)
				-F $fragmentsize
				$params
		CMD
		# -s may be set to auto, to skip alignment::inferstrandness detour and star version conflicts by serate star conda env
		# but I encountered errors yelling unable to determine strandness AND alignment::star checks for star index version
		# -k "\$(ls '\$CONDA_PREFIX/var/lib/arriba/known_fusions_${genomeversion}_'*.gz)"
		# -p "\$(ls '\$CONDA_PREFIX/var/lib/arriba/protein_domains_${genomeversion}_'*.gff3)"

		if [[ $version -lt 2 ]]; then
			# NOTE: arriba v1 also removes chr prefix. v2.0 does not
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				tail -n +2 "$odir/$b.tsv"
			CMD
				perl -F'\t' -lane '
					BEGIN{
						print "gene1\tgene2\tjunction_count\tspanning_pairs\tbreakpoint1\tbreakpoint2\teffect"
					}
					$F[21]="NA" if $F[21] eq ".";
					$F[4]="chr$F[4]"; $F[5]="chr$F[5]";
					for my $g1 (split/,/,$F[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$F[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($g1,$g2,$F[11]+$F[12],$F[13],$F[4],$F[5],$F[21])) : join("\t",($g2,$g1,$F[11]+$F[12],$F[13],$F[5],$F[4],$F[21]));
						}
					}
				'
			CMD
				tee >(awk -F '\t' 'NR==1 || \$3+\$4>=10' > "$outdir/$b.tsv") > "$outdir/$b.full.tsv" | cat
			CMD
		else
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				tail -n +2 "$odir/$b.tsv"
			CMD
				perl -F'\t' -lane '
					BEGIN{
						print "gene1\tgene2\tjunction_count\tspanning_pairs\tbreakpoint1\tbreakpoint2\teffect"
					}
					$F[15]="NA" if $F[15] eq ".";
					for my $g1 (split/,/,$F[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$F[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($g1,$g2,$F[9]+$F[10],$F[11],$F[4],$F[5],$F[15])) : join("\t",($g2,$g1,$F[9]+$F[10],$F[11],$F[5],$F[4],$F[15]));
						}
					}
				'
			CMD
				tee >(awk -F '\t' 'NR==1 || \$3+\$4>=10' > "$outdir/$b.tsv") > "$outdir/$b.full.tsv" | cat
			CMD
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c arriba -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function fusions::join(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-i <fusionsdir> | path to
			-o <outdir>     | path to
			-d <tool>       | name identical to subdir in fusionsdir. parameter can be used multiple times
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads fusionsdir outdir tmpdir="${TMPDIR:-/tmp}"
	declare -a tools
	while getopts 'S:s:t:i:o:d:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); fusionsdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir" ;;
			d)	tools+=("$OPTARG");;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage


	echo -e "sample tool fusion gene1 gene2 junction_count spanning_pairs breakpoint1 breakpoint2 effect" | sed 's/ /\t/g' > "$outdir/fusions.tsv"
	head -1 "$outdir/fusions.tsv" > "$outdir/fusions.full.tsv"
	local s farr fsfus
	while read -r s; do
		for tool in "${tools[@]}"; do
			case $tool in
				arriba) farr=$(find -L "$fusionsdir/arriba/" -type f -name "$s*.fusions.tsv" | grep .);;
				starfusion) fsfus=$(realpath -s "$fusionsdir/starfusion/$s/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect");;
				*) commander::printerr "$tool not supported for joining gene fusion predictions"; return 1;;
				# f="$fusionsdir/fusioncatcher/$(ls "$fusionsdir/fusioncatcher" | awk -v s=$s 's~$1')/final-list_candidate-fusion-genes.txt"
				# full: tail -n +2 "$f" | awk -F '\t' -v OFS='\t' -v s=$s -v t=$t '{gsub(":(+|-)$","",$9); gsub(":(+|-)$","",$10); if($6+$5>=10){print s,t,$1,$2,$6,$5,"chr"$9,"chr"$10,$16}}'
				# filtered: tail -n +2 "$f" | awk -F '\t' -v OFS='\t' -v s=$s -v t=$t '{gsub(":(+|-)$","",$9); gsub(":(+|-)$","",$10); print s,t,$1,$2,$6,$5,"chr"$9,"chr"$10,$16}'
			esac
		done

		if [[ $(head -1 "$farr" | tr '\t' '\n' | wc -l) -lt 30 ]]; then
			# arriba v1
			{
				tail -n +2 "$farr" | perl -slane '
					$F[21]="NA" if $F[21] eq ".";
					$F[4]="chr$F[4]"; $F[5]="chr$F[5]";
					for my $g1 (split/,/,$F[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$F[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($s,"arriba",$g1,$g2,$F[11]+$F[12],$F[13],$F[4],$F[5],$F[21])) : join("\t",($g2,$g1,$F[11]+$F[12],$F[13],$F[5],$F[4],$F[21]));
						}
					}
				' -- -s=$s
				tail -n +2 "$fsfus" | perl -slane '
					$F[27]="NA" if $F[27] eq ".";
					@g=split/--/,$F[0]; $F[7]=~s/:[+-]$//g; $F[10]=~s/:[+-]$//g;
					for my $g1 (split/,/,$g[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$g[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($s,"starfusion",$g1,$g2,$F[1],$F[2],$F[7],$F[10],$F[27])) : join("\t",($s,"starfusion",$g2,$g1,$F[1],$F[2],$F[10],$F[7],$F[27]));
						}
					}
				' -- -s=$s
			} | helper::sort -t $threads -k3,3 -k4,4 -k7,7 -k8,8 | tee >(awk '$5+$6>=10' >> "$outdir/fusions.tsv") >> "$outdir/fusions.full.tsv"
		else
			# arriba v2
			{
				tail -n +2 "$farr" | perl -slane '
					$F[15]="NA" if $F[15] eq ".";
					for my $g1 (split/,/,$F[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$F[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($s,"arriba",$g1,$g2,$F[9]+$F[10],$F[11],$F[4],$F[5],$F[15])) : join("\t",($s,"arriba",$g2,$g1,$F[9]+$F[10],$F[11],$F[5],$F[4],$F[15]));
						}
					}
				' -- -s=$s
				tail -n +2 "$fsfus" | perl -slane '
					$F[27]="NA" if $F[27] eq ".";
					@g=split/--/,$F[0]; $F[7]=~s/:[+-]$//g; $F[10]=~s/:[+-]$//g;
					for my $g1 (split/,/,$g[0]){ $g1=~s/\(\d+\)$//;
						for my $g2 (split/,/,$g[1]){ $g2=~s/\(\d+\)$//;
							next if $g1 eq $g2;
							print $g1 lt $g2 ? join("\t",($s,"starfusion",$g1,$g2,$F[1],$F[2],$F[7],$F[10],$F[27])) : join("\t",($s,"starfusion",$g2,$g1,$F[1],$F[2],$F[10],$F[7],$F[27]));
						}
					}
				' -- -s=$s
			} | helper::sort -t $threads -k3,3 -k4,4 -k7,7 -k8,8 | tee >(awk '$5+$6>=10' >> "$outdir/fusions.tsv") >> "$outdir/fusions.full.tsv"
		fi
	done < <(find -L "$fusionsdir/arriba/star/" -type f -name "*.sorted.bam" -exec basename {} .sorted.bam \; | grep .)

	head -1 "$outdir/fusions.tsv" > "$outdir/fusions.merged.tsv"
	tail -n +2 "$outdir/fusions.tsv" | perl -lane 'if($F[0] eq $o[0] && $F[2] eq $o[2] && $F[3] eq $o[3]){@o=@F if $F[4]+$F[5] > $o[4]+$o[5]}else{print join"\t",@o if @o; @o=@F}END{print join"\t",@o}' >> "$outdir/fusions.merged.tsv"
	head -1 "$outdir/fusions.full.tsv" > "$outdir/fusions.full.merged.tsv"
	tail -n +2 "$outdir/fusions.full.tsv" | perl -lane 'if($F[0] eq $o[0] && $F[2] eq $o[2] && $F[3] eq $o[3]){@o=@F if $F[4]+$F[5] > $o[4]+$o[5]}else{print join"\t",@o if @o; @o=@F}END{print join"\t",@o}' >> "$outdir/fusions.full.merged.tsv"

	return 0
}
