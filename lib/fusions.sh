#! /usr/bin/env bash
# (c) Konstantin Riege

fusions::starfusion(){
	odir=fusions

	#as of v1.8 and still true for v1.9, STAR 2.7.2b is requiered

	#NEEDS CTAT_genome_lib

	# running kickoff mode
	for bam in mapped; do
		odir="$odir/($basename $bam .bam)"
		mkdir -p "$odir"
		while [[ ${bam%.*} != $bam ]];
			do bam=${bam%.*}
			[[ -e $bam.bam ]] && f=$bam.bam
		done
		bam=$f
		ln -sfnr $bam $odir/Aligned.out.bam
		find $(dirname $bam) -maxdepth 1 -regextype posix-extended -regex ".*$(basename $bam .bam)\.(_|SJ|Log|Chi).*" -exec bash -c "ln -sfnr {} $odir/\$(basename {} | sed 's/$(basename $bam .bam)\.//')" \;

		#propably, this is not necessary, since STAR-Fusion detects previously finished pipeline steps an we can go on with "from scratch"
		commander::makecmd -a cmd -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			cd "$odir"
		CMD
			STAR-Fusion
			--genome_lib_dir /misc/paras/data/genomes/CTAT_genome_lib_GRCh38_gencode_v33_Apr062020
			--CPU $THREADS
			-J $odir/Chimeric.out.junction
			--output_dir "$odir"
			--FusionInspector validate
			--examine_coding_effect
		CMD

	done

	# running from scratch
	odir="$odir/($basename $i)"
	mkdir -p "$odir"
	declare -a cmd=()
	# export CTAT_GENOME_LIB=/misc/paras/data/genomes/CTAT_genome_lib_GRCh38_gencode_v33_Apr062020
	commander::makecmd -a cmd -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		cd "$odir"
	CMD
		STAR-Fusion
		--genome_lib_dir /misc/paras/data/genomes/CTAT_genome_lib_GRCh38_gencode_v33_Apr062020
		--CPU $THREADS
		--left_fq "$i"
		--right_fq "$j"
		--output_dir "$odir"
		--FusionInspector validate
		--examine_coding_effect
	CMD

	#todo filter against arriba blacklist
	#since we do not limit alignMatesGapMax/alignIntronMax, output may contain commonly observed read-through fusions due to polymerase misses STOP
	return 0
}

fusions::arriba(){
	# compatible with CTAT_genome_lib

	# Detection of chimeric reads must be enabled in STAR by specifying the parameter --chimSegmentMin.
	# In addition, the parameter --chimOutType WithinBAM must be specified to cause STAR to report chimeric reads as supplementary alignments in
	# the main output file Aligned.out.sam. Old versions of STAR (or when STAR is run with --chimOutType SeparateSAMold) wrote supplementary
	# alignments to a separate file named Chimeric.out.sam. Arriba is compatible with this mode of use (see parameter -c).

	# Arriba STAR kickoff needs WithinBAM SoftClip or SeparateSAMold instead of Junctions (STAR-Fusion) <- try outputting both
	# related question (Junctions WithinBAM) @ https://github.com/alexdobin/STAR/issues/685
	# you are using --chimMultimapNmax 10, which is a "new" chimeric detection algorithm, allowing detection of multimapping chimeras.
	# Presently, it can only output Chimeric.out.junction file and cannot output chimeric alignments within the BAM file.
	# You can use --chimOutType Junctions WithinBAM with --chimMultimapNmax 0 ("old" unique chimeric detection) and without --peOverlap* options.

	# additional or differently setupped paremeter compared to STAR-Fusion
	# --chimScoreJunctionNonGTAG 0 (-4)

	#Arriba compared to other tools:
	# Most STAR-based fusion detection tools only consider chimeric alignments as evidence for gene fusions and are blind to focal deletions,
	# hence. As a workaround, these tools recommend reducing the value of the parameter --alignIntronMax. But this impairs the quality of alignment,
	# because it reduces the scope that STAR searches to find a spliced alignment.

	conda activate py2
	db=$CONDA_PREFIX/var/lib/arriba/grch38_gencode.blacklist

	# -x path/to/bam (WithinBAM) only or -x and -c path/to/Chimeric.out.sam (SeparateSAMold)
	# optinal -s yes|no|reverse (strandness 1|0|2)
	STAR ... | arriba \
    -a $genome -g $gtf -b /path/to/blacklist.tsv.gz \
    -x path/to/bam -c path/to/Chimeric.out.sam \
    -o fusions.tsv -O fusions.discarded.tsv \
    -T -P
    -s $( case $strandness in 0) echo "no";; 1) echo "yes";; *) echo "reverse";; esac)
	-F $fragmentsize

	return 0
}

fusions::fusioncatcher(){
	# for current v1.20, STAR 2.7.2b is requiered

	#Please run download-human-db.sh to download the Human Ensembl database to
	#${FC_DB_PATH}=/misc/paras/data/programs/rippchen/conda/envs/fusioncatcher/share/fusioncatcher-1.20/db/current
	#If you have a custom database, please use "fusioncatcher" with the option "--data".

	return 0
}

