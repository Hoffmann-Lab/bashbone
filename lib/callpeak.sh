#! /usr/bin/env bash
# (c) Konstantin Riege

callpeak::mkreplicates() { #for peak calling only
	[[ $norepl ]] && return 0
	[[ ${#ridx[@]} -gt 0 ]] && echo ":INFO: generating pools" || echo ":INFO: generating pseudo-replicates"

	cmd1=()
	cmd2=()
	if [[ $ridx ]]; then # pool replicates: m[N1 N2 T1 T2 R1 R2] -> m[N1 N2 T1 T2 R1 R2 P1 P2]
		for r in ${mapper[@]}; do
			declare -n m=$r
			j=$[${#m[@]}-1]
			mkdir -p $tmpdir/$r
			for i in ${!nidx[@]}; do
				tf=${m[${tidx[$i]}]}
				rf=${m[${ridx[$i]}]}
				o=${tf%.*}+$(basename ${rf%.*}).pseudopool.bam
				tmp1=$tmpdir/$r/$(basename $tf)
				tmp2=$tmpdir/$r/$(basename $rf)
				# does not work in rare cases - reason unknown
				#cmd1+=("samtools view -@ $threads -b -s 0.5 $tf | samtools reheader -P <(echo) - > $tmp\0")
				#cmd2+=("samtools view -@ $threads -b -s 0.5 $rf | samtools reheader -P <(echo) - >> $tmp\0")
				#cmd3+=("samtools reheader -P <(samtools view -H $tf) $tmp > $o\0")
				cmd1+=("samtools view -@ $threads -b -s 0.5 $tf > $tmp1\0")
				cmd1+=("samtools view -@ $threads -b -s 0.5 $rf > $tmp2\0")
				cmd2+=("samtools cat $tmp1 $tmp2 > $o\0")
				m+=($o)
				pidx+=($((++j)))
			done
		done
		
		if [[ $nridx ]]; then
			# create fully pooled bams and manipulate indices for nidx based peak calling
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 P1 P2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 P1 P2 NFP1 FP1 NFP2 FP2]
			# n  1 2   3 4   11 13
			# nr 3 4
			# t  5 6   5 6   5  6
			# r  7 8   7 8   7  8
			# p  9 10  9 10  12 14
			# -> idr: 1 vs 9 + 1 vs 5 + 1 vs 7
			# 			2 vs 10 + 2 vs 6 + 2 vs 8
			# 		  3 vs 9 + 3 vs 5 + 3 vs 7
			# 			4 vs 10 + 4 vs 6 + 4 vs 8
			# 		  11 vs 12 + 11 vs 5 + 11 vs 7
			# 			13 vs 14 + 13 vs 6 + 13 vs 8

			nidx+=(${nridx[@]})
			tidx+=(${tidx[@]} ${tidx[@]})
			ridx+=(${ridx[@]} ${ridx[@]})
			pidx+=(${pidx[@]})
			addindex=true
			for r in ${mapper[@]}; do
				declare -n m=$r
				j=$[${#m[@]}-1]
				mkdir -p $tmpdir/$r
				for i in ${!nridx[@]}; do
					nf=${m[${nidx[$i]}]}
					nrf=${m[${nridx[$i]}]}
					o=${nf%.*}+$(basename ${nrf%.*}).pool.bam
					cmd2+=("samtools cat $nf $nrf > $o\0")
					m+=($o)
					$addindex && nidx+=($((++j)))

					o=${tf%.*}+$(basename ${rf%.*}).pool.bam
					cmd2+=("samtools cat $tf $rf > $o\0")
					m+=($o)
					$addindex && pidx+=($((++j)))
				done
				addindex=false
			done
		fi

		thr=1
	else # make pseudo-replicates from pseudo-pool: # m[N1 N2 P1 P2] -> m[N1 N2 P1 P2 T1 R1 T2 R2]
		instances=0
		for r in ${mapper[@]}; do
			declare -n m=$r
			instances=$[instances+${#m[@]}]
		done
		[[ thr=$[threads/instances] -eq 0 ]] && thr=1

		for r in ${mapper[@]}; do
			declare -n m=$r
			j=$[${#m[@]}-1]
			mkdir -p $tmpdir/$r
			for i in ${!nidx[@]}; do
				pf=${m[${pidx[$i]}]}
				tidx+=($((++j)))
				ridx+=($((++j)))
				o=${pf%.*}.pseudoreplicate # greped during mappingstats
				tmp=$tmpdir/$r/$(basename $o)
				# -u not equals --output-fmt SAM => a small compression level reduces stream amount which finally make decompression faster (optimum -l 3)
				# tested with 223G bam (~1m30 vs ~1m) and 650M bam (~1m30 vs ~1m)
				cmd1+=("samtools collate -@ $thr -n $thr -u -l 3 -O $pf $tmp | samtools view -@ $thr - | split --numeric-suffixes=1 --additional-suffix=.bam -a 1 -l $[($(samtools flagstat -@ $thr $pf | head -1 | cut -d ' ' -f 1)+1)/2] --filter='cat <(samtools view -H $pf) - | samtools view -@ $thr -b > \$FILE' - $o\0")
				m+=(${o}1.bam)
				m+=(${o}2.bam)
			done
		done

		thr=$threads
	fi

	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	if [[ ! $Srep ]]; then
		{	echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P $thr -I {} bash -c {} && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}
		} || return 1
	fi

	return 0
}

callpeak::macs() {
	[[ $nomacs ]] && return 0
	echo ":INFO: peak calling - macs"

	cmd0=()
	cmd1=()
	cmd2=()
	cmd3=()
	cmd4=()
	for r in ${mapper[@]}; do
		declare -n m=$r
		odir=$outdir/peaks/$r/macs
		tmp=$tmpdir/$r
		mkdir -p $odir $tmp
		poolidr=()
		for i in ${!nidx[@]}; do
			on=$(basename ${m[${nidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			on=${on%.*}
			op=$(basename ${m[${pidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			op=${op%.*}
			ot=$(basename ${m[${tidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			ot=${ot%.*}
			or=$(basename ${m[${ridx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			or=${or%.*}

			# TODO check newer versions for fix - still failes as of Jan 2019 v2.1.2 
			# to avoid errors in MACS2.IO.Parser.BAMParser
			cmd0+=("samtools view -h ${m[${pidx[$i]}]} > ${m[${pidx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${tidx[$i]}]} > ${m[${tidx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${ridx[$i]}]} > ${m[${ridx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${nidx[$i]}]} > ${m[${nidx[$i]}]}.sam\0")

			cmd1+=("macs2 callpeak -t ${m[${pidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$op.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1\0")
			cmd1+=("macs2 callpeak -t ${m[${pidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$op.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift 0 --extsize $fragmentsize\0")
			cmd1+=("macs2 callpeak -t ${m[${tidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$ot.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1\0")
			cmd1+=("macs2 callpeak -t ${m[${tidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$ot.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift -0 --extsize $fragmentsize\0")
			cmd1+=("macs2 callpeak -t ${m[${ridx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$or.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1\0")
			cmd1+=("macs2 callpeak -t ${m[${ridx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $on-vs-$or.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift -0 --extsize $fragmentsize\0")

			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$on-vs-$op.model_peaks.narrowPeak $odir/$on-vs-$op.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$on-vs-$op.merged_peaks.narrowPeak\0")
			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$on-vs-$ot.model_peaks.narrowPeak $odir/$on-vs-$ot.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$on-vs-$ot.merged_peaks.narrowPeak\0")
			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$on-vs-$or.model_peaks.narrowPeak $odir/$on-vs-$or.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$on-vs-$or.merged_peaks.narrowPeak\0")

			cmd3+=("idr --samples $odir/$on-vs-$ot.merged_peaks.narrowPeak $odir/$on-vs-$or.merged_peaks.narrowPeak --peak-list $odir/$on-vs-$op.merged_peaks.narrowPeak --input-file-type narrowPeak --output-file $tmp/$on-vs-$op.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

			cmd4+=("cut -f 1-10 $tmp/$on-vs-$op.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$on-vs-$op.idr.narrowPeak\0")

			poolidr+=($odir/$on-vs-$op.merged_peaks) # if a normal replicate is given, then run idr on pseudopools vs pool
			if [[ ${#poolidr[@]} -eq 3 ]]; then
				o=$on-VERSUS-$(basename ${poolidr[0]} .merged_peaks)+$(basename ${poolidr[1]} .merged_peaks)
				cmd3+=("idr --samples ${poolidr[0]}.narrowPeak ${poolidr[1]}.narrowPeak --peak-list ${poolidr[2]}.narrowPeak --input-file-type narrowPeak --output-file $tmp/$o.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

				cmd4+=("cut -f 1-10 $tmp/$o.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$o.idr.narrowPeak\0")
				poolidr=()
			fi
		done
	done

	echo -n ${cmd0[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd3[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd4[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	# IDR can fail with warning: ValueError: Peak files must contain at least 20 peaks post-merge, Merged peaks were written to the output file
	if [[ ! $Smacs ]]; then
		{	echo -ne ${cmd0[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py3 && \
			echo -ne ${cmd3[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}; true && \
			conda activate py2 && \
			echo -ne ${cmd4[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}
		} || return 1
	fi

	return 0
}

callpeak::macs_m6a() {
	[[ $nomacs ]] && return 0
	echo ":INFO: peak calling - macs"

	cmd0=()
	cmd1=()
	cmd2=()
	cmd3=()
	cmd4=()
	for r in ${mapper[@]}; do
		declare -n m=$r
		odir=$outdir/peaks/$r/macs
		tmp=$tmpdir/$r
		mkdir -p $odir $tmp
		for i in ${!nidx[@]}; do
			op=$(basename ${m[${pidx[$i]}]})
			op=${op%.*}
			ot=$(basename ${m[${tidx[$i]}]})
			ot=${ot%.*}
			or=$(basename ${m[${ridx[$i]}]})
			or=${or%.*}

			# TODO check newer versions for fix - still failes as of Jan 2019 v2.1.2 
			# to avoid errors in MACS2.IO.Parser.BAMParser
			cmd0+=("samtools view -h ${m[${pidx[$i]}]} > ${m[${pidx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${tidx[$i]}]} > ${m[${tidx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${ridx[$i]}]} > ${m[${ridx[$i]}]}.sam\0")
			cmd0+=("samtools view -h ${m[${nidx[$i]}]} > ${m[${nidx[$i]}]}.sam\0")

			cmd1+=("macs2 callpeak -t ${m[${pidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $op.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --mfold 3 500 --bw $fragmentsize\0")
			cmd1+=("macs2 callpeak -t ${m[${pidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $op.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift -30 --extsize $[fragmentsize-30]\0")
			cmd1+=("macs2 callpeak -t ${m[${tidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $ot.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --mfold 3 500 --bw $fragmentsize\0")
			cmd1+=("macs2 callpeak -t ${m[${tidx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $ot.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift -30 --extsize $[fragmentsize-30]\0")
			cmd1+=("macs2 callpeak -t ${m[${ridx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $or.model --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --mfold 3 500 --bw $fragmentsize\0")
			cmd1+=("macs2 callpeak -t ${m[${ridx[$i]}]}.sam -c ${m[${nidx[$i]}]}.sam -f SAM -g hs --outdir $odir -n $or.nomodel --tempdir $tmp -B --SPMR --keep-dup all -q 0.05 --verbose 1 --nomodel --shift -30 --extsize $[fragmentsize-30]\0")

			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$op.model_peaks.narrowPeak $odir/$op.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$op.narrowPeak\0")
			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$ot.model_peaks.narrowPeak $odir/$ot.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$ot.narrowPeak\0")
			cmd2+=("bedtools merge -c 5,7,8,9 -o collapse -i <(sort -k1,1V -k2,2n -k3,3n $odir/$or.model_peaks.narrowPeak $odir/$or.nomodel_peaks.narrowPeak) | perl -M'List::Util qw(max)' -lane 'print join(\"\\\t\",@F[0..2],\"merged_peak_\".(++\$x),max(split/,/,\$F[3]),\".\",max(split/,/,\$F[4]),max(split/,/,\$F[5]),max(split/,/,\$F[6]),\"-1\")' > $odir/$or.narrowPeak\0")

			cmd3+=("idr --samples $odir/$ot.narrowPeak $odir/$or.narrowPeak --peak-list $odir/$op.narrowPeak --input-file-type narrowPeak --output-file $tmp/$op.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

			cmd4+=("cut -f 1-10 $tmp/$op.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$op.idr.narrowPeak\0")
		done
	done

	echo -n ${cmd0[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd3[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd4[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	if [[ ! $Smacs ]]; then
		{	echo -ne ${cmd0[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py3 && \
			echo -ne ${cmd3[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py2 && \
			echo -ne ${cmd4[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}
		} || return 1
	fi

	return 0
}

callpeak::gem() {
	[[ $nogem ]] && return 0
	echo ":INFO: peak calling - gem"

	if [[ ! $Sgem ]]; then
		mkdir -p $tmpdir/genome
		for r in ${mapper[@]}; do
			declare -n m=$r
			samtools view -H ${m[${pidx[0]}]} | sed -rn '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t\2/p}' > $tmpdir/genome/chr.info
			break
		done
		perl -lane 'if($_=~/^>(\S+)/){close F; open F,">'$tmpdir'/genome/$1.fa"; print F $_;}else{print F $_}; END{close F;}' $genome
	fi

	jmem=$memory
	[[ $jmem -gt $mem ]] && jmem=$mem
	jgct=$[(3+5*threads/8)/instances]
    [[ $jgct -eq 0 ]] && jgct=1
    jmgct=$[jgct/4]
    [[ $jmgct -eq 0 ]] && jmgct=1

	cmd1=()
	cmd2=()
	cmd3=()
	for r in ${mapper[@]}; do
		declare -n m=$r
		odir=$outdir/peaks/$r/gem		
		tmp=$tmpdir/$r/gem
		mkdir -p $tmp
		pooldir=()
		for i in ${!nidx[@]}; do
			on=$(basename ${m[${nidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			on=${on%.*}
			op=$(basename ${m[${pidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			op=${op%.*}
			ot=$(basename ${m[${tidx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			ot=${ot%.*}
			or=$(basename ${m[${ridx[$i]}]} | sed -r 's/\.(sorted|unique|rmdup)//g')
			or=${or%.*}
			
			# --k_min 4 --k_max 13 and replace GPS_events to GEM_events
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$on-vs-$op --expt ${m[${pidx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $insdir/bin/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP && cp $odir/$on-vs-$op/$on-vs-$op.GPS_events.narrowPeak $odir/$on-vs-$op.narrowPeak\0")
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$on-vs-$ot --expt ${m[${tidx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $insdir/bin/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP && cp $odir/$on-vs-$ot/$on-vs-$ot.GPS_events.narrowPeak $odir/$on-vs-$ot.narrowPeak\0")
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$on-vs-$or --expt ${m[${ridx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $insdir/bin/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP && cp $odir/$on-vs-$or/$on-vs-$or.GPS_events.narrowPeak $odir/$on-vs-$or.narrowPeak\0")

			cmd2+=("idr --samples $odir/$on-vs-$ot.narrowPeak $odir/$on-vs-$or.narrowPeak --peak-list $odir/$on-vs-$op.narrowPeak --input-file-type narrowPeak --output-file $tmp/$on-vs-$op.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

			cmd3+=("cut -f 1-10 $tmp/$on-vs-$op.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$on-vs-$op.idr.narrowPeak\0")

			poolidr+=($odir/$on-vs-$op)
			if [[ ${#poolidr[@]} -eq 3 ]]; then
				o=$on-VERSUS-$(basename ${poolidr[0]})+$(basename ${poolidr[1]})
				cmd3+=("idr --samples ${poolidr[0]}.narrowPeak ${poolidr[1]}.narrowPeak --peak-list ${poolidr[2]}.narrowPeak --input-file-type narrowPeak --output-file $tmp/$o.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

				cmd4+=("cut -f 1-10 $tmp/$o.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$o.idr.narrowPeak\0")				
				poolidr=()
			fi
		done
	done

	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd3[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	if [[ ! $Sgem ]]; then
		{	echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P 1 -I {} bash -c {} && \
			conda activate py3 && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py2 && \
			echo -ne ${cmd3[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}
		} || return 1
	fi

	return 0
}

callpeak::gem_m6a() {
	[[ $nogem ]] && return 0
	echo ":INFO: peak calling - gem"

	if [[ ! $Sgem ]]; then
		mkdir -p $tmpdir/genome
		for r in ${mapper[@]}; do
			declare -n m=$r
			samtools view -H ${m[${pidx[0]}]} | sed -rn '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t\2/p}' > $tmpdir/genome/chr.info
			break
		done
		perl -lane 'if($_=~/^>(\S+)/){close F; open F,">'$tmpdir'/genome/$1.fa"; print F $_;}else{print F $_}; END{close F;}' $genome
	fi

	jmem=$memory
	[[ $jmem -gt $mem ]] && jmem=$mem
	jgct=$[(3+5*threads/8)/instances]
    [[ $jgct -eq 0 ]] && jgct=1
    jmgct=$[jgct/4]
    [[ $jmgct -eq 0 ]] && jmgct=1

	cmd1=()
	cmd2=()
	cmd3=()
	for r in ${mapper[@]}; do
		declare -n m=$r
		odir=$outdir/peaks/$r/gem		
		tmp=$tmpdir/$r/gem
		for i in ${!nidx[@]}; do
			op=$(basename ${m[${pidx[$i]}]})
			op=${op%.*}
			ot=$(basename ${m[${tidx[$i]}]})
			ot=${ot%.*}
			or=$(basename ${m[${ridx[$i]}]})
			or=${or%.*}
			
			mkdir -p $tmp $odir/$op $odir/$ot $odir/$or
			# --k_min 4 --k_max 13 and replace GPS_events to GEM_events
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$op --expt ${m[${pidx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $RIPPCHEN/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP --smooth $[fragmentsize/2] && cp $odir/$op/$op.GPS_events.narrowPeak $odir/$op.narrowPeak\0")
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$ot --expt ${m[${tidx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $RIPPCHEN/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP --smooth $[fragmentsize/2] && cp $odir/$ot/$ot.GPS_events.narrowPeak $odir/$ot.narrowPeak\0")
			cmd1+=("gem -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jmgct -Djava.io.tmpdir=$tmpdir --t $threads --genome $tmpdir/genome --g $tmpdir/genome/chr.info --out $odir/$or --expt ${m[${ridx[$i]}]} --ctrl ${m[${nidx[$i]}]} --f SAM --nrf --d $RIPPCHEN/gem/Read_Distribution_default.txt --s 2400000000 --q $(echo 0.05 | awk '{print -log($1)/log(10)}') --outNP --smooth $[fragmentsize/2] && cp $odir/$or/$or.GPS_events.narrowPeak $odir/$or.narrowPeak\0")

			cmd2+=("idr --samples $odir/$ot.narrowPeak $odir/$or.narrowPeak --peak-list $odir/$op.narrowPeak --input-file-type narrowPeak --output-file $tmp/$op.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")

			cmd3+=("cut -f 1-10 $tmp/$op.idr | sort -k1,1V -k2,2n -k3,3n | uniq > $odir/$op.idr.narrowPeak\0")
		done
	done

	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd3[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	if [[ ! $Sgem ]]; then
		{	echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P 1 -I {} bash -c {} && \
			conda activate py3 && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py2 && \
			echo -ne ${cmd3[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {}
		} || return 1
	fi

	return 0
}

___template() {
	[[ $notemplate ]] && return 0
	echo ":INFO: template"
# (sequence-)tag means read
# cutting ends are mate pairs
# narrowpeak file format
# 1 chrom - Name of the chromosome (or contig, scaffold, etc.).
# 2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# 3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined aschromStart=0, chromEnd=100, and span the bases numbered 0-99.
# 4 name - Name given to a region (preferably unique). Use '.' if no name is assigned.
# 5 score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
# 6 strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
# 7 signalValue - Measurement of overall (usually, average) enrichment for the region.
# 8 pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
# 9 qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
# 10 peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

	cmd1=()
	cmd2=()
	cmd3=()
	for r in ${mapper[@]}; do
		declare -n m=$r #not necessary here
		mkdir -p $tmpdir/$r
		for i in ${!nidx[@]}; do
			declare -n sn=$r${nidx[$i]}
			declare -n sp=$r${pidx[$i]}
			declare -n st=$r${tidx[$i]}
			declare -n sr=$r${ridx[$i]}
			
			op=$outdir/peaks/macs/$(basename ${m[${pidx[$i]}]})
			op=${op%.*}
			ot=$outdir/peaks/macs/$(basename ${m[${tidx[$i]}]})
			ot=${ot%.*}
			or=$outdir/peaks/macs/$(basename ${m[${ridx[$i]}]})
			or=${or%.*}.narrowPeak
			tmp=$tmpdir/$r
			for i in ${!sn[@]}; do
				cmd1+=("-t ${sp[$i]} -c ${sn[$i]}\0")
				cmd1+=("-t ${st[$i]} -c ${sn[$i]}\0")
				cmd1+=("-t ${sr[$i]} -c ${sn[$i]}\0")
			done
			cmd2+=("cat $tmp/$(basename $op)*.narrowPeak > $op.narrowPeak\0")
			cmd2+=("cat $tmp/$(basename $ot)*.narrowPeak > $ot.narrowPeak\0")
			cmd2+=("cat $tmp/$(basename $or)*.narrowPeak > $or.narrowPeak\0")

			cmd3+=("idr --samples $ot.narrowPeak $or.narrowPeak --peak-list $op.narrowPeak --input-file-type narrowPeak --output-file $op.idr --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR\0")
		done
	done

	echo -n ${cmd1[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd2[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'
	echo -n ${cmd3[@]} | sed 's/\\0\s*/\n/g' | awk '$1=":CMD: "$1'

	if [[ ! $Stemplate ]]; then
		{	echo -ne ${cmd1[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py3 && \
			echo -ne ${cmd2[@]} | /usr/bin/time -v xargs -0 -P $threads -I {} bash -c {} && \
			conda activate py2
		} || return 1
	fi

	return 0
}
