#! /usr/bin/env bash
# (c) Konstantin Riege

compile::_parse(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
		usage:
			-i <path>    | installation base
			-t <threads> | number of
			-u <upgrade> | true/false conda envs
		EOF
		return 1
	}

	local OPTIND arg mandatory
	declare -n _insdir_parse _threads_parse _upgrade_parse
	while getopts 'r:s:c:i:t:u:' arg; do
		case $arg in
			r) ((++mandatory)); _insdir_parse="$OPTARG";;
			s) ((++mandatory)); _threads_parse=$OPTARG;;
			c) _upgrade_parse=$OPTARG;;
			i) ((++mandatory)); _insdir_parse="$OPTARG";;
			t) ((++mandatory)); _threads_parse=$OPTARG;;
			u) _upgrade_parse=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	return 0
}

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::tools -i "$insdir" -t $threads
	compile::conda -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads
	compile::java -i "$insdir" -t $threads
	compile::trimmomatic -i "$insdir" -t $threads
	compile::sortmerna -i "$insdir" -t $threads
	compile::segemehl -i "$insdir" -t $threads
	compile::starfusion -i "$insdir" -t $threads
	compile::preparedexseq -i "$insdir" -t $threads
	compile::revigo -i "$insdir" -t $threads
	compile::gem -i "$insdir" -t $threads
	compile::m6aviewer -i "$insdir" -t $threads
	compile::idr -i "$insdir" -t $threads
	compile::newicktopdf -i "$insdir" -t $threads
	compile::ssgsea -i "$insdir" -t $threads
	compile::bgztail -i "$insdir" -t $threads
	compile::mdless -i "$insdir" -t $threads

	return 0
}

compile::bashbone() {
	local insdir threads version src=$(dirname $(readlink -e $0))
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads "$@"
	source $src/lib/version.sh
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	return 0
}

compile::tools() {
	local insdir threads i src=$(dirname $(readlink -e $0))

	commander::printinfo "installing statically pre-compiled tools"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	mkdir -p $insdir/latest
	for i in $(find "$src/tools" -mindepth 1 -maxdepth 1 -type d); do
		cp -r "$i" "$insdir/$(basename "$i")"
		ln -sfn "$insdir/$(basename "$i")/bin" "$insdir/latest/$(basename "$i" | cut -d '-' -f 1)"
	done

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::tools -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads -u true

	return 0
}

compile::conda(){
	local tmpdir
	_cleanup::compile::conda(){
		rm -rf "$tmpdir"
	}

	local insdir threads url n bin
	commander::printinfo "installing conda"
	compile::_parse -r insdir -s threads "$@"
	url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	wget -q -O "$insdir/miniconda.sh" "$url"
	mkdir -p "$insdir/conda"
	bash "$insdir/miniconda.sh" -b -u -f -p "$insdir/conda"
	rm -f "$insdir/miniconda.sh"

	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	conda update -y conda

	# perl from conda-forge is compiled with threads, perl from bioconda not (recently removed) - thus there is an old perl-threaded version
	# perl + modules and R for scripts and R based bashbone functions (freeze latest r-base version for which all packages are available)
	# to avoid R downgrades due to modules R built versions, compile them manually below (ggpubr requires nlopt)
	# non-oracle java for m6aviewer
	# htslib htseq for prepare_dexseq and some r packages
	# bcftools vcflib vt for vcf normalization
	# ucsc-facount khmer for peak calling required effective genome size estimation
	# ghostscript for ps2pdf
	commander::printinfo "setup conda base env"
	conda install -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
		gcc_linux-64 gxx_linux-64 gfortran_linux-64 \
		glib pkg-config make automake cmake \
		bzip2 pigz pbzip2 \
		wget curl ghostscript dos2unix \
		sra-tools entrez-direct \
		datamash samtools bedtools ucsc-facount khmer \
		htslib htseq bcftools vcflib vt vcftools \
		perl perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-dbi perl-db-file perl-xml-parser perl-bioperl perl-bio-eutilities \
		java-jdk \
		nlopt r-base=4.0.2
	cpanm Switch

	tmpdir="$insdir/tmp"
	mkdir -p "$tmpdir"

	# basics
	declare -a cmd1
	commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			install.packages(c("BiocManager","devtools","codetools"),
				repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
		'
	CMD
	# bioconductor
	declare -a cmd2
	commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			BiocManager::install(c("biomaRt","BiocParallel","genefilter","DESeq2","DEXSeq","clusterProfiler","TCGAutils","TCGAbiolinks","survminer","impute","preprocessCore","GO.db","AnnotationDbi"),
				ask=F, Ncpus=$threads, clean=T, destdir="$tmpdir");
		'
	CMD
	# cran - needs to be last since WGCNA depends on bioconductor packages impute,...
	# R-Forge.r does not complaine about knapsack not being compatible with R>=4
	declare -a cmd3
	commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			install.packages(c("reshape2","WGCNA","dplyr","tidyverse","ggpubr","ggplot2","gplots","RColorBrewer","svglite","pheatmap","treemap","data.table"),
				repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
			install.packages(c("knapsack"), repos="http://R-Forge.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
		'
	CMD
	# github
	declare -a cmd4
	commander::makecmd -a cmd4 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			devtools::install_github("andymckenzie/DGCA", upgrade="never", force=T, clean=T, destdir="$tmpdir");
		'
	CMD
	commander::runcmd -t 1 -a cmd1
	commander::runcmd -t 1 -a cmd2
	commander::runcmd -t 1 -a cmd3
	commander::runcmd -t $threads -a cmd4

	commander::printinfo "conda clean up"
	conda clean -y -a
	conda deactivate

	return 0
}

compile::conda_tools() {
	local insdir threads upgrade=false tool n bin doclean=false
	declare -A envs

	compile::_parse -r insdir -s threads -c upgrade "$@"
	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r tool; do
		envs[$tool]=true
	done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

	# better do not predefine python version. if tool recipe depends on earlier version, conda installs an older or the oldest version (freebayes)

	# arriba 2.x , successor of 1.2 (arriba=1.2) has new star parameters incompatible with star < 2.7.6
	for tool in fastqc cutadapt rcorrector star bwa rseqc subread htseq arriba picard bamutil macs2 peakachu diego gatk4 freebayes varscan igv intervene raxml metilene; do
		n=${tool/=*/}
		n=${n//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} || {
			doclean=true

			commander::printinfo "setup conda $n env"
			conda create -y -n $n #python=3
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		}
		# link commonly used base binaries into env
		for bin in perl bgzip samtools bcftools bedtools vcfsamplediff; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	done
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

	# manual setup of requirements from bioconda meta.yaml (see compile::starfusion) due to non-latest installation via conda
	# note: recent star is not compatible with CTAT plug-n-play genome index as of CTAT for star-fusion v1.9 (star indexer 2.7.1a)
	tool=star-fusion
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $n env"
		conda create -y -n $n #python=3
		# propably enought: perl perl-set-intervaltree perl-carp perl-carp-assert perl-db-file perl-io-gzip perl-json-xs perl-uri \
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			perl perl-file-path perl-getopt-long perl-set-intervaltree perl-carp perl-carp-assert perl-data-dumper perl-findbin perl-db-file perl-io-gzip perl-json-xs perl-uri perl-list-moreutils perl-list-util perl-storable \
			igv-reports star gmap bowtie bbmap samtools blast
	}
	for bin in perl bgzip samtools bcftools bedtools vcfsamplediff; do
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done

	# customized env setupus

	tool=vardict
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $n env"
		conda create -y -n $n #python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool vardict-java readline=6
	}
	for bin in perl bgzip samtools bcftools bedtools vcfsamplediff; do
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done

	tool=snpeff
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $n env"
		conda create -y -n $n #python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool snpsift
	}
	for bin in perl bgzip samtools bcftools bedtools vcfsamplediff; do
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done

	tool=platypus-variant
	n=platypus
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $n env"
		conda create -y -n $n #python=2
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
	}
	for bin in perl bgzip samtools bcftools bedtools vcfsamplediff; do
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done

	# this is a pipeline itself with own genome and databases and thus will not be part of bashbone
	# commander::printinfo "setup conda fusion-catcher env"
	# tool=fusioncatcher
	# n=${tool//[^[:alpha:]]/}
	# conda create -y -n $n python=2
	# conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
	# for bin in perl samtools bedtools; do
	# 	conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	# done
	# conda activate $n
	# commander::printinfo "downloading databases"
	# download-human-db.sh
	# rm -f $FC_DB_PATH/*.tar.gz* # env variable
	# conda deactivate

	$doclean && {
		commander::printinfo "conda clean up"
		conda clean -y -a
	}

	conda deactivate
	return 0
}

compile::java() {
	local insdir threads url version

	commander::printinfo "installing java"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url="https://download.oracle.com/otn-pub/java/jdk/15.0.1%2B9/51f4f36ad4ef43e39d0dfdbaf6549e32/jdk-15.0.1_linux-x64_bin.tar.gz"
	wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O $insdir/java.tar.gz $url
	version=$(echo $url | perl -lane '$_=~/jdk-([^-_]+)/; print $1')
	tar -xzf $insdir/java.tar.gz -C $insdir
	rm $insdir/java.tar.gz
	mkdir -p $insdir/latest
	ln -sfn $(ls -vd $insdir/jdk-*/bin | tail -1) $insdir/latest/java

	return 0
}

compile::_javawrapper() {
	local java=java
	[[ $3 ]] && java="$3"
	cat <<- EOF > "$1" || return 1
		#!/usr/bin/env bash
		java="$java"
		[[ \$JAVA_HOME && -e "\$JAVA_HOME/bin/java" ]] && java="\$JAVA_HOME/bin/java"
		declare -a jvm_mem_args jvm_prop_args pass_args
		for arg in \$@; do
			case \$arg in
				-D*) jvm_prop_args+=("\$arg");;
				-XX*) jvm_prop_args+=("\$arg");;
				-Xm*) jvm_mem_args+=("\$arg");;
				*) pass_args+=("\$arg");;
			esac
		done
		[[ ! \$jvm_mem_args ]] && jvm_mem_args+=("-Xms1024m") && jvm_mem_args+=("-Xmx4g")
		exec "\$java" "\${jvm_mem_args[@]}" "\${jvm_prop_args[@]}" -jar "$2" "\${pass_args[@]}"
	EOF
	chmod 755 "$1"

	return 0
}

compile::trimmomatic(){
	# conda trimmomatic wrapper is written in python and thus cannot handle process substitutions
	local insdir threads url

	commander::printinfo "installing trimmomatic"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='http://www.usadellab.org/cms/?page=trimmomatic'
	url='http://www.usadellab.org/cms/'$(curl -s $url | grep Version | grep -oE '[^"]+Trimmomatic-[0-9]+\.[0-9]+\.zip' | sort -Vr | head -1)
	wget -q $url -O $insdir/trimmomatic.zip
	unzip -o -d $insdir $insdir/trimmomatic.zip
	rm $insdir/trimmomatic.zip
	cd $(ls -dv $insdir/Trimmomatic-*/ | tail -1)
	mkdir -p $insdir/latest bin
	compile::_javawrapper bin/trimmomatic $(readlink -e trimmomatic-*.jar) $insdir/latest/java/java
	ln -sfn $PWD/bin $insdir/latest/trimmomatic

	return 0
}

compile::sortmerna() {
	local insdir threads url

	commander::printinfo "installing sortmerna"
	echo $insdir
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='https://github.com/biocore/sortmerna/archive/2.1.tar.gz'
	wget -q $url -O $insdir/sortmerna.tar.gz
	tar -xzf $insdir/sortmerna.tar.gz -C $insdir
	rm $insdir/sortmerna.tar.gz
	cd $(ls -dv $insdir/sortmerna-*/ | tail -1)
	make clean || true
	./configure --prefix=$PWD
	make -j $threads
	make install -i # ignore errors caused by --prefix=$PWD
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/sortmerna
	commander::printinfo "indexing databases"
	for i in rRNA_databases/*.fasta; do
		o=index/$(basename $i .fasta)-L18
		echo -ne "bin/indexdb_rna --ref $i,$o -m 4096 -L 18\0"
	done | xargs -0 -P $threads -I {} bash -c {}

	return 0
}

compile::segemehl() {
	local insdir threads url

	commander::printinfo "installing segemehl"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/'
	url="$url"$(curl -s $url | grep -oE 'segemehl-[0-9\.]+\.tar\.gz' | sort -Vr | head -1)
	wget -q $url -O $insdir/segemehl.tar.gz
	tar -xzf $insdir/segemehl.tar.gz -C $insdir
	rm $insdir/segemehl.tar.gz
	cd $(ls -dv $insdir/segemehl-*/ | tail -1)
	export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
	make clean || true
	make -j $threads all
	mkdir -p bin
	mv *.x bin
	touch bin/segemehl bin/haarz
	chmod 755 bin/*
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/segemehl
	cat <<- 'EOF' > $insdir/latest/segemehl/segemehl
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
		l=$(pkg-config --variable=libdir htslib)
		[[ $l ]] && export LD_LIBRARY_PATH=$l
		$(cd $(dirname $0) && echo $PWD)/segemehl.x $*
	EOF
	cat <<- 'EOF' > $insdir/latest/segemehl/haarz
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
		l=$(pkg-config --variable=libdir htslib)
		[[ $l ]] && export LD_LIBRARY_PATH=$l
		$(cd $(dirname $0) && echo $PWD)/haarz.x $*
	EOF

	return 0
}

compile::starfusion() {
	local insdir threads url

	# conda recipe has either star > 2.7.0f (v1.6) or star > 2.5 (>=v1.8, plus python compatibility issues)
	commander::printinfo "installing starfusion"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='https://github.com/'$(curl -s https://github.com/STAR-Fusion/STAR-Fusion/releases | grep -oE 'STAR-Fusion/\S+STAR-Fusion-v[0-9\.]+\.FULL\.tar\.gz' | sort -Vr | head -1)
	wget -q $url -O $insdir/starfusion.tar.gz
	tar -xzf $insdir/starfusion.tar.gz -C $insdir
	rm $insdir/starfusion.tar.gz
	cd $(ls -dv $insdir/STAR-Fusion-*/ | tail -1)
	mkdir -p $insdir/latest
	ln -sfn $PWD $insdir/latest/starfusion

	return 0
}

compile::preparedexseq() {
	local insdir threads

	commander::printinfo "installing dexseq"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	cd $insdir
	rm -rf Subread_to_DEXSeq
	git clone https://github.com/vivekbhr/Subread_to_DEXSeq.git
	cd Subread_to_DEXSeq
	mkdir -p bin
	mv *.py bin
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/subreadtodexseq

	return 0
}

compile::revigo() {
	local insdir threads

	commander::printinfo "installing revigo"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	cd $insdir
	rm -rf revigo
	git clone https://gitlab.leibniz-fli.de/kriege/revigo.git
	cd revigo
	mkdir bin
	compile::_javawrapper bin/revigo $(readlink -e RevigoStandalone.jar) $insdir/latest/java/java
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/revigo

	return 0
}

compile::gem() {
	local insdir threads url version

	commander::printinfo "installing gem"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='https://groups.csail.mit.edu/cgs/gem/download/'
	url="$url"$(curl -s $url | grep -oE "gem.v[0-9\.]+\.tar\.gz" | sort -Vr | head -1)
	version=$(basename $url | sed -E 's/gem.v([0-9\.]+)\.tar\.gz/\1/')
	wget -q $url -O $insdir/gem.tar.gz
	tar -xzf $insdir/gem.tar.gz -C $insdir
	mv $insdir/gem $insdir/gem-$version
	rm $insdir/gem.tar.gz
	cd $insdir/gem-$version
	mkdir -p bin
	wget -q -O bin/Read_Distribution_default.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt
	wget -q -O bin/Read_Distribution_CLIP.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt
	wget -q -O bin/Read_Distribution_ChIP-exo.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_ChIP-exo.txt
	compile::_javawrapper bin/gem $(readlink -e gem.jar) $insdir/latest/java/java
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/gem

	return 0
}

compile::m6aviewer() {
	local insdir threads url version

	commander::printinfo "installing m6aviewer"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='http://dna2.leeds.ac.uk/m6a/m6aViewer_1_6_1.jar'
	version="1.6.1"
	rm -rf $insdir/m6aviewer-$version
	mkdir -p $insdir/m6aviewer-$version/bin
	cd $insdir/m6aviewer-$version
	wget -q $url -O m6aviewer.jar
	compile::_javawrapper bin/m6aviewer $(readlink -e m6aviewer.jar)
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/m6aviewer

	return 0
}

compile::idr() {
	local insdir threads url

	commander::printinfo "installing idr"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='https://github.com/kundajelab/idr'
	url="$url/"$(curl -s $url/tags | grep -oE "archive\/[0-9\.]+\.tar\.gz" | sort -Vr | head -1)
	wget -q $url -O $insdir/idr.tar.gz
	tar -xzf $insdir/idr.tar.gz -C $insdir
	rm $insdir/idr.tar.gz
	cd $(ls -vd $insdir/idr-*/ | tail -1)
	pip install numpy matplotlib
	python setup.py install
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/idr

	return 0
}

compile::newicktopdf(){
	local insdir threads url

	commander::printinfo "installing idr"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	url='ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/njplot/newicktopdf'
	mkdir -p $insdir/newicktopdf
	wget -q $url -O $insdir/newicktopdf/newicktopdf
	chmod 755 $insdir/newicktopdf/newicktopdf
	mkdir -p $insdir/latest
	ln -sfn $insdir/newicktopdf $insdir/latest/newicktopdf
}

compile::ssgsea() {
	local insdir threads

	commander::printinfo "installing ssgsea"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base
	cd $insdir
	rm -rf ssGSEA-gpmodule
	git clone https://github.com/GSEA-MSigDB/ssGSEA-gpmodule.git
	cd ssGSEA-gpmodule
	mv src bin
	chmod 755 bin/*
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/ssgseagpmodule

	git clone https://github.com/broadinstitute/ssGSEA2.0.git
	rm -rf ssGSEA-broad
	mv ssGSEA2.0 ssGSEA-broad
	cd ssGSEA-broad
	find . -type f -exec dos2unix {} \;
	find . -type f -name "*.R" -exec chmod 755 {} \;
	mkdir -p $insdir/latest
	ln -sfn $PWD $insdir/latest/ssgseabroad

	return 0
}

compile::bgztail() {
	local insdir threads url

	commander::printinfo "installing bgztail"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base

	url='https://github.com/'$(curl -s https://github.com/circulosmeos/bgztail/releases | grep -oE 'circulosmeos/bgztail/\S+v[0-9\.]+\.tar\.gz' | sort -Vr | head -1)
	wget -q $url -O $insdir/bgztail.tar.gz
	tar -xzf $insdir/bgztail.tar.gz -C $insdir
	rm $insdir/bgztail.tar.gz
	cd $(ls -dv $insdir/bgztail-*/ | tail -1)
	mkdir -p bin
	mv bgztail bin
	chmod 755 bin/bgztail
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/bgztail

	return 0
}

compile::mdless() {
	local insdir threads url

	commander::printinfo "installing mdless"
	compile::_parse -r insdir -s threads "$@"
	source $insdir/conda/bin/activate base

	url='https://github.com/'$(curl -s https://github.com/ttscoff/mdless/releases | grep -oE 'ttscoff/mdless/\S+\/[0-9\.]+\.tar\.gz' | sort -Vr | head -1)
	wget -q $url -O $insdir/mdless.tar.gz
	tar -xzf $insdir/mdless.tar.gz -C $insdir
	rm $insdir/mdless.tar.gz
	cd $(ls -dv $insdir/mdless-*/bin | tail -1)
	sed -i 's@require@$LOAD_PATH.unshift File.expand_path("../../lib", __FILE__)\nrequire@' mdless
	mkdir -p $insdir/latest
	ln -sfn $PWD $insdir/latest/mdless

	return 0
}