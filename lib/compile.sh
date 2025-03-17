#! /usr/bin/env bash
# (c) Konstantin Riege

function compile::_parse(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
		usage:
			-i <path>      | installation base
			-t <threads>   | number of
			-g <useconfig> | true/false
			-u <upgrade>   | true/false conda envs
		EOF
		return 1
	}

	local OPTIND arg mandatory
	declare -n _insdir_parse _threads_parse _upgrade_parse _useconfig_parse
	while getopts 'r:s:c:f:i:t:u:g:' arg; do
		case $arg in
			# declare references
			r) ((++mandatory)); _insdir_parse="$OPTARG";;
			s) ((++mandatory)); _threads_parse=$OPTARG;;
			c) _upgrade_parse=$OPTARG;;
			f) ((++mandatory)); _useconfig_parse=$OPTARG;;

			# from outside, set variables
			i) ((++mandatory)); _insdir_parse="$OPTARG";;
			t) ((++mandatory)); _threads_parse=$OPTARG;;
			u) _upgrade_parse=$OPTARG;;
			g) ((++mandatory)); _useconfig_parse=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	return 0
}

function compile::all(){
	compile::tools "$@" # needs to be first
	compile::bashbone "$@"
	compile::conda "$@"
	compile::conda_tools "$@"
	compile::java "$@"
	compile::trimmomatic "$@"
	# compile::segemehl "$@" # do not use. see comment in function
	compile::starfusion "$@"
	compile::preparedexseq "$@"
	compile::revigo "$@"
	compile::gem "$@"
	compile::m6aviewer "$@"
	compile::matk "$@"
	compile::newicktopdf "$@"
	compile::ssgsea "$@"
	compile::gztool "$@"
	compile::mdless "$@"
	# compile::moose "$@"

	return 0
}

function compile::lite(){
	compile::tools "$@" # needs to be first
	compile::bashbone "$@"
	compile::gztool "$@"
	compile::mdless "$@"

	return 0
}

function compile::bashbone(){
	local insdir threads cfg version src="$(dirname "$(dirname "$(readlink -e "$0")")")"
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$src/lib/version.sh"
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	rm -f "$insdir/bashbone-$version/scripts/"+(setup|test).sh
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	return 0
}

function compile::tools(){
	local insdir threads cfg i src="$(dirname "$(dirname "$(readlink -e "$0")")")"
	declare -a mapdata

	commander::printinfo "installing pre-compiled tools"
	compile::_parse -r insdir -s threads -f cfg "$@"
	mkdir -p "$insdir/latest"
	mapfile -t mapdata < <(find -L "$src/tools" -mindepth 1 -maxdepth 1 -type f -name "*.tar.gz")
	for i in "${mapdata[@]}"; do
		rm -rf "$insdir/$(basename "$i" .tar.gz)"
		tar -xzf "$i" -C "$insdir"
		ln -sfn "$insdir/$(basename "$i" .tar.gz)/bin" "$insdir/latest/$(basename "$i" | cut -d '-' -f 1)"
	done
	echo -e 'will cite' | "$insdir/latest/parallel/parallel" --citation &> /dev/null || true

	# BASHBONE_TOOLSDIR set to INSDIR/insdir in setup.sh
	_bashbone_setpath

	return 0
}

function compile::upgrade(){
	compile::tools "$@" # needs to be first
	compile::bashbone "$@"
	compile::conda_tools -u true "$@"

	return 0
}

function compile::conda_old(){
	local insdir threads cfg url
	commander::printinfo "installing conda"
	compile::_parse -r insdir -s threads -f cfg "$@"
	url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
	wget -q --show-progress --progress=bar:force -O "$insdir/miniconda.sh" "$url"
	mkdir -p "$insdir/conda"
	bash "$insdir/miniconda.sh" -b -u -f -p "$insdir/conda"
	rm -f "$insdir/miniconda.sh"

	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	conda update -y conda

	# as of 2023 fails with latest conda
	# conda install -y --override-channels -c conda-forge mamba 'python_abi=*=*cp*'
	# conda env config vars set MAMBA_NO_BANNER=1
	# conda config --set changeps1 False
	# alternative use conda with mamba solver
	conda install -y --override-channels -c conda-forge conda-libmamba-solver
	conda config --set solver libmamba

	commander::printinfo "conda clean up"
	# mamba clean -y -a
	conda clean -y -a

	return 0
}

function compile::conda(){
	local insdir threads cfg url
	commander::printinfo "installing conda"
	compile::_parse -r insdir -s threads -f cfg "$@"
	url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
	wget -q --show-progress --progress=bar:force -O "$insdir/miniconda.sh" "$url"
	mkdir -p "$insdir/conda"
	bash "$insdir/miniconda.sh" -b -u -f -p "$insdir/conda"
	rm -f "$insdir/miniconda.sh"

	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	conda env config vars set MAMBA_NO_BANNER=1
	# conda config --set changeps1 False
	# mamba update -y mamba # better dont for sake of miniconda.sh -u and especially, do not update conda!

	commander::printinfo "conda clean up"
	mamba clean -y -a

	return 0
}

function compile::conda_tools(){
	local insdir threads upgrade=false cfg tool n star_version bin doclean=false src="$(dirname "$(dirname "$(readlink -e "$0")")")" f
	declare -A envs
	compile::_parse -r insdir -s threads -c upgrade -f cfg "$@"
	$upgrade && commander::printinfo "validating conda environments" || commander::printinfo "installing conda environments"

	local tmpdir="$insdir/tmp"
	mkdir -p "$tmpdir" "$insdir/config"
	declare -a cmdchk

	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r n; do
		commander::makecmd -a cmdchk -s '|' -o "$tmpdir/$n.yaml" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- 'CMD'
			mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c $([[ $n == "bashbone" ]] && echo r || echo defaults)
		CMD
			grep -vi -e '^prefix:' -e certifi
		CMD
			grep -vE -- '-\s+idr='
		CMD
			sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/'
		CMD
			sed -n '/variables:/,$!p'
		CMD
		envs[$n]=true
	done < <(mamba info -e | awk -v prefix="^$insdir" '$NF ~ prefix {print $1}')
	commander::runcmd -i $threads -a cmdchk

	# setup commonly used tools in bashbone source
	# new: use own bashbone environment to be able to (re)create it from yaml because conda env update -n base --file base.yaml has conflicts
	n=bashbone
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true
		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			# manually remove to be compatible with different version: create had --force in newer version replace by -y
			# as of v1.5.9, non-existing envs cause error
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			# r and main channel is included in conda’s “defaults” channel built by Anaconda Inc.
			# mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults -r -c anaconda \
			# -> too many channels and too many tools slow down and/or break the solver
			# -> may add --strict-channel-priority ,but there is no strict option available for creating env from config file
			# => just keep r channel (bashbone only)
			# yet unresolved bgzip performance drop from 1.17+ https://github.com/samtools/htslib/issues/1767
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c r \
				gcc_linux-64 gxx_linux-64 gfortran_linux-64 \
				glib pkg-config make automake cmake \
				zlib bzip2 pbzip2 \
				git wget curl ghostscript dos2unix \
				sra-tools entrez-direct cgpbigwig ucsc-bigwiginfo wiggletools \
				datamash samtools bedtools ucsc-facount khmer \
				libdeflate "htslib<=1.16" htseq bcftools vcflib vt vcftools \
				perl perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-xml-parser perl-dbi perl-db-file "perl-bioperl>=1.7" perl-bio-eutilities \
				java-jdk \
				nlopt "r-base>=4" \
				r-biocmanager r-devtools r-codetools r-argparser \
				bioconductor-biomart bioconductor-biocparallel bioconductor-genefilter bioconductor-deseq2 bioconductor-dexseq bioconductor-clusterprofiler bioconductor-tcgautils r-r.utils \
				r-survminer bioconductor-impute bioconductor-preprocesscore bioconductor-go.db bioconductor-annotationdbi bioconductor-annotationforge bioconductor-enrichplot bioconductor-rrvgo \
				r-reshape2 r-wgcna r-dplyr r-tidyverse r-ggpubr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-treemap r-data.table r-ggridges r-ashr r-dendextend

			# python stuff
			declare -a cmd1
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
				python -m pip install rapidgzip==0.14.2 indexed-bzip2
			CMD
			# this commit adds sparse gztool-with-lines index from stdin - will be released as 0.14.3 or when lz4 support is implemented as 0.15
			# python3 -m pip install --force-reinstall 'git+https://github.com/mxmlnkn/indexed_bzip2.git@1d2d1a0db5204ec56d7f40f4821d81100737de1f#egginfo=rapidgzip&subdirectory=python/rapidgzip'
		fi

		# perl stuff
		declare -a cmd1
		commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
			cpanm Switch Tree::Simple
		CMD

		# R basics
		echo "rm -rf '$tmpdir'" >> "$BASHBONE_CLEANUP"
		declare -a cmd2
		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(which unzip)");
				Sys.setenv(TAR="$(which tar)");
				install.packages(c("BiocManager","devtools","codetools"),
					repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
			'
		CMD

		# # bioconductor
		# declare -a cmd3
		# commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
		# 	Rscript - <<< '
		# 		options(unzip="$(which unzip)");
		# 		Sys.setenv(TAR="$(which tar)");
		# 		BiocManager::install(c("biomaRt","BiocParallel","genefilter","DESeq2","DEXSeq","clusterProfiler","TCGAutils","TCGAbiolinks","survminer","impute","preprocessCore","GO.db","AnnotationDbi"),
		# 			ask=F, Ncpus=$threads, clean=T, destdir="$tmpdir");
		# 	'
		# CMD

		# # cran - needs to be last since WGCNA depends on bioconductor packages impute,...
		# # R-Forge.r does not complaine about knapsack not being compatible with R>=4
		# declare -a cmd4
		# commander::makecmd -a cmd4 -s '&&' -c {COMMANDER[0]}<<- CMD
		# 	Rscript - <<< '
		# 		options(unzip="$(which unzip)");
		# 		Sys.setenv(TAR="$(which tar)");
		# 		install.packages(c("reshape2","WGCNA","dplyr","tidyverse","ggpubr","ggplot2","gplots","RColorBrewer","svglite","pheatmap","treemap","data.table"),
		# 			repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
		# 		install.packages(c("knapsack"), repos="http://R-Forge.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
		# 	'
		# CMD

		# # github
		# declare -a cmd5
		# commander::makecmd -a cmd5 -s '&&' -c {COMMANDER[0]}<<- CMD
		# 	Rscript - <<< '
		# 		options(unzip="$(which unzip)");
		# 		Sys.setenv(TAR="$(which tar)");
		# 		devtools::install_github("andymckenzie/DGCA", upgrade="never", force=T, clean=T, destdir="$tmpdir");
		# 	'
		# CMD

		# as of 2022, conda can install all r-packages without conflicts. instead manual compilation causes troubles
		# bioconductor-tcgabiolinks conda package is outdated and errornous. Apr 2022 tcga db changed way to access, thus latest tcgabiolinks from git required
		cmd2=()
		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(which unzip)");
				Sys.setenv(TAR="$(which tar)");
				install.packages(c("knapsack"), repos="http://R-Forge.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
				devtools::install_github("andymckenzie/DGCA", upgrade="never", force=T, clean=T, destdir="$tmpdir");
				devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI.data", upgrade="never", force=T, clean=T, destdir="$tmpdir");
				devtools::install_github("BioinformaticsFMRP/TCGAbiolinks", upgrade="never", force=T, clean=T, destdir="$tmpdir");
			'
		CMD

		commander::runcmd -c bashbone -i 1 -a cmd1
		commander::runcmd -c bashbone -i 1 -a cmd2

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c r | grep -vi -e "^prefix:" -e certifi | grep -vE -- '-\s+idr=' | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"
	}

	# better do not predefine python version. if tool recipe depends on earlier version, conda installs an older or the oldest version (freebayes)
	for tool in fastqc cutadapt rcorrector star bwa bbmap rseqc subread htseq picard bamutil fgbio macs2 genrich peakachu diego gatk4 freebayes varscan igv intervene deeptools raxml metilene dupsifter umi_tools methyldackel idr clust seacr gopeaks homer kmc; do
		n=${tool/=*/}
		n=${n//[^[:alpha:]]/}
		[[ $tool == "bwa" ]] && tool+=" bwa-mem2"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
		else
			$upgrade && ${envs[$n]:=false}
		fi || {
			doclean=true

			commander::printinfo "setup conda $n env"
			if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
				mamba env remove -y -n $n || true
				mamba env create -n $n --file "$src/config/$n.yaml"
			else
				mamba create -y -n $n
				mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool
			fi

			mkdir -p "$insdir/config"
			mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

			for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
				mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
			done
		}
	done
	# necessary fixes
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl"
	sed -i 's/gzip -cd/rapidgzip -kcd -P 4/' "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl"

	cat <<-EOF > "$insdir/conda/envs/seacr/bin/SEACR.sh"
		#!/usr/bin/env bash
		exec "\$(realpath -s "\$(dirname "\$0")")/$(basename "$insdir/conda/envs/seacr/bin/"SEACR_*.sh)" "\$@"
	EOF
	chmod 755 "$insdir/conda/envs/seacr/bin/SEACR.sh"

	# tool=gopeaks
	# n=${tool/=*/}
	# n=${n//[^[:alpha:]]/}
	# if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
	# 	$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	# else
	# 	$upgrade && ${envs[$n]:=false}
	# fi || {
	# 	doclean=true

	# 	commander::printinfo "setup conda $n env"
	# 	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
	#		mamba env remove -y -n $n || true
	# 		mamba env create -n $n --file "$src/config/$n.yaml"
	# 	else
	# 		mamba create -y -n $n
	# 		mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults gopeaks pandas numpy pybedtools scipy pysam
	# 	fi

	# 	git clone https://github.com/gartician/summit.git "$insdir/conda/envs/gopeaks/src"
	# 	echo '#!/usr/bin/env python' > "$insdir/conda/envs/gopeaks/bin/summit.py"
	# 	cat "$insdir/conda/envs/gopeaks/src/summit/summit.py" >> "$insdir/conda/envs/gopeaks/bin/summit.py"
	# 	chmod 755 "$insdir/conda/envs/gopeaks/bin/summit.py"

	# 	mkdir -p "$insdir/config"
	# 	mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

	# 	for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
	# 		mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	# 	done
	# }

	star_version=$(mamba list -n star -f star | tail -1 | awk '{print $2}')
	# manual setup of requirements from bioconda meta.yaml (see compile::starfusion) due to non-latest installation via conda
	# note: recent star (star indexer 2.7.1a) is not compatible with CTAT plug-n-play genome index for star-fusion v1.9 (star indexer 2.4.1)
	tool=star-fusion
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n #python=3
			# propably enought: perl perl-set-intervaltree perl-carp perl-carp-assert perl-db-file perl-io-gzip perl-json-xs perl-uri \
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults \
				perl perl-file-path perl-getopt-long perl-set-intervaltree perl-carp perl-carp-assert perl-data-dumper perl-findbin perl-db-file perl-io-gzip perl-json-xs perl-uri perl-list-moreutils perl-list-util perl-storable \
				igv-reports "star=$star_version" gmap bowtie bbmap samtools blast
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	# tool="sortmerna<3"
	# n=${tool/=*/}
	# n=${n//[^[:alpha:]]/}
	# if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
	# 	$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	# else
	# 	$upgrade && ${envs[$n]:=false}
	# fi || {
	# 	doclean=true

	# 	commander::printinfo "setup conda $n env"
	# 	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
	# 		mamba env remove -y -n $n || true
	# 		mamba env create -n $n --file "$src/config/$n.yaml"
	# 	else
	# 		mamba create -y -n $n
	# 		mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults "$tool"
	# 	fi

	# 	git clone --depth 1 https://github.com/biocore/sortmerna "$insdir/conda/envs/sortmerna/src"
	# 	mv "$insdir/conda/envs/sortmerna/src/data/rRNA_databases/" "$insdir/conda/envs/sortmerna/"
	# 	mkdir -p "$insdir/conda/envs/sortmerna/rRNA_databases/index"
	# 	rm -rf "$insdir/conda/envs/sortmerna/src"

	# 	declare -a cmdidx
	# 	for f in "$insdir/conda/envs/sortmerna/rRNA_databases/"*.fasta; do
	# 		commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD
	# 			indexdb_rna
	# 			--ref "$f","$insdir/conda/envs/sortmerna/rRNA_databases/index/$(basename "$f" .fasta)-L18"
	# 			-m 4096
	# 			-L 18
	# 		CMD
	# 	done
	# 	commander::runcmd -c sortmerna -i $threads -a cmdidx

	# 	mkdir -p "$insdir/config"
	# 	mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

	# 	for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
	# 		mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	# 	done
	# }

	tool="sortmerna"
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	# conda version of 4.3.7 core dumps! compile manually!
	if [[ -e "$src/config/$n.yaml-noexist" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml-noexist" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults gcc_linux-64=11.4 gxx_linux-64=11.4 glib pkg-config make automake cmake conda-build
			#mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults sortmerna
		fi

		mkdir -p "$insdir/conda/envs/$n/rRNA_databases/index" "$insdir/conda/envs/$n/src"
		wget -q --show-progress --progress=bar:force "https://github.com/biocore/$n/releases/download/v4.3.4/database.tar.gz" -O "$insdir/conda/envs/$n/rRNA_databases/database.tar.gz"
		tar -xzf "$insdir/conda/envs/$n/rRNA_databases/database.tar.gz" -C "$insdir/conda/envs/$n/rRNA_databases"
		chmod 666 "$insdir/conda/envs/$n/rRNA_databases/"*.fasta
		rm -f "$insdir/conda/envs/$n/rRNA_databases/database.tar.gz"

		declare -a cmdidx
		local tmp="$(mktemp -d -p "${TMPDIR:-/tmp}" cleanup.XXXXXXXXXX.sortmerna)"
		rm -f "$insdir/conda/envs/$n/conda-bld/linux-64/$n"*

		commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD {COMMANDER[7]}<<- CMD
			cd "$insdir/conda/envs/$n/src"
		CMD
			git init
		CMD
			git remote add origin https://github.com/conda-forge/staged-recipes
		CMD
			git fetch origin pull/25438/head
		CMD
			git reset --hard FETCH_HEAD
		CMD
			conda build --no-anaconda-upload recipes/sortmerna
		CMD
			mamba install -y --override-channels -c conda-forge -c bioconda -c defaults "$insdir/conda/envs/$n/conda-bld/linux-64/sortmerna-4.3.7"*.tar.bz2
		CMD
			sortmerna
				--ref "$insdir/conda/envs/$n/rRNA_databases/smr_v4.3_fast_db.fasta"
				--index 1
				-L 18
				-m 4096
				--threads $threads
				--workdir "$tmp/"
				--idx-dir "$insdir/conda/envs/$n/rRNA_databases/index"
		CMD
		# path can be sortmerna-4.3.7.beta.1-py312hde9dbc4_0.tar.bz2 or sortmerna-4.3.7-py312hde9dbc4_0.tar.bz2

		# commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD
		#	sortmerna
		#		--ref "$insdir/conda/envs/$n/rRNA_databases/smr_v4.3_fast_db.fasta"
		#		--index 1
		#		-L 18
		#		-m 4096
		#		--threads $threads
		#		--workdir "$tmp/"
		#		--idx-dir "$insdir/conda/envs/$n/rRNA_databases/index"
		# CMD
		commander::runcmd -c sortmerna -i $threads -a cmdidx

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	# arriba 2.x , successor of 1.2 (arriba=1.2) has new star parameters incompatible with star < 2.7.6
	tool=arriba
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults "$tool$(awk -F '.' '{if ($1>=2 && $2>=7){print ">=2"}else{print "<2"}}' <<< $star_version)'" "star=$star_version"
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=bwameth
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool bwa-mem2
		fi
		# get latest functions from pull requensts like support for bwa-mem2 and report of supplementary/split alignments
		# use absolute path, because on some machines sometimes bwa-mem2 aborts with error: prefix is too long
		curl -s "https://raw.githubusercontent.com/brentp/bwa-meth/master/bwameth.py" | \
			sed -E -e '/--threads/{s/$/\n    p.add_argument("-s", "--score", type=int, default=40)/}' \
			-e '/threads=args.threads/{s/$/\n            score=args.score,/}' \
			-e 's/threads=1,/threads=1, score=40,/' \
			-e 's@"bwa-mem2 index@"\\\"" + os.path.dirname(__file__) + "/bwa-mem2\\\" index@' \
			-e 's@bwa-mem2 mem -T 40 -B 2 -L 10 -CM@ \\\"" + os.path.dirname(__file__) + f"/bwa-mem2\\\" mem -T {score} -a -B 2 -L 10 -C -Y@' \
		> "$insdir/conda/envs/bwameth/bin/bwameth.py"
		# by removal of -M splits/chimeric reads are marked as supplementary (which is the way to go!).
		# -Y: apply soft-clipping instead of hard clipping to keep sequence info in bam (can be changed via)
		# squeeze in score parameter to control bwa minoutscore

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=vardict
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool vardict-java readline=6
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=snpeff
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n #python=3
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool snpsift
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=platypus-variant
	n=platypus
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=danpos
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"

		rm -rf "$insdir/DANPOS3"
		mkdir -p "$insdir/DANPOS3"
		cd "$insdir/DANPOS3"
		git init
		git remote add origin https://github.com/sklasfeld/DANPOS3
		git fetch origin 665bdd115f77a0e4c4f2f07f53f42b34490d703c
		git reset --hard FETCH_HEAD
		# chmod 755 *.py
		cd - > /dev/null

		sed -i -E 's/,\s*lower_tail\s*=\s*False,\s*log_bool\s*=\s*True\s*//' "$insdir/DANPOS3/reads.py"

		for x in $(grep -nF 'while( (0-functions.div(float(functions.ppois(height,m.item()).split()[-1]),log(10))) < pheight):height+=1' "$insdir/DANPOS3/wig.py" | cut -d : -f 1); do
			sed -i -E "$x,$x{s/^(\s*).+/\1while( (0-functions.div(float(str(functions.ppois(height, m.item())).split('.')[-1]), log(10))) < pheight):height+=1/}" "$insdir/DANPOS3/wig.py"
		done
		for x in $(grep -nF 'pvl=functions.div(float(functions.ppois(v.item(),m.item()).split()[-1]),log(10))' "$insdir/DANPOS3/wig.py" | cut -d : -f 1); do
			sed -i -E "$x,$x{s/^(\s*).+/\1pvl=functions.div(float(str(functions.ppois(v.item(),m.item())).split('.')[-1]),log(10))/}" "$insdir/DANPOS3/wig.py"
		done
		for x in $(grep -nF 'while( (0-functions.div(float(functions.ppois(height,m.item())).split()[-1],log(10))) < pheight):height+=1' "$insdir/DANPOS3/wig.py" | cut -d : -f 1); do
			sed -i -E "$x,$x{s/^(\s*).+/\1while( (0-functions.div(float(str(functions.ppois(height,m.item())).split('.')[-1]),log(10))) < pheight):height+=1/}" "$insdir/DANPOS3/wig.py"
		done
		for x in $(grep -nF 'while( (0-float(functions.ppois(height,m.item())/log(10).split()[-1])) < pcut):height+=1' "$insdir/DANPOS3/wig.py" | cut -d : -f 1); do
			sed -i -E "$x,$x{s@^(\s*).+@\1while( (0-float(functions.ppois(height,m.item())/float(str(log(10)).split()[-1]))) < pcut):height+=1@}" "$insdir/DANPOS3/wig.py"
		done

		mkdir -p "$insdir/latest"
		ln -sfn "$insdir/DANPOS3" "$insdir/latest/danpos"

		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			# ucsc from v400 does not allow stdin anymore
			mamba install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults wiggletools samtools scipy "r-base>=4" $(awk '!/^\s*#/{print $1}' "$insdir/latest/danpos/requirements.txt")
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=salmon
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
		$upgrade && diff "$src/config/$n.yaml" "$tmpdir/$n.yaml" &> /dev/null
	else
		$upgrade && ${envs[$n]:=false}
	fi || {
		doclean=true

		commander::printinfo "setup conda $n env"
		if [[ -e "$src/config/$n.yaml" ]] && $cfg; then
			mamba env remove -y -n $n || true
			mamba env create -n $n --file "$src/config/$n.yaml"
		else
			mamba create -y -n $n
			mamba install -n $n --override-channels -c conda-forge -c bioconda -c defaults $n snakemake docopt pandas r-tidyverse r-scales r-WriteXLS r-BiocManager
		fi

		mkdir -p "$insdir/config"
		mamba env export -n $n --no-builds --override-channels -c conda-forge -c bioconda -c defaults | grep -vi -e "^prefix:" -e certifi | sed -E 's/=+/==/; s/- r-base==.*/- r-base>=4/' | sed -n '/variables:/,$!p' > "$insdir/config/$n.yaml"

		for bin in perl bgzip rapidgzip ibzip2 samtools bcftools bedtools vcfsamplediff bg2bw bwcat; do
			mamba list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/envs/bashbone/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done

		rm -rf "$insdir/SalmonTE"
		mkdir -p "$insdir/SalmonTE"
		cd "$insdir/SalmonTE"
		git init
		git remote add origin https://github.com/hyunhwan-jeong/SalmonTE
		git fetch origin 0a405190b7b5796b646814f14870b268bdf33249
		git reset --hard FETCH_HEAD
		cd - > /dev/null
		sed -i 's/quant --minAssignedFrags/quant --discardOrphansQuasi --minAssignedFrags/' "$insdir/SalmonTE/snakemake/Snakefile.paired"
		mv "$insdir/SalmonTE/salmon/linux/bin/salmon" "$insdir/SalmonTE/salmon/linux/bin/salmon.old"
		ln -sfnr "$insdir/conda/envs/$n/bin/salmon" "$insdir/SalmonTE/salmon/linux/bin/salmon"
		ln -sfn "$insdir/SalmonTE" "$insdir/latest/salmon"
	}

	# this is a pipeline itself with own genome and databases and thus will not be part of bashbone
	# commander::printinfo "setup conda fusion-catcher env"
	# tool=fusioncatcher
	# n=${tool//[^[:alpha:]]/}
	# conda create -y -n $n
	# conda install -n $n -y --override-channels -c conda-forge -c bioconda -c defaults $tool
	# for bin in perl samtools bedtools; do
	# 	conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	# done
	# conda activate $n
	# commander::printinfo "downloading databases"
	# download-human-db.sh
	# rm -f $FC_DB_PATH/*.tar.gz* # env variable

	$doclean && {
		commander::printinfo "conda clean up"
		mamba clean -y -a
	}

	return 0
}

function compile::java(){
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/java.url"
	commander::printinfo "installing java"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='https://download.oracle.com/otn-pub/java/jdk/15.0.1%2B9/51f4f36ad4ef43e39d0dfdbaf6549e32/jdk-15.0.1_linux-x64_bin.tar.gz'
		url='https://download.oracle.com/java/17/archive/jdk-17.0.12_linux-x64_bin.tar.gz'
		url=$(wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O - 'https://www.oracle.com/java/technologies/downloads' | grep -oE 'https[^"]+latest/jdk-[0-9.]+_linux-x64_bin.tar.gz' | sort -Vr | head -1)
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/java.url"
	wget -q --show-progress --progress=bar:force --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" "$url" -O "$insdir/java.tar.gz"
	tar -xzf "$insdir/java.tar.gz" -C "$insdir"
	rm "$insdir/java.tar.gz"
	mkdir -p "$insdir/latest"
	ln -sfn "$(ls -vd "$insdir/jdk-"*/bin | tail -1)" $insdir/latest/java

	return 0
}

function compile::_javawrapper(){
	local java=java
	[[ $3 ]] && java="$3"
	cat <<- EOF > "$1" || return 1
		#!/usr/bin/env bash
		java="$java"
		[[ \$JAVA_HOME && -e "\$JAVA_HOME/bin/java" ]] && java="\$JAVA_HOME/bin/java"
		declare -a jvm_mem_args jvm_prop_args pass_args
		for arg in "\$@"; do
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

function compile::trimmomatic(){
	# conda trimmomatic wrapper is written in python and thus cannot handle process substitutions
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/trimmomatic.url"
	commander::printinfo "installing trimmomatic"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='http://www.usadellab.org/cms/?page=trimmomatic'
		url='http://www.usadellab.org/cms/'$(curl -s "$url" | grep Version | grep -oE '[^"]+Trimmomatic-[0-9]+\.[0-9]+\.zip' | head -1)
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/trimmomatic.url"

	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/trimmomatic.zip"
	unzip -o -d "$insdir" "$insdir/trimmomatic.zip"
	rm "$insdir/trimmomatic.zip"
	cd "$(ls -dv "$insdir/Trimmomatic-"*/ | tail -1)"
	mkdir -p "$insdir/latest" bin
	compile::_javawrapper bin/trimmomatic "$(readlink -e trimmomatic-*.jar)" "$insdir/latest/java/java"
	ln -sfn "$PWD/bin" "$insdir/latest/trimmomatic"
	cd - > /dev/null

	return 0
}

function compile::segemehl(){
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/segemehl.url"
	commander::printinfo "installing segemehl"
	compile::_parse -r insdir -s threads -f cfg "$@"

	# attention: official source does not allow multi-member gzip/bgzip files and has a haarz bug
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/'
		url='http://legacy.bioinf.uni-leipzig.de/Software/segemehl/'
		url="$url"$(curl -s "$url" | grep -oE 'downloads/segemehl-[0-9\.]+\.tar\.gz' | head -1)
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/segemehl.url"

	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/segemehl.tar.gz"
	tar -xzf "$insdir/segemehl.tar.gz" -C "$insdir"
	rm "$insdir/segemehl.tar.gz"
	cd "$(ls -dv "$insdir/segemehl-"*/ | tail -1)"
	export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
	make clean || true
	make -j $threads all

	mkdir -p bin
	mv *.x bin
	touch bin/segemehl bin/haarz
	chmod 755 bin/*
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/segemehl"
	cd - > /dev/null

	cat <<- 'EOF' > "$insdir/latest/segemehl/segemehl"
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
		l="$(pkg-config --variable=libdir htslib)"
		[[ $l ]] && export LD_LIBRARY_PATH="$l"
		unset MALLOC_ARENA_MAX
		exec "$(realpath -se "$(dirname "$0")")/segemehl.x" "$@"
	EOF
	cat <<- 'EOF' > "$insdir/latest/segemehl/haarz"
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
		l="$(pkg-config --variable=libdir htslib)"
		[[ $l ]] && export LD_LIBRARY_PATH="$l"
		unset MALLOC_ARENA_MAX
		exec "$(realpath -se "$(dirname "$0")")/haarz.x" "$@"
	EOF

	return 0
}

function compile::starfusion(){
	# conda recipe has either star > 2.7.0f (v1.6) or star > 2.5 (>=v1.8, plus python compatibility issues)
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/starfusion.url"
	commander::printinfo "installing starfusion"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url=$(curl -s https://api.github.com/repos/STAR-Fusion/STAR-Fusion/releases | grep -E 'browser_download_url.*\.FULL\.tar\.gz"$' | head -1 | sed -E 's/.*"([^"]+)"$/\1/')
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/starfusion.url"

	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/starfusion.tar.gz"
	tar -xzf "$insdir/starfusion.tar.gz" -C "$insdir"
	rm "$insdir/starfusion.tar.gz"
	cd "$(ls -dv "$insdir/STAR-Fusion"-*/ | tail -1)"
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD" "$insdir/latest/starfusion"
	cd - > /dev/null

	return 0
}

function compile::preparedexseq(){
	local insdir threads cfg
	commander::printinfo "installing dexseq"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	rm -rf "$insdir/Subread_to_DEXSeq"
	mkdir -p "$insdir/Subread_to_DEXSeq"
	cd "$insdir/Subread_to_DEXSeq"
	git init
	git remote add origin https://github.com/vivekbhr/Subread_to_DEXSeq
	git fetch origin 23ccd2ff4ff605622d3f3575130b31b6761a0990
	git reset --hard FETCH_HEAD
	mkdir -p bin
	mv *.py bin
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/subreadtodexseq"
	cd - > /dev/null

	return 0
}

function compile::revigo(){
	local insdir threads cfg
	commander::printinfo "installing revigo"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	rm -rf "$insdir/revigo"
	git clone --depth 1 https://gitlab.leibniz-fli.de/kriege/revigo "$insdir/revigo"
	cd "$insdir/revigo"
	mkdir -p bin
	compile::_javawrapper bin/revigo "$(readlink -e RevigoStandalone.jar)" "$insdir/latest/java/java"
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/revigo"
	cd - > /dev/null

	return 0
}

function compile::gem(){
	local insdir threads cfg url version src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/gem.url"
	commander::printinfo "installing gem"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='https://groups.csail.mit.edu/cgs/gem/download/'
		url="$url"$(curl -s "$url" | grep -oE "gem.v[0-9\.]+\.tar\.gz" | tail -1)
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/gem.url"

	version=$(basename "$url" | sed -E 's/gem.v([0-9\.]+)\.tar\.gz/\1/')
	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/gem.tar.gz"
	tar -xzf "$insdir/gem.tar.gz" -C "$insdir"
	mv "$insdir/gem" "$insdir/gem-$version"
	rm "$insdir/gem.tar."gz
	cd "$insdir/gem-$version"
	mkdir -p bin
	wget -q --show-progress --progress=bar:force -O bin/Read_Distribution_default.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt
	wget -q --show-progress --progress=bar:force -O bin/Read_Distribution_CLIP.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt
	wget -q --show-progress --progress=bar:force -O bin/Read_Distribution_ChIP-exo.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_ChIP-exo.txt
	compile::_javawrapper bin/gem "$(readlink -e gem.jar)" "$insdir/latest/java/java"
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/gem"
	cd - > /dev/null

	return 0
}

function compile::m6aviewer(){
	local insdir threads cfg url version fallback=false src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/m6aviewer.url"
	commander::printinfo "installing m6aviewer"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		# url='http://dna2.leeds.ac.uk/m6a/m6aViewer_1_6_1.jar'
		url='https://gitlab.leibniz-fli.de/kriege/m6aviewer/-/raw/main/m6aviewer.jar'
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/m6aviewer.url"

	version="1.6.1"
	mkdir -p "$insdir/m6aviewer-$version"
	wget -q --show-progress --progress=bar:force --timeout=2 --waitretry=2 --tries=2 "$url" -O "$insdir/m6aviewer-$version/m6aviewer.jar" || fallback=true
	if $fallback; then
		rm -rf "$insdir/m6aviewer-$version"
		git clone --depth 1 https://gitlab.leibniz-fli.de/kriege/m6aviewer "$insdir/m6aviewer-$version"
	fi
	cd "$insdir/m6aviewer-$version"
	mkdir -p bin
	compile::_javawrapper bin/m6aviewer "$(readlink -e m6aviewer.jar)"
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/m6aviewer"
	cd - > /dev/null

	return 0
}

function compile::matk(){
	local insdir threads cfg url version src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/matk.url"
	commander::printinfo "installing matk" # deeprip
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='http://matk.renlab.org/download/MATK-1.0.jar'
		url='https://github.com/lazyky/MATK_backup/releases/download/v0.1dev/MATK-1.0.jar'
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/matk.url"

	version="1.0"
	mkdir -p "$insdir/matk-$version"
	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/matk-$version/matk.jar"
	cd "$insdir/matk-$version"
	mkdir -p bin
	compile::_javawrapper bin/matk "$(readlink -e matk.jar)"
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/m6aviewer"
	cd - > /dev/null

	return 0
}

function compile::newicktopdf(){
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/newicktopdf.url"
	commander::printinfo "installing newicktopdf"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='ftp://pbil.univ-lyon1.fr/pub/mol_phylogeny/njplot/newicktopdf'
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/newicktopdf.url"

	mkdir -p "$insdir/newicktopdf"
	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/newicktopdf/newicktopdf"
	chmod 755 "$insdir/newicktopdf/newicktopdf"
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/newicktopdf" "$insdir/latest/newicktopdf"
}

function compile::ssgsea(){
	local insdir threads cfg url dir src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/ssgsea.url"
	commander::printinfo "installing ssgsea"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url=$(curl -s https://api.github.com/repos/GSEA-MSigDB/ssGSEA-gpmodule/releases | grep -E 'browser_download_url.*\.zip' | head -1 | sed -E 's/.*"([^"]+)"$/\1/')
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/ssgsea.url"

	dir="$insdir/$(basename "$url" .zip)"
	mkdir -p "$dir/bin"
	cd "$dir"
	wget -q --show-progress --progress=bar:force "$url" -O ssgsea.zip
	unzip -o ssgsea.zip
	rm ssgsea.zip
	chmod 755 *.R
	mv *.R bin
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/ssgsea"
	cd - > /dev/null

	# git clone https://github.com/broadinstitute/ssGSEA2.0.git
	# rm -rf ssGSEA-broad
	# mv ssGSEA2.0 ssGSEA-broad
	# cd ssGSEA-broad
	# find . -type f -exec dos2unix "{}" \;
	# find . -type f -name "*.R" -exec chmod 755 "{}" \;
	# mkdir -p "$insdir/latest"
	# ln -sfn "$PWD" "$insdir/latest/ssgseabroad"

	return 0
}

function compile::gztool(){
	local insdir threads cfg url version src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/gztool.url"
	commander::printinfo "installing gztool"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url=$(curl -s https://api.github.com/repos/circulosmeos/gztool/releases | grep -E 'browser_download_url.*linux\.x86_64"$' | head -1 | sed -E 's/.*"([^"]+)"$/\1/')
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/gztool.url"

	version="$(basename "$(dirname "$url")")"
	version=${version/v/}
	mkdir -p "$insdir/gztool-$version/bin"
	cd "$insdir/gztool-$version"
	wget -q --show-progress --progress=bar:force "$url" -O bin/gztool
	chmod 755 bin/gztool
	cat <<-'EOF' > bin/gztail
	#! /usr/bin/env bash
	[[ $# -lt 1 || "$1" =~ (^-h|[[:space:]]-h) ]] && {
	    echo "usage: $(basename "$0") [<#>] <file.gz>"
	    exit
	}
	if [[ $2 ]]; then
	    n=$1
	    n=$((--n))
	    f="$2"
	else
	    n=9
	    f="$1"
	fi
	p="$(realpath -se "$(dirname "$0")")"
	l=$("$p/gztool" -l "$f" |& sed -nE 's/.*\s+lines\s+:\s+([0-9]+).*/\1/p')
	"$p/gztool" -v 0 -L $((l-n)) "$f"
	EOF
	chmod 755 bin/gztail
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD/bin" "$insdir/latest/gztool"
	cd - > /dev/null

	return 0
}

function compile::mdless(){
	local insdir threads cfg url src="$(dirname "$(dirname "$(readlink -e "$0")")")/config/mdless.url"
	commander::printinfo "installing mdless"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	if [[ -e "$src" ]] && $cfg; then
		url=$(cat "$src")
	else
		url='https://github.com/ttscoff/mdless'
		url="$url/"$(curl -s "$url/tags" | grep -oE "archive\S+\/[0-9\.]+\.tar\.gz" | head -1)
	fi
	mkdir -p "$insdir/config"
	echo "$url" > "$insdir/config/mdless.url"

	wget -q --show-progress --progress=bar:force "$url" -O "$insdir/mdless.tar.gz"
	tar -xzf "$insdir/mdless.tar.gz" -C "$insdir"
	rm "$insdir/mdless.tar.gz"
	cd "$(ls -dv "$insdir/mdless"-*/bin | tail -1)"
	sed -i 's@require@$LOAD_PATH.unshift(File.expand_path("../../lib", __FILE__))\nrequire@' mdless
	mkdir -p "$insdir/latest"
	ln -sfn "$PWD" "$insdir/latest/mdless"
	cd - > /dev/null

	return 0
}

function compile::moose(){
	local insdir threads cfg
	commander::printinfo "installing moose2"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$insdir/conda/bin/activate" bashbone
	rm -rf "$insdir/moose2"
	mkdir -p "$insdir/moose2"
	cd "$insdir/moose2"
	git init
	git remote add origin https://github.com/grabherr/moose2
	git fetch origin a14ba3ef463c479c8037d1bb986f3cb6b5525da5
	git reset --hard FETCH_HEAD
	sed -iE 's/\s*-static//' Makefile
	export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
	make clean || true
	make -j $threads
	mkdir -p bin
	find . -maxdepth 1 -type f -executable -exec mv {} bin \;
	ln -sfn "$PWD/bin" "$insdir/latest/moose"
	cd - > /dev/null

	return 0
}
