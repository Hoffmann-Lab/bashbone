# Bashbone

A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

## For developers

- Write code in your favorite language as Here-documents for later orchestrated execution
- Add object-oriented programming (oop) like syntactic sugar to bash variables and arrays to avoid complex parameter-expansions, variable-expansions and brace-expansions
- Execute commands in parallel on your machine or submit them as jobs to a sun grid engine (SGE)
- Infer number of parallel instances according to targeted memory consumption or targeted threads per instance
- Get a full bash error stack trace upon SIGERR

## For users

- Easily design multi-threaded pipelines to perform NGS related tasks
- For model and non-model organisms
- Availability of many best-practice parameterized and run-time tweaked software wrappers
- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks

- For paired-end and single-end derived raw sequencing or prior mapped read data
  - RNA-Seq protocols (RNA, RIP, m6A, ..)
  - DNA-Seq protocols (WGS, ChIP, Chip-exo, ATAC, CAGE, Quant, Cut&Tag, ..)
  - Bisulfite converted DNA-Seq protocols (WGBS, RRBS)
- Data quality anlysis and preprocessing
  - adapter and poly-mono/di-nucleotide clipping
  - quality trimming
  - error correction
  - artificial rRNA depletion
- Read alignment and post-processing
  - knapsack problem based slicing of alignment files for parallel task execution
  - sorting, filtering,
  - UMI based deduplication or removal of optical and PCR duplicates
  - generation of pools and pseudo-replicates
  - read group modification, split N-cigar reads, left-alignment and base quality score recalibration
- Gene fusion detection
- Methyl-C calling and prediction of differentially methylated regions
- Expression analysis
  - Read quantification, TPM and Z-score normalization and heatmap plotting
  - Inference of strand specific library preparation methods
  - Inference of differential expression as well as clusters of co-expression
  - Detection of differential splice junctions and differential exon usage
  - Gene ontology (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering
- Implementation of Encode3 best-practice ChIP-Seq Peak calling
  - Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related IP-Seq data
  - Inference of effective genome sizes
- Variant detection from DNA or RNA sequencing experiments
  - Integration of multiple solutions for germline and somatic calling
  - VCF normalization
  - Tree reconstruction from homozygous sites
- SSGSEA and survival analysis from TCGA cancer expression data
- Genome and SRA data retrieval
  - Genome to transcriptome conversion
  - Data visualization via IGV batch processing

# License

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download

This will download you a copy which includes the latest developments

```bash
git clone --recursive https://github.com/Hoffmann-Lab/bashbone
```

To check out the latest release (irregularly compiled) do

```bash
cd bashbone
git checkout $(git describe --tags)
```

# Quick start (without installation)

Please see installation section to get all third-party tools set-upped and subsequently all functions to work properly.

Load the library and list available functions. Each function comes with a usage.

```bash
source ./activate.sh
bashbone -h
bashbone -l
```

## Developers centerpiece

To create commands, print them and to execute them in parallel (according to your resources), utilize

```bash
commander::makecmd
commander::printcmd
commander::runcmd
# or
commander::qsubcmd
```

### Example

Showcase of the builtin file descriptors array `COMMANDER` with separator (`-s '|'`) based concatenation of interpreted and non-interpreted Here-Documents. Each concatenation will be stored as a single command in a de-referenced array (`-a cmds`) Optionally, an output file can be defined (`-o <path/to/file>`). `BASHBONE_ERROR` can be used to alter the default error message.

```bash
declare -a cmds
threads=4
for i in sun mercury venus earth mars jupiter saturn uranus neptune pluto; do
	commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		sleep $((RANDOM%10+1));
		echo $i
	CMD
		awk '{
			print "hello "$1
		}'
	CMD
done
commander::printcmd -a cmds
BASHBONE_ERROR="hello world failed"
commander::runcmd -v -b -t $threads -a cmds
unset BASHBONE_ERROR
```

### OOP bash

A small excerpt of possible member functions. See `bashbone -a` helper underscore functions for a full list.

```bash
declare -a x
helper::addmemberfunctions -v x
x.push "hello"
x.push "world"
x.print
x.substring 2 4
x.print
x.get -1
x.pop
x.print
x.shift
x.print
<x.*>
```

# Installation

## Full installation of all third party tools used in bashbone functions

```bash
./setup.sh -i all -d <path/to/installation>
source <path/of/installation/latest/bashbone/activate.sh>
bashbone -h
```

## Upgrade to a newer release (sources only)

```bash
./setup.sh -i upgrade -d <path/of/installation>
```

## Update tools

The setup routine will always install the latest software via conda, which can be updated by running the related setup functions again.

```bash
./setup.sh -i conda_tools -d <path/of/installation>
```

Trimmomatic, segemehl, STAR-Fusion, GEM and mdless will be installed next to the conda environments. If new releases are available, they will be automatically fetched and installed upon running the related setup functions again.

```bash
./setup.sh -i trimmomatic,segemehl,starfusion,gem,mdless -d <path/of/installation>
```

# Usage

To load bashbone execute
```bash
source <path/of/installation/latest/bashbone/activate.sh>
bashbone -h
```

In order to get all function work properly, enable bashbone to use conda environments. Conda can be disabled analogously.
```bash
bashbone -c
bashbone -s
```

Shortcut:

```bash
source <path/of/installation/latest/bashbone/activate.sh> -c true
bashbone -s
```

## Retrieve SRA datasets

Use the enclosed script to fetch sequencing data from SRA

```bash
source <path/of/installation/latest/bashbone/activate.sh> -c true
sra-dump.sh -h
```

## Retrieve genomes

Use the enclosed script to fetch human hg19/hg38 or mouse mm9/mm10 genomes and annotations. Plug-n-play CTAT genome resource made for gene fusion detection and shipped with STAR index can be selected optionally.

```bash
source <path/of/installation/latest/bashbone/activate.sh> -c true
dlgenome.sh -h
```

## Sample info file

In order to perform desired comparative tasks, some functions require a sample info file.

Assume this input:
<br>

| Treatment   | Replicate 1       | Replicate 2       |
| :---        | :---:             | :---:             |
| wild-type   | path/to/wt1.fq    | path/to/wt2.fq    |
| treatment A | path/to/trA_1.fq  | path/to/trA_2.fq  |
| treatment B | path/to/trB.n1.fq | path/to/trB.n2.fq |

And this desired output (N=2 vs N=2 each):

- wt_vs_A
- wt_vs_B
- A_vs_B

Then the info file should consist of:

- At least 4 tab-separated columns (`<name>`, `<main-factor>`, `single-end|paired-end`, `<replicate>`)
- Optionally, additional factors
- First column needs to consist of unique prefixes of input fastq basenames which can be expand to full file names

|        |     |            |     |        |
| ---    | --- | ---        | --- | ---    |
| wt1    | wt  | single-end | N1  | female |
| wt2    | wt  | single-end | N2  | male   |
| trA_1  | A   | single-end | N1  | female |
| trA_2  | A   | single-end | N2  | male   |
| trB.n1 | B   | single-end | N1  | female |
| trB.n2 | B   | single-end | N2  | male   |


## Adapter sequences

Sequences can be found in the Illumina Adapter Sequences Document (<https://www.illumina.com/search.html?q=Illumina Adapter Sequences Document>) or Illumina Adapter Sequences HTML (<https://support-docs.illumina.com/SHARE/adapter-sequences.htm>) and the resource of Trimmomatic (<https://github.com/usadellab/Trimmomatic/tree/main/adapters>), FastQC respectively (<https://github.com/s-andrews/FastQC/blob/master/Configuration>).

The following excerpt is independent of the indexing type, i.e. single, unique dual (UD) or combinatorial dual (CD).

Nextera (Transposase Sequence), TruSight, AmpliSeq, stranded total/mRNA Prep, Ribo-Zero Plus: CTGTCTCTTATACACATCT

TruSeq (Universal) Adapter with A prefix due to 3' primer A-tailing : AGATCGGAAGAGC

TruSeq full length DNA & RNA R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq full length DNA MethC R1: AGATCGGAAGAGCACACGTCTGAAC R2: AGATCGGAAGAGCGTCGTGTAGGGA

TruSeq Small RNA: TGGAATTCTCGGGTGCCAAGG

## Example

Tiny example pipeline to perform gene fusion detection and differential expression analyses. Check out further pre-compiled pipelines for peak calling from *IP-Seq experiments and differential expression- and ontology analysis from RNA-Seq data (rippchen) or for multiple variant calling options from Exome-Seq/WG-Seq or RNA-Seq data following GATK best-practices in an optimized, parallelized fashion (muvac).

```bash
source <path/of/installation/latest/bashbone/activate.sh> -c true
genome=<path/to/fasta>
genomeidx=<path/to/segemehl.idx>
gtf=<path/to/gtf>
declare -a comparisons=($(ls <path/to/sample-info/files*>))
declare -a fastqs=($(ls <path/to/fastq-gz/files*>))
declare -a adapters=(AGATCGGAAGAGC)
declare -a qualdirs
fragmentsize=200
threads=4
memory=64000
preprocess::fastqc -t $threads -o results/qualities/raw -p /tmp -1 fastqs
preprocess::add4stats -r qualdirs -a results/qualities/raw -1 fastqs
preprocess::trimmomatic -t $threads -1 fastqs
preprocess::add4stats -r qualdirs -a results/qualities/trimmed -1 fastqs
preprocess::rmpolynt -t $threads -o results/polyntclipped -1 fastqs
preprocess::add4stats -r qualdirs -a results/qualities/polyntclipped -1 fastqs
preprocess::cutadapt -t $threads -o results/adapterclipped -a adapters -1 fastqs
preprocess::add4stats -r qualdirs -a results/qualities/adapterclipped -1 fastqs
preprocess::rcorrector -t $threads -o results/corrected -p /tmp -1 fastqs
preprocess::sortmerna -t $threads -o results/rrnafiltered -p /tmp -1 fastqs
preprocess::add4stats -r qualdirs -a results/qualities/rrnafiltered -1 fastqs
preprocess::qcstats -i qualdirs -o results/stats -p /tmp -1 fastqs
fusions::arriba -t $threads -g $genome -a $gtf -o results/fusions -p /tmp -f $fragmentsize -1 fastqs
fusions::starfusion -t $threads -g $genome -g $gtf -o results/fusions -1 fastqs
declare -a mapped
alignment::segemehl -t $threads -g $genome -x $genomeidx -o results/mapped -1 fastqs -r mapped
alignment::add4stats -r mapped
alignment::postprocess -j uniqify -t $threads -p /tmp -o results/mapped -r mapped
alignment::add4stats -r mapped
alignment::postprocess -r sort -t $threads -p /tmp -o results/mapped -r mapped
alignment::postprocess -r index -t $threads -p /tmp -o results/mapped -r mapped
alignment::qcstats -t $threads -o results/stats -r mapped
declare -A strandness
alignment::inferstrandness -t $threads -g $gtf -p /tmp -r mapped -x strandness
quantify::featurecounts -t $threads -p /tmp -g $gtf -o results/counted -r mapped -x strandness
quantify::tpm -t $threads -g $gtf -o results/counted -r mapped
expression::diego -t $threads -g $gtf -c comparisons -i results/counted -p /tmp -o results/diffexonjunctions -r mapped -x strandness
expression::deseq -t $threads -g $gtf -c comparisons -i results/counted -o results/diffgenes -r mapped
cluster::coexpression -t $threads -M $memory -g $gtf -b protein_coding -f 02 -i results/counted -o results/coexpressed -r mapped
```

# Third-party software

## In production

| Tool | Source | DOI |
| ---  | ---    | --- |
| Arriba        | <https://github.com/suhrig/arriba/>                                 | NA |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BWA           | <https://github.com/lh3/bwa>                                        | 10.1093/bioinformatics/btp324 |
| BWA-mem2      | <https://github.com/bwa-mem2/bwa-mem2>                              | 10.1109/IPDPS.2019.00041 |
| BWA-meth      | <https://github.com/brentp/bwa-meth>                                | arXiv:1401.1129 |
| BCFtools      | <http://www.htslib.org/doc/bcftools.html>                           | 10.1093/bioinformatics/btr509 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| bgztail       | <https://github.com/circulosmeos/bgztail>                           | NA |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler>        | 10.1089/omi.2011.0118 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 |
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656  |
| fgbio         | <http://fulcrumgenomics.github.io/fgbio/>                           | NA  |
| freebayes     | <https://github.com/ekg/freebayes>                                  | arXiv:1207.3907 |
| GATK          | <https://github.com/broadinstitute/gatk>                            | 10.1101/gr.107524.110 <br> 10.1038/ng.806 |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/kundajelab/idr>                                 | 10.1214/11-AOAS466 |
| IGV           | <http://software.broadinstitute.org/software/igv>                   | 10.1038/nbt.1754 |
| Intervene     | <https://github.com/asntech/intervene>                              | 10.1186/s12859-017-1708-7 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| metilene      | <https://www.bioinf.uni-leipzig.de/Software/metilene/>              | 10.1101/gr.196394.115 |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Platypus      | <https://rahmanteamdevelopment.github.io/Platypus>                  | 10.1038/ng.3036 |
| RAxML         | <https://cme.h-its.org/exelixis/web/software/raxml/index.html>      | 10.1093/bioinformatics/btl446 |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| ReSeqC        | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| SAMtools      | <http://www.htslib.org/doc/samtools.html>                           | 10.1093/bioinformatics/btp352 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| STAR-Fusion   | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>                   | 10.1101/120295 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| UMI-tools     | <https://github.com/CGATOxford/UMI-tools>                           | 10.1101/gr.209601.116 |
| VarDict       | <https://github.com/AstraZeneca-NGS/VarDict>                        | 10.1093/nar/gkw227 |
| VarScan       | <http://dkoboldt.github.io/varscan>                                 | 10.1101/gr.129684.111 |
| vcflib        | <https://github.com/vcflib/vcflib>                                  | NA |
| Vt            | <https://genome.sph.umich.edu/wiki/Vt>                              | 10.1093/bioinformatics/btv112 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |

## In preparation

| Tool | Source | DOI |
| ---  | ---    | --- |
| HISAT2          | <https://daehwankimlab.github.io/hisat2>                   | 10.1038/nmeth.3317 |
| SnpEff          | <https://pcingola.github.io/SnpEff>                        | 10.4161/fly.19695 |


# Supplementary information

Bashbone functions can be executed in parallel instances and thus are able to be submitted as jobs into a queuing system like a Sun Grid Engine (SGE). This could be easily done by using scripts written via here-documents or via the bashbone builtin `commander::qsubcmd`. The latter makes use of array jobs, which enables to wait for completion of all jobs, handle single exit codes and amend used resources via `qalter -tc <instances> <jobname>`.

```bash
source <path/of/installation/latest/muvac/activate.sh>
declare -a cmds=()
for i in *R1.fastq.gz; do
	j=${i/R1/R2}
	sh=job_$(basename $i .R1.fastq.gz)
	commander::makecmd -a cmd1 -c {COMMANDER[0]}<<- CMD
		bashbone -h
	CMD
done
commander::qsubcmd -r -l h="<hostname>|<hostname>" -p <env> -t <threads> -i <instances> -n <jobname> -o <logdir> -a cmds
# analogously: echo job.\$SGE_TASK_ID.sh | qsub -sync n -pe <env> <threads> -t 1-<#jobs> -tc <instances> -l h="<hostname>|<hostname>" -S /bin/bash -N <jobname> -o <logfile> -j y -V -cwd
```

In some cases a glibc pthreads bug (<https://sourceware.org/bugzilla/show_bug.cgi?id=23275>) may cause pigz failures (`internal threads error`) and premature termination of tools leveraging on it e.g. Cutadapt. One can circumvent this by e.g. making use of an alternative pthreads library via `LD_PRELOAD`

```bash
source <path/of/installation/latest/bashbone/activate.sh>
LD_PRELOAD=/lib64/noelision/libpthread.so.0 bashbone -h
LD_PRELOAD=/gsc/biosw/src/glibc-2.32/lib/libpthread.so.0 bashbone -h
```

# Closing remarks

Bashbone is a continuously developed library and actively used in my daily work. As a single developer it may take me a while to fix errors and issues. Feature requests cannot be handled so far, but I am happy to receive pull request.
