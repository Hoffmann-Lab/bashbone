# Bashbone

A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

For developers
- Write one-liners in your favorite language as Here-documents for later execution
- Add object-oriented programming (oop) like syntactic sugar to bash variables and arrays to avoid complex parameter-expansions, variable-expansion and brace-expansions
- Execute commands in parallel on your machine or submit them as jobs to a sun grid engine (SGE)
- Infer number of parallel instances according to targeted memory consumption or targeted thread per instance

For users
- Easily design multi-threaded pipelines to perform NGS related tasks
- For model AND non-model organisms
- Availability of many best-practice parameterized software wrappers
- Most software related parameters will be inferred directly from your data so that all functions have a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline
- Profit from many run-time tweaks

Covered Workflows for paired-end and single-end derived raw or prior mapped read data
- Data preprocessing (quality check, adapter clipping, quality trimming, error correction, artificial rRNA depletion)
- Read alignment and postprocessing
	- knapsack problem based slicing of alignment files for parallel task execution
	- sorting, filtering, unique alignment extraction, removal of optical duplicates
	- generation of pools and pseudo-replicates
	- read group modification, split N-cigar reads, left-alignment and base quality score recalibration
- Read quantification, TPM and Z-score normalization (automated inference of strand specific library preparation methods)
- Inference of differential expression as well as co-expression clusters
- Detection of differential splice junctions and differential exon usage
- Gene onotlogy (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering
- Free implementation of Encode3 best-practice ChIP-Seq Peak calling (automated inference of effective genome sizes)
- Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related *IP-Seq data
- Germline and somatic variant detection from DNA or RNA sequencing experiments plus VCF normalization

# License
The whole project is licensed under the GPL v3 (see LICENSE file for details)
**except** the the third-party tools set-upped during installation
Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download

```
git clone --recursive https://github.com/koriege/bashbone.git
git checkout $(git describe --tags)
```

# Quick start without installation - i.e. most functions will not work

Load the library and list available functions
Each function comes with a usage
```
source activate.sh
bashbone
```

## developers centerpiece

To create commands, print them and to execute them in parallel (according to your resources), utilize
```
commander::makecmd
commander::printcmd
commander::runcmd
# or
commander::qsubcmd
```

Example using the builtin array COMMANDER
```
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
commander::runcmd -v -b -t $threads -a cmds

```

oop bash
```
x=()
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

```
setup -i all -d <path/to/installation>
source <path/of/installation/activate.sh>
bashbone
```

## Update to a new release
```
setup -i upgrade -d <path/of/installation>
source <path/of/installation/activate.sh>
bashbone
```

# Usage

In order to perform desired comparisons, some functions require a sample info file
Assume this input:
- path/to/wt1.fq path/to/wt2.fq
- path/to/trA_1.fq path/to/trA_2.fq 
- path/to/trB.n1.fq path/to/trB.n2.fq
And this desired output:
- wt_vs_A wt_vs_b A_vs_B (N=2 vs N=2 each)
Then the info file should have:
- at least 4 columns (<name>, <main-factor>, single-end|paired-end, <replicate>) plus optional, additional factors
- unique prefixes of input fastq basenames in the first column
```
wt1      wt   single-end   N1   PE   female
wt2      wt   single-end   N2   SE   male
trA_1    A    single-end   N1   PE   female
trA_2    A    single-end   N2   PE   male
trB.n1   B    single-end   N1   SE   female
trB.n2   B    single-end   N2   PE   male
```

Example for differential gene expression analysis
```
source <path/of/installation/activate.sh> -c true
genome=<path/to/fasta>
genomeidx=<path/to/segemhl.idx>
gtf=<path/to/gtf>
comparisons=($(ls <path/to/sample-info/files*>))
fastqs=($(ls <path/to/fastq-gz/files*>))
adapters=(AGATCGGAAGAGC)
threads=4
memory=8000
memPerInstance=2000
preprocess::fastqc -t $threads -o results/qualities/raw -p /tmp -1 fastqs
preprocess::trimmomatic -t $threads -1 fastqs
preprocess::fastqc -t $threads -o results/qualities/trimmed -p /tmp -1 fastqs
preprocess::cutadapt -t $threads -o results/clipped -a adapters -1 fastqs
preprocess::fastqc -t $threads -o results/qualities/clipped -p /tmp -1 fastqs
preprocess::rcorrector -t $threads -o results/corrected -p /tmp -1 fastqs
preprocess::sortmerna -t $threads -o results/rrnafiltered -p /tmp 
preprocess::fastqc -t $threads -o results/qualities/rrnafiltered -p /tmp -1 fastqs
qualdirs=($(ls -d results/qualities/*/))
preprocess::qcstats -i qualdirs -o results/stats -p /tmp -1 fastqs
mapped=()
alignment::segemehl -t $threads -g $genome -x $genomeidx -o results/mapped -1 fastqs -r mapped
alignment::add4stats -r mapper
alignment::postprocess -j uniqify -t $threads -p /tmp -o results/mapped -r mapped
alignment::add4stats -r mapper
alignment::postprocess -r sort -t $threads -p /tmp -o results/mapped -r mapped
alignment::postprocess -r index -t $threads -p /tmp -o results/mapped -r mapped
alignment::bamstats -t $threads -o results/stats -r mapped
quantify::featurecounts -t $threads -p /tmp -g $gtf -o results/counted -r mapped
quantify::tpm -t $threads -g $gtf -o results/counted -r mapped
expression::deseq -t $threads -g $gtf -c comparisons -i results/counted -o results/deseq -r mapped
```

# Third-party software

| Tool        | Source                                              | doi/pmid |
| ----------- | --------------------------------------------------- | -------- |
| segemehl    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |

