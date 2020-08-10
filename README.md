# Bashbone
---

A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

## For developers

- Write one-liners in your favorite language as Here-documents for later execution
- Add object-oriented programming (oop) like syntactic sugar to bash variables and arrays to avoid complex parameter-expansions, variable-expansions and brace-expansions
- Execute commands in parallel on your machine or submit them as jobs to a sun grid engine (SGE)
- Infer number of parallel instances according to targeted memory consumption or targeted threads per instance

## For users

- Easily design multi-threaded pipelines to perform NGS related tasks
- For model AND non-model organisms
- Availability of many best-practice parameterized and run-time tweaked software wrappers
- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks

- For paired-end and single-end derived raw sequencing or prior mapped read data
- Data preprocessing (quality check, adapter clipping, quality trimming, error correction, artificial rRNA depletion)
- Read alignment and post-processing
  - knapsack problem based slicing of alignment files for parallel task execution
  - sorting, filtering, unique alignment extraction, removal of optical duplicates
  - generation of pools and pseudo-replicates
  - read group modification, split N-cigar reads, left-alignment and base quality score recalibration
- Read quantification, TPM and Z-score normalization (automated inference of strand specific library preparation methods)
- Inference of differential expression as well as co-expression clusters
- Detection of differential splice junctions and differential exon usage
- Gene ontology (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering
- Free implementation of Encode3 best-practice ChIP-Seq Peak calling (automated inference of effective genome sizes)
- Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related *IP-Seq data
- Germline and somatic variant detection from DNA or RNA sequencing experiments plus VCF normalization

# License
---

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download
---

```bash
git clone --recursive https://github.com/koriege/bashbone.git
git checkout $(git describe --tags)
```

# Quick start (without installation)
Please see installation section to get all third-party tools set-upped and subsequently all functions to work properly.

Load the library and list available functions. Each function comes with a usage.

```bash
source activate.sh
bashbone
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

Example using the builtin file descriptors array `COMMANDER` with seperator (`-s '|'`) based concatenation of interpreted and non-interpreted Here-Documents. Each concatenation will be stored as a single command in a de-referenced array (`-a cmds`) Optionally, an output file can be defined (`-o <path/to/file>`).

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
commander::runcmd -v -b -t $threads -a cmds

```

### OOP bash

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
---

```
setup -i all -d <path/to/installation>
source <path/of/installation/activate.sh>
bashbone
```

## Update to a newer release

```
setup -i upgrade -d <path/of/installation>
source <path/of/installation/activate.sh>
bashbone
```

# Usage
---

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
- wt_vs_b
- A_vs_B 

Then the info file should consist of:

- At least 4 columns (`<name>`, `<main-factor>`, `single-end|paired-end`, `<replicate>`)
- Optionally, additional factors
- Unique prefixes of input fastq basenames in the first column which expand to the full file name

|        |     |            |     |        |
| ---    | --- | ---        | --- | ---    |
| wt1    | wt  | single-end | N1  | female |
| wt2    | wt  | single-end | N2  | male   |
| trA_1  | A   | single-end | N1  | female |
| trA_2  | A   | single-end | N2  | male   |
| trB.n1 | B   | single-end | N1  | female |
| trB.n2 | B   | single-end | N2  | male   |


## Example

Tiny example pipeline to perform differential gene expression analysis

```bash
source <path/of/installation/activate.sh> -c true
genome=<path/to/fasta>
genomeidx=<path/to/segemhl.idx>
gtf=<path/to/gtf>
comparisons=($(ls <path/to/sample-info/files*>))
fastqs=($(ls <path/to/fastq-gz/files*>))
adapters=(AGATCGGAAGAGC)
threads=4
memoryPerInstance=2000
preprocess::fastqc -t $threads -o results/qualities/raw -p /tmp -1 fastqs
preprocess::trimmomatic -t $threads -1 fastqs
preprocess::fastqc -t $threads -o results/qualities/trimmed -p /tmp -1 fastqs
preprocess::cutadapt -t $threads -o results/clipped -a adapters -1 fastqs
preprocess::fastqc -t $threads -o results/qualities/clipped -p /tmp -1 fastqs
preprocess::rcorrector -t $threads -o results/corrected -p /tmp -1 fastqs
preprocess::sortmerna -t $threads -m $memoryPerInstance -o results/rrnafiltered -p /tmp 
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
---

## Production
| Tool        | Source                                              | doi/pmid |
| --- | --- | --- |
| segemehl    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| STAR    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Trimmomatic    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| SortMeRNA    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Rcorrector    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| DEXSeq    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| WGCNA    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| DGCA    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| REVIGO    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| GEM    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| IDR    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Macs2    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| fastqc    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| khmer    | <https://khmer.readthedocs.io/en/v2.1.1/index.html> | 10.12688/f1000research.6924.1 |
| bcftools    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| GATK    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| vcflib    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| vt    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| bamutil    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| bedtools    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Picard    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| featureCounts    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| DESeq2    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Cutadapt    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| ReSeqC    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| HTSeq    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| DIEGO    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| rngtools    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| gseabase    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| clusterprofiler    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |

## To implement
| Tool        | Source                                              | doi/pmid |
| --- | --- | --- |
| STAR-Fusion    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| BWA    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Hisat2    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Freebayes    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Varscan    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| Platypus    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| VarDict    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |
| snpeff    | <http://www.bioinf.uni-leipzig.de/Software/segemehl/> | 19750212, 24512684 |




