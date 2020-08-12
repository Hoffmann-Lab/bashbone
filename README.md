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

## In production

| Tool | Source | DOI |
| ---  | ---    | --- |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BCFtools      | <http://www.htslib.org/doc/bcftools.html>                           | 10.1093/bioinformatics/btr509 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 |
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656  |
| GATK          | <https://github.com/broadinstitute/gatk>                            | 10.1101/gr.107524.110 <br> 10.1038/ng.806 |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/kundajelab/idr>                                 | 10.1214/11-AOAS466 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| ReSeqC        | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| vcflib        | <https://github.com/vcflib/vcflib>                                  | NA |
| Vt            | <https://genome.sph.umich.edu/wiki/Vt>                              | 10.1093/bioinformatics/btv112 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |


## In preparation

| Tool | Source | DOI |
| ---  | ---    | --- |
| BWA             | <https://github.com/lh3/bwa>                               | 10.1093/bioinformatics/btp324 |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler> | NA |
| freebayes       | <https://github.com/ekg/freebayes>                         | arXiv:1207.3907 |
| HISAT2          | <https://daehwankimlab.github.io/hisat2>                   | 10.1038/nmeth.3317 |
| Platypus        | <https://rahmanteamdevelopment.github.io/Platypus>         | 10.1038/ng.3036 |
| SnpEff          | <https://pcingola.github.io/SnpEff>                        | 10.4161/fly.19695 |
| STAR-Fusion     | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>          | 10.1101/120295 |
| VarDict         | <https://github.com/AstraZeneca-NGS/VarDict>               | 10.1093/nar/gkw227 |
| VarScan         | <http://dkoboldt.github.io/varscan>                        | 10.1101/gr.129684.111 |






