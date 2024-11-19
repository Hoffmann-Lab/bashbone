# Bashbone

A bash and biobash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

# Outline

- [Bashbone](#bashbone)
  - [For developers - bash library](#for-developers---bash-library)
  - [For users - biobash library](#for-users---biobash-library)
    - [Covered Tasks](#covered-tasks)
- [License](#license)
- [Download](#download)
- [Bash library usage (without full installation)](#bash-library-usage-without-full-installation)
  - [Do's and don'ts](#dos-and-donts)
  - [Quick start](#quick-start)
  - [Developers centerpiece](#developers-centerpiece)
    - [Example](#example)
    - [Cleanup](#cleanup)
  - [Helper functions](#helper-functions)
    - [OOP bash](#oop-bash)
    - [Multithreaded implementations](#multithreaded-implementations)
    - [Misc](#misc)
  - [Extending the library](#extending-the-library)
    - [Local cleanup](#local-cleanup)
    - [Global cleanup](#global-cleanup)
- [Installation](#installation)
  - [Full installation of all third party tools used in bashbones biobash library](#full-installation-of-all-third-party-tools-used-in-bashbones-biobash-library)
  - [Upgrade to a newer release (sources only)](#upgrade-to-a-newer-release-sources-only)
  - [Update tools](#update-tools)
- [Biobash library usage (requires installation)](#biobash-library-usage-requires-installation)
  - [Enclosed scripts](#enclosed-scripts)
    - [Retrieve SRA datasets](#retrieve-sra-datasets)
    - [Retrieve genomes](#retrieve-genomes)
    - [Merge/collate features](#mergecollate-features)
    - [Generate pseudo replicates](#generate-pseudo-replicates)
    - [Index and random access flat files](#index-and-random-access-flat-files)
  - [Sample info file](#sample-info-file)
  - [Adapter sequences](#adapter-sequences)
  - [Pipeline example](#pipeline-example)
- [Third-party software](#thirdparty-software)
- [Supplementary information](#supplementary-information)
- [Closing remarks](#closing-remarks)

## For developers - bash library
[&#x25B2; back to top](#bashbone)

- Get a full bash error stack trace in interactive shells or within scripts
- Write command line code in your favorite programming language via Here-documents for later orchestrated execution
- Add object-oriented programming (oop) like syntactic sugar to bash variables and arrays to avoid complex parameter-expansions, variable-expansions and brace-expansions
- Execute commands in parallel on your machine or submit them as jobs to a workflow manager like sun grid engine (SGE) and log stdout, stderr and exit codes per job
- Benchmark runtime and memory usage
- Infer number of parallel instances according to targeted memory consumption or targeted threads per instance
- Log execution of bash functions at different verbosity levels
- Extend the library by custom bash functions which will inherit
  - Stack trace
  - Termination of all function related (sub-)processes, including asynchronous background jobs upon error/exit or when reaching prompt-command (interactive shell)
  - Removal of temporary files created via `mktemp` and execution of custom cleanup commands upon error/exit or when reaching prompt-command (interactive shell)
- Profit from helper functions that implement
  - Multi-threaded joining of multiple files
  - Multi-threaded sorting
  - Multi-threaded de-compression
  - Multi-threaded compression plus indexing for random access by byte offset or line number without noticeable overhead
  - Multi-threaded application of commands on an compressed and indexed or flat file on per-line or per-chunk basis

## For users - biobash library
[&#x25B2; back to top](#bashbone)

- Get a full bash error stack trace in interactive shells or within scripts
- Easily design multi-threaded pipelines to perform NGS related tasks
- Use many best-practice parameterized and heavily run-time tweaked software wrappers
- Most software related parameters will be inferred directly from input data, so that all functions require just a minimal set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)

### Covered tasks
[&#x25B2; back to top](#bashbone)

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
  - sorting, filtering
  - UMI based de-duplication or removal of optical and PCR duplicates
  - generation of pools and pseudo-replicates
  - read group modification, split N-cigar reads, left-alignment and base quality score recalibration
- Gene fusion detection
- Methyl-C calling and prediction of differentially methylated regions
- Expression analysis
  - Read quantification (also from quasi-mappings), TPM and Z-score normalization and heatmap plotting
  - Inference of strand specific library preparation methods
  - Inference of differential expression as well as clusters of co-expression
  - Detection of differential splice junctions and differential exon usage
  - Gene ontology (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering and visualizations
- Implementation of ENCODE v3 best-practice ChIP-Seq Peak calling
  - Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related IP-Seq data
  - Inference of effective genome sizes
- Variant detection from DNA or RNA sequencing experiments
  - Integration of multiple solutions for germline and somatic calling
  - VCF normalization
  - Tree reconstruction from homozygous sites
- ssGSEA and survival analysis from TCGA cancer expression data
- Genome and SRA data retrieval
  - Genome to transcriptome conversion
  - Data visualization via IGV batch processing

# License
[&#x25B2; back to top](#bashbone)

The whole project is licensed under the GPL v3 (see LICENSE file for details), **except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download
[&#x25B2; back to top](#bashbone)

This will download you a copy which includes the latest developments

```bash
git clone --recursive https://github.com/Hoffmann-Lab/bashbone
```

To check out the latest release (**irregularly compiled**) do

```bash
cd bashbone
git checkout $(git describe --tags)
```

# Bash library usage (without full installation)
[&#x25B2; back to top](#bashbone)

## Do's and don'ts
[&#x25B2; back to top](#bashbone)

When used, in a script, bashbone is meant to be sourced at the very top to handle positional arguments and to re-execute (`-r true`) the script under its own process group id in order to take care of proper termination (`-a "$@"`). It will enable error stack tracing and sub-process handling globally by setting traps for `EXIT` `ERR` `RETURN` `INT`. So, don't override them. In case your script intends to spawn deamons use `setsid` or disable bashbone first.

```bash
#!/usr/bin/env bash
source <path/to/bashbone>/activate.sh -r true -a "$@"
# do stuff
# now spawn deamons
setsid deamon1 &
bashbone -x
deamon2 &
```

Please note, that error tracing in bash is circumvented by using `||` or '&&' constructs. Therefore, avoid them in any context of function calls.

```bash
#!/usr/bin/env bash
source <path/to/bashbone>/activate.sh -r true -a "$@"

function myfun(){
  cat file_not_found
  echo "error ignored. starting time consuming calculation now."
}
# DON'T !
myfun || echo "failed with code $?"
```

## Quick start
[&#x25B2; back to top](#bashbone)

To get all third-party tools set-upped and subsequently all biobash bashbone functions to work properly, **see also**

- [Installation](#installation)
- [Biobash library usage (requires installation)](#biobash-library-usage-requires-installation)

For a lite installation that gets you the minimum required tools (`GNU parallel`, `gztool`, `mdless`) in order to make use of developer functions, execute

```bash
scripts/setup.sh -i lite -d <path/to/installation>
```

To see the usage, do

```bash
scripts/setup.sh -h

DESCRIPTION
Bashbone setup routine

SYNOPSIS
setup.sh -i [all|upgrade] -d [path]

OPTIONS
-i | --install [lite|all|upgrade] : install into given directory
-g | --use-config                 : use supplied yaml files and URLs instead of cutting edge tools
-d | --directory [path]           : installation path
-t | --threads [value]            : threads - predicted default: 32
-l | --log [path]                 : log file - default: [-d]/install.log
-v | --verbose                    : enable verbose mode
-h | --help                       : prints this message

DEVELOPER OPTIONS
-s | --source [path,..]           : source file(s) to overload compile::[lite|all|upgrade|<tool>] functions
-i | --install [<tool>,..]        : install into given directory
```

Now load the bashbone library in an interactive terminal session. Note, that **none of your environment settings will be modified**.

```bash
source <path/of/installation>/latest/bashbone/activate.sh
```

To see the activate script usage, do

```bash
source <path/of/installation>/latest/bashbone/activate.sh -h

This is bashbone activation script.

To see lists of available options and functions, source me and execute bashbone -h

Usage:
-h              | this help
-l <legacymode> | true/false let commander inerts line breaks, thus crafts one-liners from makecmd here-documents
                  default: false
-i <path>       | to installation root <path>/latest
                  default: inferred from script location
                  hint: run activation from source code, indeed enables basic functions, but will fail on executing tools
-p <path>       | to temporary directory. default: $TMPDIR or /tmp
-c <activate>   | true/false conda from [-i]/conda/bin and prepend bashbone installed tools to PATH
                  default: false
-x <string>     | a function or command string to be called upon EXIT signal. this function or command will receive the exit code as last positional argument
-s <path>       | extend/overload bashbone functions by libraries at <path>/lib and via -c option also add <path>/scripts to PATH
-r <re-execute> | true/false the script that implements activate.sh under a new process group to not bubble up INT/TERM or kill a pipeline upon exit
-a <optarg>     | use as last option in case bashbone is sourced in a script wich makes use of positional arguments (optargs)
```

Listing available functions, scripts, access to this documentation, and deactivation can be performed via

```bash
bashbone -h

Bashbone v1.4.0

Is a bash/biobash library for workflow and pipeline design within, but not restricted to, the scope of Next Generation Sequencing (NGS) data analyses.

Usage:
-h        | this help
-l        | enable or disable legacy mode
-r        | open readme
-c        | activate bashbone conda environment or deactivate conda
-s        | list bashbone scripts
-f <name> | list all functions or if <name> is supplied, execute it
-d        | list bashbone developer functions
-e        | list installed tools and versions
-x        | exit bashbone i.e. remove functions and disable conda. within scripts also restores environment i.e. traps and shell options

```

Each function, which can be auto-completed when using bash v5 or later, comes with its own usage. Without full installation, only developer functions will work and can be listed this way

```bash
bashbone -d

commander::makecmd              commander::print                 commander::printcmd  commander::printerr
commander::printinfo            commander::qalter                commander::qstat     commander::qsubcmd
commander::runalter             commander::runcmd                commander::runstat   commander::warn
configure::instances_by_memory  configure::instances_by_threads  configure::jvm       configure::memory_by_instances
helper::addmemberfunctions      helper::basename                 helper::capply       helper::cat
helper::index                   helper::isarray                  helper::ishash       helper::join
helper::lapply                  helper::makecatcmd               helper::multijoin    helper::pgzip
helper::ps2pdf                  helper::sort                     helper::vcfsort      progress::log
```

Now, for example, inspect the usage of
```bash
helper::sort

helper::sort wraps GNU sort in a parallelized fashion as drop-in replacement
-f <infile>    | optional. default: stdin
               | path to input file
-o <outfile>   | optional. default: stdout
               | path to output file
-t <threads>   | mandatory
               | number of threads
-M <maxmemory> | optional. default: all available
               | amount of memory in megabytes to allocate
```


To unset bashbone functions in your interactive shell or to revert changes made to the environment when used in a script, do

```bash
bashbone -x
```

When bashbone should be used within a script, which may uses positional arguments, hand them over to bashbone activation script by adding this snipped to the very top of your script

```bash
#! /usr/bin/env bash
source <path/of/installation>/latest/bashbone/activate.sh -r true -a "$@"
```

## Developers centerpiece
[&#x25B2; back to top](#bashbone)

To create commands, print, and to execute them in parallel (with optional benchmarking of runtime and memory consumption), further inspect and modify task concurrency, utilize

```bash
commander::makecmd
commander::printcmd
commander::runcmd
commander::runstat
commander::runalter

# truncated usage of commander::makecmd

easily compose and store command in an array for later execution

-a <cmds>          | mandatory. array to append the command onto
-v <variable>      | optional. env variable(s) necessary for command. can be used multiple times
-s <separator>     | optional. character to separate command chunks, when file descriptor array COMMANDER is used. default: '\n'
-m                 | force multi line support in case of activated bashbone legacy mode, which inerts line breaks
-o <outfile>       | optional. redirect stdout to output file
-O <outfile>       | optional. append stdout to output file
-c <cmd|{fd[0]}..> | ALWAYS LAST OPTION. command line string(s) and or here-documents or file descriptor array COMMANDER


# usage commander::runcmd

orchestrates the execution of commands stored within an array by utilizing GNU parallel

-v             | optional. toggle on verbose mode which will print the commands prior to their execution
-b             | optional. toggle on benchmarking run time and peak memory consumption
-c <env>       | optional. activate a bashbone conda environment prior to command execution
-i <instances> | optional. number of parallel instances. see also commander::runalter. default: 1
-t <instances> | optional. obsolete synonym for -i
-s <idx[:idx]> | optional. restrict the execution of commands by array index start or range. default: 1
-n <name>      | optional. prefix for logs and job scripts
-o <path>      | optional. output directory for logs, job scripts and exit codes
-r             | optional. override existing logs. default: append
-a <cmds>      | mandatory. array of commands. see also commander::makecmd

example 1:
cmds=("echo job1" "echo job2" "echo job3")
cmds+=("sleep 2; echo job4")
commander::runcmd -v -b -i 2 -a cmds

example 2:
commander::runcmd -i 2 -n current -o ~/jobs -r -a cmds
```

Parallel task execution is also possible by submitting commands as an array job to the workload manager Sun Grid Engine or Son of Grid Engine fork (SGE) via the in-line replacement functions

```bash
commander::qsubcmd
commander::qstat
commander::qalter
```

To monitor any bash function and redirect all stdout and stderr into a log file while defining what is going to be printed to the terminal via verbosity levels (`-v [1..3]`), use

```bash
progress::log -v 1 -o <logfile> -f commander::runcmd
```

To infer a suitable number of parallel instances for a local execution of tasks given the current state of CPU and memory resources, call bashbone configuration functions. Settings for memory and threads per instance and number of instances as well as garbage collection threads and concurrency for java jvm, matching the current machine (unless the `-d` switch is used to perform a dry run), will be returned.

The execution of applications known to have a large memory footprint should be parameterized according to

```bash
configure::instances_by_memory
configure::memory_by_instances
configure::jvm
```

Task parallelization settings for lightweight applications, can be inferred via

```bash
configure::instances_by_threads
```

### Example
[&#x25B2; back to top](#bashbone)

For a later orchestrated, parallel execution of commands, they need to be stored beforehand, which typically calls in for escaping special characters to not get them interpreted by the shell. `commander::makecmd` offers two methods to solve this by using interpreted (`EOF`) or non-interpreted (`'EOF'`) Here-documents and file descriptor arrays.

```bash
for i in sun mercury venus earth mars jupiter saturn uranus neptune pluto; do
  echo "$i" | awk '{print "hello "$1}'
done
# store
declare -a cmds
for i in sun mercury venus earth mars jupiter saturn uranus neptune pluto; do
  cmds+=("echo \"$i\" | awk '{print \"hello \"\$1}'")
done
```

**Solution 1** relies on variable reference(s) via `-v <var>` and an uninterpreted Here-document. The content of each Here-document will be stored as a single command string in a de-referenced array (`-a cmds`).

```bash
declare -a cmds
for i in sun mercury venus earth mars jupiter saturn uranus neptune pluto; do
  commander::makecmd -a cmds -v i -c <<-'EOF'
    echo "$i" | awk '{print "hello "$1}'
  EOF
done
```

**Solution 2** showcase the builtin file descriptors array `COMMANDER` from which the content of interpreted or non-interpreted Here-documents will be read and concatenated using the pipe symbol as separator via `-s '|'`. Each concatenation will be stored as a single command string in a de-referenced array (`-a cmds`). Optionally, an output file can be defined (`-o <path/to/file>`).

```bash
declare -a cmds
for i in sun mercury venus earth mars jupiter saturn uranus neptune pluto; do
	commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<-EOF {COMMANDER[1]}<<-'EOF'
		echo $i
	EOF
		awk '{
			print "hello "$1
		}'
	EOF
done
```

The stored commands can now be printed and executed in parallel instances. Optionally, stdout and sterr of each instance as well as its exit code can be re-directed into log files via `-o <logdir> -n <prefix>`. 

```bash
commander::printcmd -a cmds
instances=4
commander::runcmd -v -b -i $instances -a cmds
```

### Cleanup
[&#x25B2; back to top](#bashbone)

When commands are crafted for later execution, bashbone allows for multiple cleanup strategies to be applied on exit (success or failure). A **reversely** executed, temporary script, which is accessible through the `$BASHBONE_CLEANUP` variable. An overloaded `mktemp` function, that automatically adds removal commands to the cleanup script. An after all called `_on_exit` function, which holds the job scripts exit code in the first positional argument.

```bash
source <path/to/bashbone>/activate.sh

commander::makecmd -a cmds -c <<-'EOF'
  function _on_exit(){
    rm -r /tmp/tmptest3
    echo "exit code is $1"
  }
  mkdir -p /tmp/tmptest1
  echo "rm -r /tmp/tmptest1" >> "$BASHBONE_CLEANUP"
  mktemp -p /tmp/tmptest2
  mkdir -p /tmp/tmptest3
  ls /tmp/tmptest*
EOF

commander::runcmd -a cmds -i 1

ls: cannot access '/tmp/tmptest*': No such file or directory
:ERROR: in job.1.sh line 8: ls /tmp/tmptest*
```

## Helper functions
[&#x25B2; back to top](#bashbone)

### OOP bash
[&#x25B2; back to top](#bashbone)

A small excerpt of possible member(-like-)functions, which help to avoid complex parameter expansions. See `bashbone -d` helper underscore functions for a full list.

```bash
declare -a arr
helper::addmemberfunctions -v arr

arr.push "hello" "world" "and" "moon"
arr.get 1 # world
arr.get -1 # moon
arr.get 0 1 # hello world
arr.get 1 -1 # world and
arr.print # hello world and moon
arr.shift # world and moon
arr.pop # world and
arr.print # world and
arr.substring 2 4 # rld d
arr.sort # d rld
arr.uc # D RLD
arr.print # D RLD
arr.length # 2
# see more available arr.* member functions via auto-completion
```

### Multithreaded implementations
[&#x25B2; back to top](#bashbone)

To execute GNU sort in a multi-threaded way, either as in-line replacement or adjusted to work for VCF files, use 

```bash
helper::sort
helper::vcfsort
```

In order to achieve fast gzip compression including indexing for random access via line number, utilize

```bash
helper::pgzip
```

When a gzip file was not created via `helper::pgzip` and needs to be indexed or if a plain text flat file should be indexed, utilize `helper::index`. The index in turn can be used for random access in order to apply a command or pipeline on a per-chunk parallelized fashion. Each process, job respectively, can be addressed by its `JOB_ID`. The following example will access the file at 10 equally distributed positions in parallel to apply the `cat` command.

```bash
helper::index
helper::capply

# example
helper::capply -i 10 -f <file[.gz]> -c 'cat > $JOB_ID.out'
```

An other method, which does not require any index, allows to apply a command or pipeline parallelized on each line of a file supplied via stdin - either directly on the stream or by creating seekable temporary files. The number of total processes executed is the number of lines divided by number of records or a fixed number of instances to which the input is streamed round robin wise. Each process can be addressed by its `JOB_ID`. Temporary file names can be addressed using `FILE` variable. The next example shows 10 parallel running processes in which the `cat` command is applied on a temporary file built from a chunk of 10000 lines. After termination of a command, a new process will be spawned handling the next chunk.

```bash
helper::lapply

# example
cat <file> | helper::lapply -f -r 10000 -t 10 -c 'cat $FILE > $JOB_ID.out'
```

Ultra fast printing of a files content, thanks to multi-threaded decompression of gzip and bzip files, including mime-type identification, can be achieved using

```
helper::cat -f <file>
```

Do a parallelized full join of multiple sorted files by unique ids in the first column(s) or of positional sorted bed files, given a separator, header and NA character can be done via

```bash
helper::multijoin <file> <file> [<file> ..]
```

### Misc
[&#x25B2; back to top](#bashbone)

To get the basename of a file with stripped off extension (two extensions in case of compressed data), and optionally report the suffix as well, use

```bash
helper::basename -f <file> -o base -e ex
echo $base
echo $ex
```

If unsure about the data type of a variable, check it via exit/return codes of

```bash
declare -a arr=(1 2 3)
declare -A hash=([a]=1 [b]=2)
helper::isarray -v arr
helper::ishash -v hash
```

## Extending the library
[&#x25B2; back to top](#bashbone)

To get a full bash error stack trace for custom bash functions also when interactive, to add proper handling of sub-process termination and to run cleanup procedures upon error or when reaching script end, prompt-command respectively (interactive shell), simply extend bashbone functions by your own library.

I.e. store a set of bash functions that **require** the `function` keyword in a `lib` directory and files with a `.sh` suffix to be sourced along with the bashbone library afterwards.

```bash
cd <path/to/custom>
mkdir lib
cat <<-'EOF' > lib/fun.sh
  function world(){
    echo "hello world"
    mars
  }

  function mars(){
    echo "hello mars"
    cmd_that_fails
  }
EOF
```

Now hand over the parent directory to `activate.sh`.

```bash
source <path/to/bashbone>/activate.sh -s <path/to/custom>
world

# will print
hello world
hello mars
cmd_that_fails: command not found
:ERROR: in <path/to/custom>/lib/fun.sh (function: mars) @ line 8: cmd_that_fails
:ERROR: in <path/to/custom>/lib/fun.sh (function: world) @ line 3: mars
```

### Local cleanup
[&#x25B2; back to top](#bashbone)

Each bashbone procedure, function from a custom extension respectively, gets a local, temporary cleanup script assigned. This bash script will be **reversely** executed upon return (success or failure) of the function. The path to the script is kept in the `$BASHBONE_CLEANUP` variable. When `mktemp` is used within the function, a removal command of this path will be automatically added to the cleanup script. For sure, custom commands can be added, too.

```bash
cd <path/to/custom>
mkdir lib
cat <<-'EOF' > lib/fun.sh
  function tmptest(){
    mkdir -p /tmp/tmptest1
    echo "rm -r /tmp/tmptest1" >> "$BASHBONE_CLEANUP"
    mktemp -p /tmp/tmptest2
    ls /tmp/tmptest*
  }
EOF

source <path/to/bashbone>/activate.sh -s <path/to/custom>
tmptest

ls: cannot access '/tmp/tmptest*': No such file or directory
:ERROR: in <path/to/custom>/lib/tmptest.sh (function: mars) @ line 4: ls /tmp/tmptest*
```

As an alternative to sourcing files from a custom library directory, functions can be defined and wrapped on the fly. Therefore, define an alias like this **before** its actual definition.

```bash
alias myfun="_bashbone_wrapper myfun"

function myfun(){
  # do stuff
}

myfun
```

### Global cleanup
[&#x25B2; back to top](#bashbone)

Upon script exit, due to failure or upon success, a custom cleanup functions can be handed over to bashbone and executed at the very end. The scripts exit code will be supplied as positional argument. Analogous to [Local cleanup](#local-cleanup), temporary files created by `mktemp` will be automatically nuked and a cleanup script, kept in the `$BASHBONE_CLEANUP` variable, is executed.

```bash
#! /usr/bin/env bash
source <path/to/bashbone>/activate.sh -x cleanup -r true -a "$@"

function cleanup(){
  rm -r /tmp/tmptest3
  ls /tmp/tmptest* # should fail
  echo "exit code is $1"
}

mkdir -p /tmp/tmptest1
echo "rm -r /tmp/tmptest1" >> "$BASHBONE_CLEANUP"
mktemp -p /tmp/tmptest2
mkdir -p /tmp/tmptest3
ls /tmp/tmptest*
```

# Installation
[&#x25B2; back to top](#bashbone)

## Full installation of all third party tools used in bashbones biobash library
[&#x25B2; back to top](#bashbone)

When using the `-g` switch (**recommended**), the setup routine will create conda environments or setups software from source according to enclosed configuration files, URLs respectively. Without `-g` switch, software is installed in latest available version, which may lead to unexpected behavior and errors. During setup, current configuration files will be written to `<path/of/installation/config>`.

```bash
# shows usage
scripts/setup.sh -h

scripts/setup.sh -g -i all -d <path/to/installation>
```

When bashbone was successfully installed, it is recommended to activate bashbone along with its conda environment to ensure all biobash functions to operate as expected. The conda environment can be stopped and (re-)started by the `bashbone -c` switch.

```
source <path/of/installation>/latest/bashbone/activate.sh -c true
```

## Upgrade to a newer release
[&#x25B2; back to top](#bashbone)

Use the `-g` switch, in order to also upgrade biobash conda environments that fail the comparison with the supplied configuration files. **Attention**: This switch will downgrade tools, if the initial installation was done for cutting edge tools i.e. without `-g`.

```bash
scripts/setup.sh -g -i upgrade -d <path/of/installation>
```

## Update tools
[&#x25B2; back to top](#bashbone)

Trimmomatic, STAR-Fusion, GEM, mdless and gztool will be installed next to the biobash conda environments. Their latest versions and download URLs will be automatically inferred.

```bash
scripts/setup.sh -i trimmomatic,starfusion,gem,mdless,gztool -d <path/of/installation>
```

# Biobash library usage (requires installation)
[&#x25B2; back to top](#bashbone)

To get all third-party tools set-upped and subsequently all biobash bashbone functions to work properly, **see also**

- [Installation](#installation)
- [Do's and don'ts](#dos-and-donts)

Load the library and list available functions run

```bash
source <path/of/installation>/latest/bashbone/activate.sh

bashbone -f

alignment::add4stats     alignment::addreadgroup      alignment::bamqc           alignment::bqsr
alignment::bulkindex     alignment::bwa               alignment::clip            alignment::clipmateoverlaps
alignment::downsample    alignment::inferstrandness   alignment::leftalign       alignment::mkreplicates
alignment::postprocess   alignment::qcstats           alignment::reorder         alignment::rmduplicates
alignment::segemehl      alignment::slice             alignment::soft2hardclip   alignment::splitncigar
alignment::star          alignment::strandsplit       alignment::tn5clip         alignment::tobed
------------------------------------------------------------------------------------------------------------
bigwig::apply            bigwig::profiles
------------------------------------------------------------------------------------------------------------
bisulfite::bwa           bisulfite::haarz             bisulfite::join            bisulfite::mecall
bisulfite::methyldackel  bisulfite::metilene          bisulfite::mspicut         bisulfite::rmduplicates
bisulfite::segemehl      
------------------------------------------------------------------------------------------------------------
cluster::coexpression    cluster::coexpression_deseq  cluster::wgcna             cluster::wgcna_deseq
------------------------------------------------------------------------------------------------------------
enrichment:go
------------------------------------------------------------------------------------------------------------
expression::deseq        expression::diego            expression::join           expression::join_deseq
------------------------------------------------------------------------------------------------------------
fusions::arriba          fusions::join                fusions::starfusion
------------------------------------------------------------------------------------------------------------
genome::indexgtf         genome::mkdict               genome::mkgodb             genome::view
------------------------------------------------------------------------------------------------------------
peaks::gem               peaks::gem_idr               peaks::genrich             peaks::genrich_idr
peaks::gopeaks           peaks::gopeaks_idr           peaks::m6aviewer           peaks::m6aviewer_idr 
peaks::macs              peaks::macs_idr              peaks::matk                peaks::matk_idr
peaks::peakachu          peaks::peakachu_idr          peaks::seacr               peaks::seacr_idr           
------------------------------------------------------------------------------------------------------------
preprocess::add4stats    preprocess::cutadapt         preprocess::dedup          preprocess::fastqc
preprocess::qcstats      preprocess::rcorrector       preprocess::rmpolynt       preprocess::sortmerna
preprocess::trimmomatic   
------------------------------------------------------------------------------------------------------------
quantify::bamcoverage    quantify::featurecounts      quantify::normalize        quantify::profiles
quantify::salmon         quantify::tpm
------------------------------------------------------------------------------------------------------------
survival::gettcga        survival::ssgsea            
------------------------------------------------------------------------------------------------------------
variants::bcftools       variants::freebayes          variants::haplotypecaller  variants::makepondb
variants::mutect         variants::panelofnormals     variants::platypus         variants::tree
variants::vardict        variants::vardict_threads    variants::varscan          variants::vcfnorm
------------------------------------------------------------------------------------------------------------
visualize::venn
```

In order to make use of biobash conda environments, which ensures all supplied scripts to work as expected, activate bashbone with conda enabled

```bash
source <path/of/installation>/latest/bashbone/activate.sh -c true
```

Or activate/deactivate conda at a later timepoint

```bash
source <path/of/installation>/latest/bashbone/activate.sh

bashbone -c
```

Or use commander functions with conda enabled (see [Developers centerpiece](#developers-centerpiece)).

```bash
commander::runcmd -c bashbone -a cmds
```

## Enclosed scripts
[&#x25B2; back to top](#bashbone)

Bashbone is shipped with a couple of scripts to be used stand alone (experimental, when bashbone is not installed and activated) or being part of the biobash functions. They comprise volcano, PCA plot and heatmap generation as well as running `DESeq2` for differential expression or `WGCNA` for co-expression analysis, an annotate.pl script to add gene names to tables and postscript files that contain gene ids and scripts for normalizing expression values (TPM, FPKM) among others. They all have their own usage and can be listed via

```bash
bashbone -s

annotate.pl           canonicals.pl    deseq2.R                 dlgenome.sh
fftool.sh             fpkm.pl          genome2transcriptome.pl  heatmap.R
id2length.pl          mergefq.sh       mergexons.pl             mergexons.sh
pca.R                 pca_deseq.R      pileup2fastq.pl          revigo.R
rrbsMspIselection.sh  shufnsplitfq.sh  sra-dump.sh              survival.R
tpm.pl                vcfix.pl         vcfixuniq.pl             vizco.R
volcano.R             wgcna.R
```

### Retrieve SRA datasets
[&#x25B2; back to top](#bashbone)

Use the enclosed script to fetch sequencing data and meta information from NCBI SRA or EBI in a parallelized fashion, given SRA or GEO accession number or to convert local sra files.

```bash
sra-dump.sh -h

DESCRIPTION
sra-dump.sh retrieves fastq data from ncbi sra based on GSM, SRR or SRX accession numbers
  - support for parallel download instances
  - if installed, sra-toolkit via ncbi gov resource is priorized over retrieval via ebi uk mirror
  - ebi uk mirror can be defined as fallback upon fastq-dump errors

VERSION
0.4.0

REQUIREMENTS
Depends on chosen options
  - esearch (for non-SRR identifiers, from eutilities https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
  - pigz (for compression threads > 1)
  - wget or curl or fastq-dump (from stra-toolkit https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/)

SYNOPSIS
sra-dump.sh [OPTIONS] [GSM|SRR|SRX|SRP] [..]
sra-dump.sh -l [OPTIONS] [SRA] [..]

OPTIONS
-d [path] : download into this directory (default: "/u/people/kriege/workspace/bash/bashbone")
-p [num]  : number of maximum parallel download instances (default: 2)
-t [num]  : pigz compression threads per fastq-dump download instance (default: no pigz)
            HINT1: single threaded gzip compression may be a bottleneck regarding download speed
-m [path] : path to temporary directory (default: "/u/people/kriege/workspace/bash/bashbone/tmp.XXXXXXXXXX.sradump")
-l        : convert local sra files to compressed fastq files
-s        : show received meta information for given accession numbers and exit
-r        : do not fetch meta information for given SRR* accession numbers
-o [file] : additionally write information for given accession numbers to file
-a [list] : fetch all, biological and technical reads, and reassign mate ids according to a comma separated list e.g. 1,3,2
-w        : unless -a, download from ebi uk mirror utilizing wget
            HINT1: outperformes fastq-dump but uses less stable connections which may requires additional runs
            HINT2: use -p 1 in a final run to ensure all files were downloaded correctly
-c        : unless -a, download from ebi uk mirror utilizing curl
            HINT1: experimental!
-f        : unless -a, use ebi uk mirror as fallback upon fastq-dump failures. requires -w or -c
-z        : unless -a, download sra file from amazon aws utilizing awscli
            HINT1: uses -p 1 due to download of multiple chunks
            HINT2: use -l in a final run to convert sra files

EXAMPLES
sra-dump.sh -p 2 -t 4 -m /tmp GSM1446883 SRR1528586 SRX663213
sra-dump.sh -l -p 2 -t 4 -m /tmp /path/to/SRR1528586.sra /path/to/SRR1528588.sra
```

### Retrieve genomes
[&#x25B2; back to top](#bashbone)

Use the enclosed script to fetch human hg19/hg38 or mouse mm10/mm11 genomes, gene and ontology annotations, orthologs, dbSNP, MSigDB, Enrichr data. The Plug-n-play CTAT genome resource, made for gene fusion detection and shipped with STAR index, can be also selected optionally.

```bash
dlgenome.sh -h

DESCRIPTION
dlgenome.sh downloads most recent human or mouse genome and annotation including gene ontology plus orthologs and optionally dbSNP

VERSION
0.7.0

SYNOPSIS
dlgenome.sh -r [hg19|hg38|mm10] -g -a -d -s -n

INPUT OPTIONS
-h | --help               : prints this message
-t | --threads [value]    : threads - predicted default: 32
-o | --out [path]         : output directory - default: /u/people/kriege/workspace/bash/bashbone
-r | --reference [string] : choose hg19 hg38 mm10 mm11
-g | --genome             : download Ensembl genome
-c | --ctat               : switch to CTAT genome and indices (~30GB)
-a | --annotation         : download Ensembl gtf
-d | --descriptions       : download Ensembl gene description and Ensembl ontology information (requires R in PATH)
-m | --msigdb             : download Ensembl gene description, MSigDB gene ontology information and other collections (requires R in PATH)
-e | --enrichr            : download Ensembl gene description and Enrichr gene ontology and other collections (requires R in PATH)
                          : NOTE: for mouse data, genes will be substituted by Ensembl orthologs if possible
-s | --dbsnp              : download Ensembl dbSNP
-n | --ncbi               : switch to NCBI dbSNP
```

The genome, using annotation information, can be converted into a STAR transcritome output and salmon compatible transcriptome, transcript-genome or transcript-chromosome.

```bash
genome2transcriptome.pl -h

DESCRIPTION

Convert a genome and its annotation into an gtf input ordered
- transcriptome - i.e. each fasta sequence is a single transcript in 5' to 3' orientation
- transcriptgenome - i.e. each fasta sequence is a single transcript in genomic orientation
- transcriptchromosome - i.e. a fake chromosome fasta sequence containing all transcripts in genomic orientation separated by 1000 Ns
(c) Konstantin Riege - konstantin{.}riege{a}leibniz-fli{.}de

EXAMPLE

a chromosomal input gtf with

pos   i                    j
.....[====....====..========]..   transcript and exons
        --    ----  ----          and CDS (optional)

will be converted into

 1                j-i
[==================]   transcript and gene
 ====|====|========    exons
   --|----|----        CDS (if CDS in input)
 uu            uuuu    UTR (if CDS in input, 3'UTR includes STOP)

PARAMETER

-h | --help      print this
-v | --verbose   report progress
-t | --tmpdir    path to temporary directory (default: $TMPDIR or /tmp)
-f | --fasta     path to fasta input file
-g | --gtf       path to gtf input file
-o | --outdir    path to output directory

REQUIREMENTS

gtf with features:
  transcript
  exon (at least one per transcript, even for ncRNAs)
  CDS (optional)
and info fields:
  gene_id
  transcript_id

RESULTS
output gtf @ outdir/transcriptome.gtf outdir/transcriptgenome.gtf outdir/transcriptchromosome.gtf
output fasta @ outdir/transcriptome.fa outdir/transcriptgenome.fa outdir/transcriptchromosome.fa

```

### Merge/collate features
[&#x25B2; back to top](#bashbone)

Use the enclosed script to merge e.g. exons of multiple transcripts of a gene with optional offsets for elongation or shrinkage.

```bash
mergexons.sh -h

usage:
mergexons.sh <gtf> <fasta|fai|bam|sam> [<exon|gene> [+|-]<offset>]

version:
0.3.0

description:
- given a gtf file, this script merges exons on transcript level and prints exons, introns, genes, and intergenic regions in gtf format to stdout
- a fasta or fasta index (fai) or any matching bam/sam file is used to infer chromosome sizes
- optionally an offset first elongates gene or exon features up (-) and/or downstream (+) by optional [+|-] offset prefix

output:
- merged exons of multiple transcripts per gene
- introns according to merged exons
- genes according to merged exons
- intergenic regions according to merged genes

requirements:
- input file paths must not include white spaces or other special characters!
- needs mergexons.pl script next to this script
- needs bedtools to be in PATH (check: ok)
- needs samtools to be in PATH in case of sam/bam input (check: ok)
```

### Generate pseudo replicates
[&#x25B2; back to top](#bashbone)

To shuffle and split NGS raw data in fastq format into two pseudo-replicates, use

```bash
shufnsplitfq.sh -h

DESCRIPTION
shufnsplitfq.sh shuffle and split single or paired fastq(.gz|.bz2) files into two pseudo-replicate fastq(.gz) files

VERSION
0.1.0

SYNOPSIS
shufnsplitfq.sh -1 <fastq> [-2 <fastq>] -o <prefix> [-p <tmpdir> -t <threads> -z]

OPTIONS
-h          | this help
-1 <path>   | input SE or first mate fastq(.gz|.bz2)
-2 <path>   | input mate pair fastq(.gz|.bz2)
-z          | compress output using pigz with [-t] threads (fallback: gzip)
-t <value>  | compression threads (default: 32)
-o <path>   | prefix of output fastq, extended by [1|2].(R1.|R2.)fastq(.gz)
-p <path>   | temp directory (default: /tmp)
```

### Index and random access flat files
[&#x25B2; back to top](#bashbone)

The following script is designed to index flat files for random access by line number

```bash
DESCRIPTION
fftool.sh flat file index or random access

VERSION
0.1.0

SYNOPSIS
fftool.sh [-i] [-f <file>] [-o <file>]

OPTIONS
-h           | this help
-i           | index input file (see -f) or from stdin (requires -o)
-f <infile>  | path to file to be indexed (see -i) or accessed (requires -f). default: stdin (requires -o)
-o <outfile> | if input comes from stdin, path to output file (mutual exclusive to -f)
-r <range>   | for random data access. format: <#|inf>L@<#>L i.e. get a certain number of lines, all until EOF by inf keyword respectively, after skipping a certain number of lines.
             | example: 5L@20L means get 5 lines after skipping 20 lines
```

## Sample info file
[&#x25B2; back to top](#bashbone)

In order to perform desired comparative tasks, some functions require a sample design info file.

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

- At least 4 tab-separated columns (`<name>`, `<main-factor>`, `NA`, `<replicate>`)
- Optionally, additional factors
- First column needs to consist of unique prefixes of input fastq basenames which can be expand to full file names

|        |     |     |     |        |
| ---    | --- | --- | --- | ---    |
| wt1    | wt  | NA  | N1  | female |
| wt2    | wt  | NA  | N2  | male   |
| trA_1  | A   | NA  | N1  | female |
| trA_2  | A   | NA  | N2  | male   |
| trB.n1 | B   | NA  | N1  | female |
| trB.n2 | B   | NA  | N2  | male   |


## Adapter sequences
[&#x25B2; back to top](#bashbone)

Adapter sequences listed below will be tested by FastQC, extracted from the reports and stored as arrays. In case of paired-end data, unknown adapter sequences will be extracted from mate overlaps utilizing BBMap.

```bash
source <path/of/installation>/latest/bashbone/activate.sh
declare -a fastq_R1=(<path/to/file> [..<path/to/file>])
declare -a fastq_R2=(<path/to/file> [..<path/to/file>])
declare -a adapter_R1 adapter_R2
preprocess::fastqc -t <threads> -o <outdir> -1 fastq_R1 [-2 fastq_R2] -a1 adapter_R1 [-a2 adapter_R2]

```

Further adapter sequences can be found in the Illumina Adapter Sequences Document (<https://www.illumina.com/search.html?q=Illumina Adapter Sequences Document>) or Illumina Adapter Sequences HTML (<https://support-docs.illumina.com/SHARE/adapter-sequences.htm>) and the resource of Trimmomatic (<https://github.com/usadellab/Trimmomatic/tree/main/adapters>), FastQC respectively (<https://github.com/s-andrews/FastQC/blob/master/Configuration>).

The following excerpt is independent of the indexing type, i.e. single, unique dual (UD) or combinatorial dual (CD).

Nextera (Transposase Sequence), TruSight, AmpliSeq, stranded total/mRNA Prep, Ribo-Zero Plus: CTGTCTCTTATACACATCT

TruSeq (Universal) Adapter with A prefix due to 3' primer A-tailing : AGATCGGAAGAGC

TruSeq full length DNA & RNA R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq full length DNA MethC R1: AGATCGGAAGAGCACACGTCTGAAC R2: AGATCGGAAGAGCGTCGTGTAGGGA

TruSeq Small RNA 3': TGGAATTCTCGGGTGCCAAGG

Ovation Methyl-Seq R1: AGATCGGAAGAGC R2: AAATCAAAAAAAC

## Pipeline example
[&#x25B2; back to top](#bashbone)

Tiny example pipeline to perform gene fusion detection and differential -and coexpression expression analyses plus gene ontology enrichment. Check out further pre-compiled pipelines for peak calling from *IP-Seq experiments and differential expression- and methylation analysis from RNA-Seq data, RRBS/WGBS respectively (rippchen) or for multiple variant calling options from Exome-Seq/WGS/RNA-Seq data including GATK best-practices in an optimized, parallelized fashion (muvac).

```bash
source <path/of/installation>/latest/bashbone/activate.sh -c true
genome=<path/to/fasta>
genomeidx=<path/to/segemehl.idx>
gtf=<path/to/gtf>
declare -a fastqs=(<path/to/fastq/files> [..<path/to/fastq/files])) # can be compressed
threads=16
memory=64000

declare -a adapter
preprocess::fastqc -t $threads -M $memory -o results/qualities/raw -a adapter -1 fastqs
preprocess::trimmomatic -t $threads -M $memory -o results/trimmed -1 fastqs
preprocess::cutadapt -t $threads -o results/adapterclipped -a adapter -1 fastqs

preprocess::rcorrector -t $threads -o results/corrected -1 fastqs
preprocess::sortmerna -t $threads -o results/rrnafiltered -1 fastqs

fragmentsize=200
fusions::arriba -t $threads -g $genome -a $gtf -o results/fusions -f $fragmentsize -1 fastqs
fusions::starfusion -t $threads -g $genome -g $gtf -o results/fusions -1 fastqs

declare -a mapped
alignment::segemehl -t $threads -g $genome -x $genomeidx -o results/mapped -1 fastqs -r mapped

alignment::postprocess -j uniqify -t $threads -o results/mapped -r mapped
alignment::postprocess -r sort -t $threads -o results/mapped -r mapped

declare -A strandness
alignment::inferstrandness -t $threads -g $gtf -r mapped -x strandness
quantify::featurecounts -t $threads -M $memory -g $gtf -o results/counted -r mapped -x strandness
quantify::tpm -t $threads -g $gtf -o results/counted -r mapped

declare -a comparisons=(<path/to/sample-info/files> [..<path/to/sample-info/files>])
expression::diego -t $threads -g $gtf -c comparisons -i results/counted -o results/diffexonjunctions -r mapped -x strandness
expression::deseq -t $threads -g $gtf -c comparisons -i results/counted -o results/diffgenes -r mapped
expression::join_deseq -t $threads -g $gtf -c comparisons -i results/counted -j results/diffgenes -o results/counted
cluster::coexpression -t $threads -M $memory -g $gtf -i results/counted -o results/coexpressed -r mapped

go=<path/to/gtf.go>
enrichment::go -t $threads -r mapper -c comparisons -l coexpressions -g $go -i results/deseq
```

# Third-party software
[&#x25B2; back to top](#bashbone)

| Tool | Source | DOI |
| ---  | ---    | --- |
| Arriba        | <https://github.com/suhrig/arriba/>                                 | NA |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BBTools       | <https://jgi.doe.gov/data-and-tools/software-tools/bbtools>         | 10.1371/journal.pone.0185056 |
| BWA           | <https://github.com/lh3/bwa>                                        | 10.1093/bioinformatics/btp324 |
| BWA-mem2      | <https://github.com/bwa-mem2/bwa-mem2>                              | 10.1109/IPDPS.2019.00041 |
| BWA-meth      | <https://github.com/brentp/bwa-meth>                                | arXiv:1401.1129 |
| BCFtools      | <http://www.htslib.org/doc/bcftools.html>                           | 10.1093/bioinformatics/btr509 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| gztool        | <https://github.com/circulosmeos/gztool>                            | NA |
| cgpBigWig     | <https://github.com/cancerit/cgpBigWig>                             | NA |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler>        | 10.1089/omi.2011.0118 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DANPOS3       | <https://github.com/sklasfeld/DANPOS3>                              | 10.1101/gr.142067.112 |
| deepTools2    | <https://deeptools.readthedocs.io/en/latest/index.html>             | 10.1093/nar/gkw257 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 <br> 10.1093/biostatistics/kxw041|
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| dupsifter     | <https://github.com/huishenlab/dupsifter>                           | 10.1093/bioinformatics/btad729 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656 |
| fgbio         | <http://fulcrumgenomics.github.io/fgbio/>                           | NA |
| freebayes     | <https://github.com/ekg/freebayes>                                  | arXiv:1207.3907 |
| GATK4         | <https://github.com/broadinstitute/gatk>                            | 10.1101/gr.107524.110 <br> 10.1038/ng.806 |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GNU Parallel  | <https://www.gnu.org/software/parallel/>                            | 10.5281/zenodo.1146014 |
| GoPeaks       | <https://github.com/maxsonBraunLab/gopeaks>                         | 10.1186/s13059-022-02707-w |
| GoSemSim      | <http://bioconductor.org/packages/release/bioc/html/GOSemSim.html>  | 10.1093/bioinformatics/btq064 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HOMER         | <http://homer.ucsd.edu/homer>                                       | 10.1016/j.molcel.2010.05.004 |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/nboley/idr>                                     | 10.1214/11-AOAS466 |
| IGV           | <http://software.broadinstitute.org/software/igv>                   | 10.1038/nbt.1754 |
| Intervene     | <https://github.com/asntech/intervene>                              | 10.1186/s12859-017-1708-7 |
| kent/UCSC utilities | <https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads> | 10.1093/bioinformatics/btq351 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| KMC           | <https://github.com/refresh-bio/KMC>                                | 10.1186/1471-2105-14-160 |
| m6aViewer     | <http://dna2.leeds.ac.uk/m6a/>                                      | 10.1261/rna.058206.116 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| MethylDackel  | <https://github.com/dpryan79/MethylDackel>                          | NA |
| metilene      | <https://www.bioinf.uni-leipzig.de/Software/metilene/>              | 10.1101/gr.196394.115 |
| moose2        | <http://grabherr.github.io/moose2/>                                 | 10.1186/s13040-017-0150-8 |
| PEAKachu      | <https://github.com/tbischler/PEAKachu>                             | NA |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Platypus      | <https://rahmanteamdevelopment.github.io/Platypus>                  | 10.1038/ng.3036 |
| pugz          | <https://github.com/Piezoid/pugz>                                   | 10.1109/IPDPSW.2019.00042 |
| rapidgzip     | <https://github.com/mxmlnkn/rapidgzip>                              | 10.1145/3588195.3592992 |
| RAxML         | <https://cme.h-its.org/exelixis/web/software/raxml/index.html>      | 10.1093/bioinformatics/btl446 |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| RSeQC         | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| RRVGO         | <https://ssayols.github.io/rrvgo>                                   | 10.17912/micropub.biology.000811 |
| Salmon        | <https://combine-lab.github.io/salmon/>                             | 10.1038/nmeth.4197 |
| SalmonTE      | <https://github.com/hyunhwan-jeong/SalmonTE>                        | 10.1142/9789813235533_0016 |
| SAMtools      | <http://www.htslib.org/doc/samtools.html>                           | 10.1093/bioinformatics/btp352 |
| SEACR         | <https://github.com/FredHutch/SEACR>                                | 10.1186/s13072-019-0287-4 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SnpEff        | <https://pcingola.github.io/SnpEff>                                 | 10.4161/fly.19695 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| STAR-Fusion   | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>                   | 10.1101/120295 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| UMI-tools     | <https://github.com/CGATOxford/UMI-tools>                           | 10.1101/gr.209601.116 |
| VarDict       | <https://github.com/AstraZeneca-NGS/VarDict>                        | 10.1093/nar/gkw227 |
| VarScan       | <http://dkoboldt.github.io/varscan>                                 | 10.1101/gr.129684.111 |
| vcflib        | <https://github.com/vcflib/vcflib>                                  | 10.1371/journal.pcbi.1009123 |
| Vt            | <https://genome.sph.umich.edu/wiki/Vt>                              | 10.1093/bioinformatics/btv112 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |
| WiggleTools   | <https://github.com/Ensembl/WiggleTools>                            | 10.1093/bioinformatics/btt737 |

# Supplementary information
[&#x25B2; back to top](#bashbone)

In some rare cases a glibc pthreads bug (<https://sourceware.org/bugzilla/show_bug.cgi?id=23275>) may cause pigz failures (`internal threads error`) and premature termination of tools leveraging on it e.g. Cutadapt and pigz. One can circumvent this by making use of a newer glibc version, an alternative pthreads library e.g. compiled without lock elision via `LD_PRELOAD` or by using the statically compiled pigz binary shipped with bashbone.

```bash
LD_PRELOAD=</path/to/no-elision/libpthread.so.0> <command>
```

# Closing remarks
[&#x25B2; back to top](#bashbone)

Bashbone is a continuously developed library and actively used in my daily work. As a single developer it may take me a while to fix errors and issues. Feature requests cannot be handled so far, but I am happy to receive pull request.
