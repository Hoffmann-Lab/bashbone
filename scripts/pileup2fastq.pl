#!/usr/bin/env perl
use v5.10;
use Getopt::Long;

sub usage {
say <<'EOF';
DESCRIPTION
Converts SAM/BAM files - i.e. primary alignments of mapped and coordiante sorted DNA-seq derived datasets (e.g. WGBS, WES) into reference based FastQ files.
Following the pre- and post-processing steps, this script preserves splitread/chimeric information from soft clipped sites and their supplementary alignments.
Please note, that insertions and deletions will be ignored.


OPTIONS
-h | --help       : this message
-l | --legacy     : enable (slower) legacy mode if samtools --version is < 1.14
-m | --mate [1|2] : append proper read mate information to read id (for paired-end data)


SYNOPSIS (requires samtools >= v1.14. See below for legacy synopsis)
# for single-end data:
samtools mpileup -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP-5 --output-QNAME --no-output-ins --no-output-ins --no-output-del --no-output-del -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -l > [fastq.tsv]

# for paired-end data - mate 1:
samtools mpileup --rf 64 -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP-5 --output-QNAME --no-output-ins --no-output-ins --no-output-del --no-output-del -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -m 1 > [fastq.R1.tsv]
# for paired-end data - mate 2:
samtools mpileup --ff 64 -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP-5 --output-QNAME --no-output-ins --no-output-ins --no-output-del --no-output-del -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -m 2 > [fastq.R2.tsv]


SYNOPSIS (LEGACY MODE)
# for single-end data:
samtools mpileup -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP --output-QNAME -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -l > [fastq.tsv]

# for paired-end data - mate 1:
samtools mpileup --rf 64 -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP --output-QNAME -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -l -m 1 > [fastq.R1.tsv]
# for paired-end data - mate 2:
samtools mpileup --ff 64 -A -q 0 -Q 0 -x -B -C 0 -d 999999 --output-BP --output-QNAME -f [sup_genome.fa] [prepped.bam] | pileup2fastq.pl -l -m 2 > [fastq.R2.tsv]


PREPROCESSING
# duplicate genome and rename chromosomes to pileup supplementary alignments first
cat <(awk '{if(/^>/){print ">sup_"substr($0,2)}else{print}}' [genome.fa]) [genome.fa] > [sup_genome.fa]
samtools faidx [sup_genome.fa]
samtools view -H [wgbs.bam] > [wgbs.header]
n=$(grep -nF SN: [wgbs.header] | head -1 | cut -d ':' -f 1)
cat <(awk -v n=$n 'NR<n' [wgbs.header]) <(sed -nE 's/SN:(.+)/SN:sup_\1/p' [wgbs.header]) <(awk -v n=$n 'NR>=n' [wgbs.header]) > [reheader.header]

# for single-end data:
cat [reheader.header] <(samtools view -F 4 -f 2048 [wgbs.bam] | awk -F '\t' -v OFS='\t' '{$3="sup_"$3; print}') | samtools view -@ [threads] -b > [sup.bam]
cat [reheader.header] <(samtools view -F 4 -F 2048 -F 256 [wgbs.bam] | perl -F'\t' -lanE 'if($F[5]=~s/(\d+)S$/$1M/){substr($F[9],-$1,$1)=~tr/ACGTacgt/VMWDvmwd/} if($F[5]=~s/^(\d+)S/$1M/){substr($F[9],0,$1)=~tr/ACGTacgt/VMWDvmwd/} say join"\t",@F') | samtools view -@ [threads] -b > [primary.bam]

# for paired-end data:
cat [reheader.header] <(samtools view -f 1 -F 4 -F 8 -f 2048 [wgbs.bam] | awk -F '\t' -v OFS='\t' '{$3="sup_"$3; print}') | samtools view -@ [threads] -b > [sup.bam]
cat [reheader.header] <(samtools view -f 1 -F 4 -F 8 -F 2048 -F 256 test.bam | perl -F'\t' -lanE 'if($F[5]=~s/(\d+)S$/$1M/){substr($F[9],-$1,$1)=~tr/ACGTacgt/VMWDvmwd/} if($F[5]=~s/^(\d+)S/$1M/){substr($F[9],0,$1)=~tr/ACGTacgt/VMWDvmwd/} say join"\t",@F') | samtools view -@ [threads] -b > [primary.bam]

# merge and index
samtools merge -@ $threads -f -c -p [prepped.bam] [sup.bam] [primary.bam]
samtools index -@ $threads [prepped.bam] [prepped.bai]


PREPROCESSING (LEGACY MODE)
# duplicate genome and rename chromosomes to pileup supplementary alignments first
cat <(awk '{if(/^>/){print ">sup_"substr($0,2)}else{print}}' [genome.fa]) [genome.fa] > [sup_genome.fa]
samtools faidx [sup_genome.fa]
samtools view -H [wgbs.bam] > [wgbs.header]
n=$(grep -nF SN: [wgbs.header] | head -1 | cut -d ':' -f 1)
cat <(awk -v n=$n 'NR<n' [wgbs.header]) <(sed -nE 's/SN:(.+)/SN:sup_\1/p' [wgbs.header]) <(awk -v n=$n 'NR>=n' [wgbs.header]) > [reheader.header]

# for single-end data:
cat [reheader.header] <(samtools view -F 4 -f 2048 [wgbs.bam] | awk -F '\t' -v OFS='\t' '{$1=$1":"length($10); $3="sup_"$3; print}') | samtools view -@ [threads] -b > [sup.bam]
cat [reheader.header] <(samtools view -F 4 -F 2048 -F 256 [wgbs.bam] | perl -F'\t' -lanE '$F[0].=":".length($F[9]); if($F[5]=~s/(\d+)S$/$1M/){substr($F[9],-$1,$1)=~tr/ACGTacgt/VMWDvmwd/} if($F[5]=~s/^(\d+)S/$1M/){substr($F[9],0,$1)=~tr/ACGTacgt/VMWDvmwd/} say join"\t",@F') | samtools view -@ [threads] -b > [primary.bam]

# for paired-end data:
cat [reheader.header] <(samtools view -f 1 -F 4 -F 8 -f 2048 [wgbs.bam] | awk -F '\t' -v OFS='\t' '{$1=$1":"length($10); $3="sup_"$3; print}') | samtools view -@ [threads] -b > [sup.bam]
cat [reheader.header] <(samtools view -f 1 -F 4 -F 8 -F 2048 -F 256 test.bam | perl -F'\t' -lanE '$F[0].=":".length($F[9]); if($F[5]=~s/(\d+)S$/$1M/){substr($F[9],-$1,$1)=~tr/ACGTacgt/VMWDvmwd/} if($F[5]=~s/^(\d+)S/$1M/){substr($F[9],0,$1)=~tr/ACGTacgt/VMWDvmwd/} say join"\t",@F') | samtools view -@ [threads] -b > [primary.bam]

# merge and index
samtools merge -@ $threads -f -c -p [prepped.bam] [sup.bam] [primary.bam]
samtools index -@ $threads [prepped.bam] [prepped.bai]


POSTPROCESSING
# for single-end data:
tr \$'\t' \$'\n' [fastq.tsv] | pigz -c -p [threads] > [fastq.gz]

# for paired-end data:
sort -k1,1 --parallel=[threads] -S [memory] -T /dev/shm  [fastq.R1.tsv] | tr \$'\t' \$'\n' | pigz -c -p [threads] > [fastq.R1.gz]
sort -k1,1 --parallel=[threads] -S [memory] -T /dev/shm  [fastq.R2.tsv] | tr \$'\t' \$'\n' | pigz -c -p [threads] > [fastq.R2.gz]


PRALLELIZATION
awk -v OFS='\t' '/^sup_/{print $1,0,$2}' samtools faidx [sup_genome.fa.fai] > [sup_regions.bed]
# for each chromosome mpileup commands can be executed in parallel instances
samtools mpileup [..] -l [sup_regions.bed] -r [chr] [prepped.bam]
# chromosomes could be further splitted into smaller chunks by -r [chr:start-end]. for sake of read integrity choose start and end postions with zero coverage


REFERENCES
(c) Konstantin Riege
EOF
	exit 0;
}

#&usage if $#ARGV == -1;
$revcom=0;
$r="";
$l="";
(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{&usage},
	'm|mate=i' => \$r,
	'l|legacy' => \$l,
) or &usage;
$r=" $r" if $r;

if($l){
	while(<>){
		chomp;
		@F=split/\t/;
		$ref=uc($F[2]);
		@bq=split//,$F[5];
		@pos=split/,/,$F[6];

		@id=split/,/,$F[7];
		@len=();
		for (@id){
			$_=~s/:(\d+)$//;
			push @len,$1;
		}

		@nfo=split//,$F[4];

		$j=0;
		for ($i=0; $i<=$#nfo; $i++){
			next if $nfo[$i]=~/[*#<>]/; # skip read del and splice junction i.e. cigar N
			if ($nfo[$i]=~/[-+]/){ # skip read ins
				$l="";
				while ($nfo[++$i]=~/(\d)/){
					$l.=$1;
				}
				$i+=$l-1;
				next;
			}
			$i+=2 if $nfo[$i] eq "^"; # skip read start

			unless(exists $read{$id[$j]}->[$pos[$j]]) { # supplementary always wins
				if ($nfo[$i]=~/[A-Z.]/){
					$read{$id[$j]}->[$pos[$j]] = $nfo[$i]=~/[ACGT.]/ ? $ref : $nfo[$i]=~tr/VMWD/ACGT/r;
				}else{
					# reverse read if necessary - do this here and not during write out, because supplementary alignments can have different orientations
					$pos[$j] = $len[$j]+1-$pos[$j];
					$read{$id[$j]}->[$pos[$j]] = $nfo[$i]=~/[acgt,]/ ? $ref=~tr/ACGT/TGCA/r : $nfo[$i]=~tr/vmwd/TGCA/r;
				}
				$qual{$id[$j]}->[$pos[$j]] = $bq[$j];
			}

			if($i+1<=$#nfo && $nfo[$i+1] eq '$'){ # return if read ends
				unless($F[0]=~/^sup_/){
					say join "\t",('@'."$id[$j]$r" , join("",$read{$id[$j]}->@*) , "+" , join("",$qual{$id[$j]}->@*)); # print if primary
					delete $read{$id[$j]}; # free memory
					delete $qual{$id[$j]};
				}
				$i++;
			}
			$j++;
		}
	}
} else {
	while(<>){
		chomp;
		@F=split/\t/;
		$ref=uc($F[2]);
		@nfo=split//,$F[4];
		@bq=split//,$F[5];
		@pos=split/,/,$F[7]; # true for samtools v1.14
		@id=split/,/,$F[6]; # true for samtools v1.14

		$j=0;
		for ($i=0; $i<=$#nfo; $i++){
			next if $nfo[$i]=~/[*#<>]/; # skip read del and splice junction i.e. cigar N
			$i+=2 if $nfo[$i] eq "^"; # skip read start

			unless(exists $read{$id[$j]}->[$pos[$j]]) { # supplementary always wins
				if ($nfo[$i]=~/[A-Z.]/){
					$read{$id[$j]}->[$pos[$j]] = $nfo[$i]=~/[ACGT.]/ ? $ref : $nfo[$i]=~tr/VMWD/ACGT/r;
				}else{
					# reverse read if necessary - do this here and not during write out, because supplementary alignments can have different orientations
					$read{$id[$j]}->[$pos[$j]] = $nfo[$i]=~/[acgt,]/ ? $ref=~tr/ACGT/TGCA/r : $nfo[$i]=~tr/vmwd/TGCA/r;
				}
				$qual{$id[$j]}->[$pos[$j]] = $bq[$j];
			}

			if($i+1<=$#nfo && $nfo[$i+1] eq '$'){ # return if read ends
				unless($F[0]=~/^sup_/){
					say join "\t",('@'."$id[$j]$r" , join("",$read{$id[$j]}->@*) , "+" , join("",$qual{$id[$j]}->@*)); # print if primary
					delete $read{$id[$j]}; # free memory
					delete $qual{$id[$j]};
				}
				$i++;
			}
			$j++;
		}
	}
}