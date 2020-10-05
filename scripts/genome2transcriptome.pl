#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;
use Bio::Index::Fasta;
use Getopt::Long;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path remove_tree);
use File::Temp;
use Try::Tiny;
use sigtrap qw(handler cleanup END);

my @tmp;
sub cleanup {
	for (@tmp){
		unlink $_ for glob "$_*";
	}
}

END {
	&cleanup();
}

sub mydie {
	say STDERR $_ for @_;
	&cleanup();
	exit 1;
}

sub usage {
print <<EOF;
DESCRIPTION
Convert a genome and its annotation into
- a transcriptome - i.e. each faste entry is a single transcript and
- a transcriptgenome - i.e. single fake chromosome containing all transcripts N-separated
(c) Konstantin Riege - konstantin{.}riege{a}leibniz-fli{.}de

an chromosomal input gtf with

pos  i                            j
.....[=======....========..=======]..   transcript and exons
                 --------  ----         and optional CDS

will be converted into

1                     j-i
[======================]   transcript as "gene"
=======|========|=======   exons
          ------|----      CDS
                     s     STOP (actually not printed)
uuuuuuuuuu            uu   UTR

PARAMETER

-h | --help      does what it says
-v | --verbose   say something
-f | --fasta     path to fasta input file
-g | --gtf       path to gtf input file
-o | --outdir    path to output directory

REQUIREMENTS
gtf with features:
  transcript
  exon (at least one per transcript, even for ncRNAs)
  CDS (optional)
and info fields:
  transcript_id
  exon_number

RESULTS
output gtf @ outdir/transcriptome.gtf outdir/transcriptgenome.gtf
output fasta @ outdir/transcriptome.fa outdir/transcriptgenome.fa

EOF
	exit 0;
}

&usage if $#ARGV == -1;
my ($verbose, $fasta, $gtf, $outdir);
(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{&usage},
	'v|verbose' => \$verbose,
	'f|fasta=s' => \$fasta,
	'g|gtf=s' => \$gtf,
	'o|odir=s' => \$outdir,
) or &usage;
&usage unless $fasta && $gtf && $outdir;

try {
	push @tmp,File::Temp->new(DIR => "/dev/shm/", SUFFIX => ".faidx")->filename;
	my $fa=Bio::Index::Fasta->new(-filename => $tmp[-1], -write_flag => 1 , -verbose => -1);
	$fa->id_parser(sub{
		my ($h)=@_;
		$h=~/^>\s*(\S+)/;
		return $1;
	});
	$fa->make_index($fasta);

	push @tmp,File::Temp->new(DIR => "/dev/shm/", SUFFIX => ".ingtf")->filename;
	my $cmd="sort -k3,3r -k1,1 -k4,4n -k5,5n '$gtf' > $tmp[-1]";
	system($cmd) == 0 or &mydie($?);
	$gtf=$tmp[-1];

	push @tmp,File::Temp->new(DIR => "/dev/shm/", SUFFIX => ".outgtf")->filename;
	my $outgtf=$tmp[-1];


	make_path($outdir);
	my (@F, $tid, %t, $eno, %eorigsta, %tseq, $seq, $subseq, %enewsta, $sta, %csta, %csto);
	my $chr="";
	open GTF,"<$gtf" or &mydie($!);
	open TGTF,">$outgtf" or &mydie($!);
	while(<GTF>){
		chomp;
		@F = split/\t/,$_;
		if($F[2] eq "transcript"){ #due to sort -k3,3r, handle transcript features first
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			$F[-1]=~s/gene_id "[^"]+"/gene_id "$tid"/;

			$t{$tid}=join"\t",@F;
			$tseq{$tid}=""; # store id and initialize sequence
		}
		if($F[2] eq "exon"){
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			next unless exists $t{$tid}; # if valid gtf, this should not occure
			$F[-1]=~s/gene_id "[^"]+"/gene_id "$tid"/;

			$seq=$fa->fetch($F[0]) unless $chr eq $F[0]; # get chr seq obj
			$chr=$F[0];
			$subseq=$seq->subseq($F[3],$F[4]);
			$F[-1]=~/exon_number "([^"]+)/;
			$eno=$1;
			$eorigsta{"$tid$eno"}=$F[3]; # keep original position to recalculate related CDS
			$F[3]=length($tseq{$tid})+1; # set new start position
			$tseq{$tid}.=$subseq; # elongate transcript sequence
			$F[4]=$F[3]+length($subseq)-1; # set new stop position
			$enewsta{"$tid$eno"}=$F[3]; # keep new position to recalculate related CDS
			$F[0]=$tid;
			say TGTF join"\t",@F;
			say join"\t",@F if $verbose;
		}
		if($F[2] eq "CDS"){
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			$F[-1]=~/exon_number "([^"]+)/;
			$eno=$1;
			next unless exists $enewsta{"$tid$eno"};
			$F[-1]=~s/gene_id "[^"]+"/gene_id "$tid"/;

			$sta=$enewsta{"$tid$eno"}+$F[3]-$eorigsta{"$tid$eno"}; # calculate new start position
			$F[4]=$sta+$F[4]-$F[3]; # set new stop position
			$F[3]=$sta; # set new start position
			$F[0]=$tid;
			say TGTF join"\t",@F;
			say join"\t",@F if $verbose;
			$csta{$tid}=$F[3] unless exists $csta{$tid}; # due to sort -k4,4 keep very first CDS start position for UTR calculation
			$csto{$tid}=$F[4]; # update until last CDS stop position for UTR calculation
		}
	}
	close GTF;

	open TFA,">".catfile($outdir,"transcriptome.fa") or &mydie($!);
	open GFA,">".catfile($outdir,"transcriptgenome.fa") or &mydie($!);
	say GFA ">chr1";
	say GFA $_ for unpack "(A100)*","N"x1000;
	my ($s, $fill);
	my $overhang=0;
	for $tid (sort {$a cmp $b} keys %t){
		@F=split/\t/,$t{$tid};
		$F[0]=$tid;
		$F[2]="gene"; # treat transcript as gene in new annotation
		$F[3]=1; # set start posotion
		$F[4]=length($tseq{$tid}); # set stop position
		say TGTF join"\t",@F;
		say join"\t",@F if $verbose;
		if(exists $csto{$tid} && $csto{$tid}<$F[4]){ # if right hand CDS present, calculate UTR
			$F[2]="5primeUTR";
			$F[3]=$csto{$tid}+1;
			if($F[6] eq "+"){
				$F[2]="3primeUTR";
				$F[3]+=3 if $F[3]+3<$F[4]; # move UTR behind stop codon
			}
			say TGTF join"\t",@F;
			say join"\t",@F if $verbose;
		}
		if(exists $csta{$tid} && $csta{$tid}>1){ # if left hand CDS present, calculate UTR
			$F[2]="5primeUTR";
			$F[3]=1;
			$F[4]=$csta{$tid}-1;
			if($F[6] eq "-"){
				$F[2]="3primeUTR";
				$F[4]-=3 if $F[4]-4>1; # move UTR behind stop codon
			}
			say TGTF join"\t",@F;
			say join"\t",@F if $verbose;
		}

		$s=$tseq{$tid};
		say TFA ">$tid";
		say TFA $_ for unpack "(A100)*",$s;

		$s="N"x$overhang.$s;
		$overhang=length($s)%1000;
		$fill=1000-$overhang;
		say GFA $_ for unpack "(A100)*",$s."N"x$fill;
	}
	close TGTF;
	close TFA;
	say GFA $_ for unpack "(A100)*","N"x$overhang;
	close GFA;

	system("sort -k1,1 -k4,4n -k5,5n $outgtf > ".catfile($outdir,"transcriptome.gtf")) == 0 or &mydie($?);
	open GTF,"<".catfile($outdir,"transcriptome.gtf") or die &mydie($!);
	open GGTF,">".catfile($outdir,"transcriptgenome.gtf") or die &mydie($!);
	my ($l, $i) = (0, 0);
	$chr="";
	while(<GTF>){
		chomp;
		@F = split/\t/,$_;
		$l=$F[4] if $F[2] eq "gene";
		$i+=1000+$l unless $chr eq $F[0];
		$chr=$F[0];
		$F[0]="chr1";
		$F[3] += $i;
		$F[4] += $i;
		say GGTF join"\t",@F;
		say join"\t",@F if $verbose;
	}
	close GTF;
	close GGTF;
} catch {
	&mydie($_);
};

exit 0;
