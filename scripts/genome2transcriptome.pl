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
use sigtrap qw(handler mydie INT TERM KILL);

my $tmp;
sub mydie {
	say STDERR $_ for @_;
	unlink $tmp->filename if defined $tmp;
	exit 1;
}

sub usage {
print <<EOF;
DESCRIPTION
Convert a genome and its annotation into a transcriptome
(c) Konstantin Riege - konstantin{.}riege{a}leibniz-fli{.}de

i.e. an input gtf with

pos  i                            j
.....[=======....========..=======]..   transcript and exons
                 --------  ----         and optional CDS

will be converted into

1                     j-i
[======================]   transcript as "genes"
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
input gtf needs to be sorted this way: 
  sort -k3,3r -k1,1 -k4,4n -k5,5n
gtf features: 
  transcript
  exon (at least one per transcript, even for ncRNAs
  CDS (optional)
gtf info fields:
  transcript_id
  exon_number

RESULTS
unsorted gtf @ outdir/transcriptome.gtf
output fasta @ outdir/transcriptome.fa

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
	$tmp=File::Temp->new(DIR => "/dev/shm/", SUFFIX => ".faidx");
	my $fa=Bio::Index::Fasta->new(-filename => $tmp->filename, -write_flag => 1 , -verbose => -1);
	$fa->id_parser(sub{
		my ($h)=@_;
		$h=~/^>\s*(\S+)/;
		return $1;
	});
	$fa->make_index($fasta);

	make_path($outdir);
	my (@F, $tid, %t, $eno, %eorigsta, %tseq, $seq, $subseq, %enewsta, $sta, %csta, %csto);
	my $chr="";
	open GTF,"<$gtf" or &mydie($!);
	open OGTF,">".catfile($outdir,"transcriptome.gtf") or &mydie($!);
	while(<GTF>){
		chomp;
		@F = split/\t/,$_;
		if($F[2] eq "transcript"){ #due to sort -k3,3r, handle transcript features first
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			$t{$tid}=join"\t",@F;
			$tseq{$tid}=""; # store id and initialize sequence
		}
		if($F[2] eq "exon"){
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			next unless exists $t{$tid}; # if valid gtf, this should not occure
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
			say OGTF join"\t",@F;
			say join"\t",@F if $verbose;
		}
		if($F[2] eq "CDS"){
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			$F[-1]=~/exon_number "([^"]+)/;
			$eno=$1;
			next unless exists $enewsta{"$tid$eno"};
			$sta=$enewsta{"$tid$eno"}+$F[3]-$eorigsta{"$tid$eno"}; # calculate new start position
			$F[4]=$sta+$F[4]-$F[3]; # set new stop position
			$F[3]=$sta; # set new start position
			say OGTF join"\t",@F;
			say join"\t",@F if $verbose;
			$csta{$tid}=$F[3] unless exists $csta{$tid}; # due to sort -k4,4 keep very first CDS start position for UTR calculation
			$csto{$tid}=$F[4]; # update until last CDS stop position for UTR calculation
		}
	}
	close GTF;

	open OFA,">$outdir/transcriptome.fa" or &mydie($!);
	for $tid (sort {$a cmp $b} keys %t){
		@F=split/\t/,$t{$tid};
		$F[2]="gene"; # treat transcript as gene in new annotation
		$F[3]=1; # set start posotion
		$F[4]=length($tseq{$tid}); # set stop position
		say OGTF join"\t",@F;
		say join"\t",@F if $verbose;
		if(exists $csto{$tid} && $csto{$tid}<$F[4]){ # if right hand CDS present, calculate UTR
			$F[2]="5primeUTR";
			$F[3]=$csto{$tid}+1;
			if($F[6] eq "+"){
				$F[2]="3primeUTR";
				$F[3]+=3 if $F[3]+3<$F[4]; # move UTR behind stop codon
			}
			say OGTF join"\t",@F;
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
			say OGTF join"\t",@F;
			say join"\t",@F if $verbose;
		}
		say OFA ">$tid";
		say OFA $_ for unpack "(A100)*",$tseq{$tid};
	}

	close OGTF;
	close OFA;
} catch {
	&mydie($_);
};

exit 0;
