#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
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
	say STDERR "ERROR $_" for @_;
	&cleanup();
	exit 1;
}

sub usage {
print <<EOF;
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
-t | --tmpdir    path to temporary directory (default: \$TMPDIR or /tmp)
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

EOF
	exit 0;
}

&usage if $#ARGV == -1;
my ($verbose, $fasta, $gtf, $outdir, $tmpdir);
(Getopt::Long::Parser->new)->getoptions(
	'h|help' => sub{&usage},
	'v|verbose' => \$verbose,
	'f|fasta=s' => \$fasta,
	't|tmpdir=s' => \$tmpdir,
	'g|gtf=s' => \$gtf,
	'o|odir=s' => \$outdir,
) or &usage;
&usage unless $gtf && $outdir;
$tmpdir = defined $ENV{TMPDIR} ? $ENV{TMPDIR} : "/tmp" unless $tmpdir;

try {
	make_path($outdir);

	open TGTF,">".catfile($outdir,"transcriptome.gtf") or &mydie($!);
	open TGGTF,">".catfile($outdir,"transcriptgenome.gtf") or &mydie($!);
	open TCGTF,">".catfile($outdir,"transcriptchromosome.gtf") or &mydie($!);

	my $fa;
	if ($fasta){
		open TFA,">".catfile($outdir,"transcriptome.fa") or &mydie($!);
		open TGFA,">".catfile($outdir,"transcriptgenome.fa") or &mydie($!);
		open TCFA,">".catfile($outdir,"transcriptchromosome.fa") or &mydie($!);

		say STDERR "reading FASTA" if $verbose;
		push @tmp,File::Temp->new(DIR => $tmpdir, SUFFIX => ".faidx")->filename;
		$fa=Bio::Index::Fasta->new(-filename => $tmp[-1], -write_flag => 1 , -verbose => -1);
		$fa->id_parser(sub{
			my ($h)=@_;
			$h=~/^>\s*(\S+)/;
			return $1;
		});
		$fa->make_index($fasta);
	}

	say STDERR "reading GTF" if $verbose;
	my %gtf;
	my %t2g;
	my @tids;
	my $tid;
	my %ti;

	# sometimes exons dont have a cds. order matters!
	open GTF,"<$gtf" or &mydie($!);
	while(<GTF>){
		chomp;
		my @F=split/\t/;
		$F[-1]=~/transcript_id "([^"]+)/;
		$tid=$1;
		if($F[2] eq "transcript"){
			$F[-1]=~/transcript_id "([^"]+)/;
			$tid=$1;
			$F[-1]=~/gene_id "([^"]+)/;
			$t2g{$tid}=$1;
			push @tids, $tid;
			push $gtf{$tid}->{transcript}->@*,$_;
		}
	}
	seek GTF, 0, 0;
	while(<GTF>){
		chomp;
		my @F=split/\t/;
		$F[-1]=~/transcript_id "([^"]+)/;
		$tid=$1;
		if($F[2] eq "exon"){
			push $gtf{$tid}->{Gexon}->@*,$_;
			if($F[6] eq "+"){
				push $gtf{$tid}->{exon}->@*,$_;
			} else {
				unshift $gtf{$tid}->{exon}->@*,$_; # reverse order for transcriptome
			}
		}
	}
	seek GTF, 0, 0;
	while(<GTF>){
		chomp;
		my @C=split/\t/;
		$C[-1]=~/transcript_id "([^"]+)/;
		$tid=$1;
		if($C[2] eq "CDS"){
			unless(exists $gtf{$tid}->{CDS}){
				$gtf{$tid}->{CDS}->@* = ("") x ($gtf{$tid}->{Gexon}->$#*+1);
				$gtf{$tid}->{GCDS}->@* = ("") x ($gtf{$tid}->{Gexon}->$#*+1);
				$ti{$tid}=0;
			}
			# pair exons and CDS (without requiring exon_number tag), because not all exons have a CDS
			my @F=split/\t/,$gtf{$tid}->{Gexon}->[$ti{$tid}];
			while( !($C[3]>=$F[3] && $C[4]<=$F[4]) ){
				$ti{$tid}++;
				@F=split/\t/,$gtf{$tid}->{Gexon}->[$ti{$tid}];
			}

			$gtf{$tid}->{GCDS}->[$ti{$tid}]=$_;
			if($C[6] eq "+"){
				$gtf{$tid}->{CDS}->[$ti{$tid}]=$_;
			} else {
				$gtf{$tid}->{CDS}->[$gtf{$tid}->{Gexon}->$#*-$ti{$tid}]=$_;
			}
		}
	}
	close GTF;

	my $chr="";
	my $seq;
	my $rev;
	my $gsta=1000;
	my $gseq="N" x 1000;
	for $tid (@tids){
		####################################################
		# transcriptome
		my $tseq="";
		my $sta=0;
		my $sto=0;
		my @out=();
		my @F;

		for my $i (0..$gtf{$tid}->{exon}->$#*){
			@F=split/\t/,$gtf{$tid}->{exon}->[$i];

			if($fasta){
				unless ($chr eq $F[0]){
					say STDERR "working on $F[0]" if $verbose;
					$seq=$fa->fetch($F[0]);
					$rev=$seq->revcom();
				}
				$chr=$F[0];

				if($F[6] eq "+"){
					$tseq.=$seq->subseq($F[3],$F[4]);
				} else {
					$tseq.=$rev->subseq($seq->length()-$F[4]+1,$seq->length()-$F[3]+1);
				}
			}
			$sta=$sto+1;
			$sto=$sta+($F[4]-$F[3]);

			push @out,join("\t",$tid,$F[1],$F[2],$sta,$sto,$F[5],"+",$F[7],$F[8]." genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";");

			# say STDERR $tid;
			if (exists $gtf{$tid}->{CDS} && $gtf{$tid}->{CDS}->[$i]){
				my $csta=0;
				my $csto=0;
				my @C=split/\t/,$gtf{$tid}->{CDS}->[$i];
				if ($C[3] == $F[3] && $C[4] == $F[4]){
					$csta=$sta;
					$csto=$sto;
				}elsif($C[4] == $F[4]){
					# + ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# + -------CDS----CDSCDS----CDS-------
					#       ^ 5'UTR
					# - ----E3E3E3----E2E2E2----E1E1E1---- (numbers in read order)
					# - ----CDS-------CDSCDS-------CDS----
					#                           ^ 5'UTR
					if ($C[6] eq "+"){
						$csto=$sto;
						$csta=$csto-($C[4]-$C[3]);
						# 5'UTR
						push @out,join("\t",$tid,$C[1],"5primeUTR",$sta,$csta-1,$C[5],"+",$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
					} else {
						$csta=$sta;
						$csto=$csta+($C[4]-$C[3]);
						# 3'UTR
						push @out,join("\t",$tid,$C[1],"3primeUTR",$csto+1,$sto,$C[5],"+",$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
					}
				} else {
					# + ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# + -------CDS----CDSCDS----CDS-------
					#                              ^ 3'UTR
					# - ----E3E3E3----E2E2E2----E1E1E1---- (numbers in read order)
					# - ----CDS-------CDSCDS-------CDS----
					#          ^ 3'UTR
					if ($C[6] eq "+"){
						$csta=$sta;
						$csto=$csta+($C[4]-$C[3]);
						# 3'UTR
						push @out,join("\t",$tid,$C[1],"3primeUTR",$csto+1,$sto,$C[5],"+",$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
					} else {
						$csto=$sto;
						$csta=$csto-($C[4]-$C[3]);
						# 5'UTR
						push @out,join("\t",$tid,$C[1],"5primeUTR",$sta,$csta-1,$C[5],"+",$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
					}
				}
				push @out,join("\t",$tid,$C[1],$C[2],$csta,$csto,$C[5],"+",$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
			}
		}
		@F=split/\t/,$gtf{$tid}->{transcript}->[0];
		$F[-1].=" genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";";
		$F[0]=$tid;
		$F[3]=1;
		$F[4]=$sto;
		$F[6]="+";
		unshift @out,join("\t",@F);
		$F[2]="gene";
		unshift @out,join("\t",@F);
		say TGTF $_ for sort {my @a=split/\t/,$a; my @b=split/\t/,$b; $a[3] <=> $b[3] || $a[4] <=> $b[4]} @out;
		if($fasta){
			say TFA ">$tid $t2g{$tid}";
			say TFA $_ for unpack "(A100)*",$tseq;
		}

		####################################################
		# transcriptgenome
		$tseq="";
		$sta=0;
		$sto=0;
		@out=();
		for my $i (0..$gtf{$tid}->{Gexon}->$#*){
			@F=split/\t/,$gtf{$tid}->{Gexon}->[$i];
			$tseq.=$seq->subseq($F[3],$F[4]) if $fasta;
			$sta=$sto+1;
			$sto=$sta+($F[4]-$F[3]);

			push @out,join("\t",$tid,$F[1],$F[2],$sta,$sto,$F[5],$F[6],$F[7],$F[8]." genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";");

			# say STDERR $tid;
			if (exists $gtf{$tid}->{GCDS} && $gtf{$tid}->{GCDS}->[$i]){
				my $csta=0;
				my $csto=0;
				my @C=split/\t/,$gtf{$tid}->{GCDS}->[$i];
				if ($C[3] == $F[3] && $C[4] == $F[4]){
					$csta=$sta;
					$csto=$sto;
				}elsif($C[4] == $F[4]){
					# + ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# + -------CDS----CDSCDS----CDS-------
					#       ^ 5'UTR
					# - ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# - -------CDS----CDSCDS----CDS-------
					#                              ^ 5'UTR
					$csto=$sto;
					$csta=$csto-($C[4]-$C[3]);
					push @out,join("\t",$tid,$C[1],$C[6] eq "+" ? "5primeUTR" : "3primeUTR",$sta,$csta-1,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
				} else {
					# + ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# + -------CDS----CDSCDS----CDS-------
					#                              ^ 3'UTR
					# - ----E1E1E1----E2E2E2----E3E3E3---- (numbers in read order)
					# - -------CDS----CDSCDS----CDS-------
					#       ^ 3'UTR
					$csta=$sta;
					$csto=$csta+($C[4]-$C[3]);
					push @out,join("\t",$tid,$C[1],$C[6] eq "+" ? "3primeUTR" : "5primeUTR",$csto+1,$sto,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
				}
				push @out,join("\t",$tid,$C[1],$C[2],$csta,$csto,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
			}
		}
		@F=split/\t/,$gtf{$tid}->{transcript}->[0];
		$F[-1].=" genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";";
		$F[0]=$tid;
		$F[3]=1;
		$F[4]=$sto;
		unshift @out,join("\t",@F);
		$F[2]="gene";
		unshift @out,join("\t",@F);
		say TGGTF $_ for sort {my @a=split/\t/,$a; my @b=split/\t/,$b; $a[3] <=> $b[3] || $a[4] <=> $b[4]} @out;
		if($fasta){
			say TGFA ">$tid $t2g{$tid}";
			say TGFA $_ for unpack "(A100)*",$tseq;
		}

		####################################################
		# transcriptchromosome
		$sta=$gsta+1;
		$sto=0;
		@out=();
		for my $i (0..$gtf{$tid}->{Gexon}->$#*){
			@F=split/\t/,$gtf{$tid}->{Gexon}->[$i];
			$gseq.=$seq->subseq($F[3],$F[4]) if $fasta;
			$sto=$sta+($F[4]-$F[3]);

			push @out,join("\t","chrT",$F[1],$F[2],$sta,$sto,$F[5],$F[6],$F[7],$F[8]." genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";");

			# say STDERR $tid;
			if (exists $gtf{$tid}->{GCDS} && $gtf{$tid}->{GCDS}->[$i]){
				my $csta=0;
				my $csto=0;
				my @C=split/\t/,$gtf{$tid}->{GCDS}->[$i];
				if ($C[3] == $F[3] && $C[4] == $F[4]){
					$csta=$sta;
					$csto=$sto;
				}elsif($C[4] == $F[4]){
					$csto=$sto;
					$csta=$csto-($C[4]-$C[3]);
					push @out,join("\t","chrT",$C[1],$C[6] eq "+" ? "5primeUTR" : "3primeUTR",$sta,$csta-1,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
				} else {
					$csta=$sta;
					$csto=$csta+($C[4]-$C[3]);
					push @out,join("\t","chrT",$C[1],$C[6] eq "+" ? "3primeUTR" : "5primeUTR",$csto+1,$sto,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
				}
				push @out,join("\t","chrT",$C[1],$C[2],$csta,$csto,$C[5],$C[6],$C[7],$C[8]." genomic_position \"$C[0]:$C[3]-$C[4]\"; genomic_strand \"$C[6]\";");
			}
		}

		@F=split/\t/,$gtf{$tid}->{transcript}->[0];
		$F[-1].=" genomic_position \"$F[0]:$F[3]-$F[4]\"; genomic_strand \"$F[6]\";";
		$F[0]="chrT";
		$F[3]=$gsta+1;
		$F[4]=$sto;
		$F[6]="+";
		unshift @out,join("\t",@F);
		$F[2]="gene";
		unshift @out,join("\t",@F);
		say TCGTF $_ for sort {my @a=split/\t/,$a; my @b=split/\t/,$b; $a[3] <=> $b[3] || $a[4] <=> $b[4]} @out;
		$gseq.="N" x 1000;
		$gsta=$sto+1000;
	}
	if($fasta){
		say TCFA ">chrT";
		say TCFA $_ for unpack "(A100)*",$gseq;
	}

	close TFA;
	close TGTF;
	close TGFA;
	close TGGTF;
	close TCFA;
	close TCGTF;
} catch {
	&mydie($_);
};

say STDERR "success" if $verbose;
exit 0;
