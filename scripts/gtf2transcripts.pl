#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use Getopt::Long;

sub help {
	print <<EOF;
Synopsis: 
gtf4dexseq.pl [options] annotation.gtf

Description:
Converts a GTF file into a transcript based GTF and/or DEXSeq readable annotation format, both suitable for featureCounts based quantification.
I.e. each transcript is treated as a gene with its related exons. Please note, no CDS or other features will be printed.

Requirements:
Input GTF needs to besorted by -k1,1 -k4,4n -k5,5n
GTF feature exon with info fields: transcript_id, gene_id, exon_number, (optional: gene_biotype)

Optional functions:
1) experimental - TODO: Better truncate exons according to CDS
UTR overlapping exons can be excluded from the output to reduce the risk of finding some of them as differntially expressed due to alternative start/stop sites.
Required feature type at column 3: a string containing 'UTR' or 'utr' e.g. five_prime_utr or 3primeUTR

2) filter GTF according to a given gene_biotype
E.g. filter for commonly poly-A transcriptsby a regular expression: (protein_coding|TEC|antisense|lincRNA|sense_intronic|sense_overlapping)
The GTF info fields will be searched for /gene_biotype "<regex>"/

Options:
-o|--outgtf <file>    (required unless -d)
-d|--outdexseq <file> (required unless -o)
-b|--biotype <regex>  (optional)
-x|--excludeutr       (optional)
-h|--help"
EOF
	exit 0;
}

my (%utr, $excludeutr, %exonnumber, %fakeid, %genesgtf, %genesdex, %exonsdex, %exonsgtf, $filedex, $filegtf, @ids);
my $biotype='';
my $fakeidc=0;

&help if $#ARGV == -1;
(Getopt::Long::Parser->new)->getoptions(
	'o|outgtf:s' => \$filegtf,
	'd|outdexseq=s' => \$filedex,
	'x|excludeutr' => \$excludeutr,
	'b|biotype:s' => sub {
		$biotype = qr/gene_biotype\s+\"$_[1]\"/;
	},
	'h|help' => \&help
) or &help;
&help unless $filegtf || $filedex;

my @in;
while(<>){
	chomp;
	my @l = split/\t/;
	if($excludeutr && $l[2]=~/utr/i){
		$utr{$l[6]}{$l[3]}{$l[4]}=1;
	} elsif ($l[2] eq "exon") {
		if($biotype){
			push @in,$_ if $l[-1]=~/gene_biotype/ && $l[-1]=~/$biotype/;
		} else {
			push @in,$_;
		}
	}
}

for (@in){
	my @l = split/\t/;
	next if exists $utr{$l[6]}{$l[3]}{$l[4]};

	$l[-1]=~/gene_id\s+\"([^\"]+)/ or die "gene_id missing - no valid GTF";
	my $gid=$1;
	my $tid;
	if ($l[-1]=~/transcript_id\s+\"([^\"]+)/){
		$tid=$1;
	} else {
		if (exists $fakeid{$gid}){
			$tid = $fakeid{$gid};
		} else {
			$fakeidc++;
			$tid = join("","T",(0)x(11-length($fakeidc)),$fakeidc);
			$fakeid{$gid}=$tid;
		}
	}

	unless (exists $exonnumber{"$tid\@$gid"}){
		$genesdex{"$tid\@$gid"} = [$l[0],"gtf4dexseq.pl","aggregate_gene",$l[3],0,".",$l[6],".","gene_id \"$tid\@$gid\""];
		$genesgtf{"$tid\@$gid"} = [$l[0],"gtf4dexseq.pl","gene",$l[3],0,".",$l[6],".","gene_id \"$tid\@$gid\""];
		push @ids, "$tid\@$gid";
	}
	$genesdex{"$tid\@$gid"}[4] = $l[4];
	$genesgtf{"$tid\@$gid"}[4] = $l[4];
	$exonnumber{"$tid\@$gid"}++;

	my $n = join("",(0)x(3-length($exonnumber{"$tid\@$gid"})),$exonnumber{"$tid\@$gid"});
	push @{$exonsdex{"$tid\@$gid"}}, [$l[0],"gtf4dexseq.pl","exonic_part",$l[3],$l[4],".",$l[6],".","transcripts \"$tid\"; exonic_part_number \"$n\"; gene_id \"$tid\@$gid\""];
	push @{$exonsgtf{"$tid\@$gid"}}, [$l[0],"gtf4dexseq.pl","exon",$l[3],$l[4],".",$l[6],".","transcript_id \"$tid\"; exon_number \"$n\"; gene_id \"$tid\@$gid\"; exon_id \"$tid\@$gid:$n\""];
}

open GTF, ">$filegtf" or die $! if $filegtf;
open DEX, ">$filedex" or die $! if $filedex;

for my $id (@ids){
	if ($#{$exonsdex{$id}} > 1){
		if($filedex){
			say DEX join("\t",@{$genesdex{$id}});
			say DEX join("\t",@$_) for @{$exonsdex{$id}};
		}
		if ($filegtf){
			say GTF join("\t",@{$genesgtf{$id}});
            say GTF join("\t",@$_) for @{$exonsgtf{$id}};
		}
	}
}

close DEX if $filedex;
close GTF if $filegtf;
