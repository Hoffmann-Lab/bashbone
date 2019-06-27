#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use Getopt::Long;

my $args = $#ARGV;

my $biotype='';
my %utr;
my $excludeutr;
my %exonnumber;
my %fakeid;
my $fakeidc=0;
my %genesgtf;
my %genesdex;
my %exonsdex;
my %exonsgtf;
my $filedex;
my $filegtf;
my @ids;

sub help {
	print <<EOF;
Synopsis: 
gtf4dexseq.pl [options] annotation.gtf

Description:
Converts a GTF file into a DEXSeq readable annotation format.

Required GTF info fields (column 9): gene_id, exon_number 
Recommended, additional info fields: transcript_id and gene_biotype to filter features for

UTR overlapping exons can be excluded from the output to reduce the risk of finding some of them as differntially expressed just by alternative start/stop sites.
Required feature type (column 3): a substring matching UTR or utr e.g. five_prime_utr

Options:
-d|--outdexseq <file> (required)
-g|--outgtf <file>    (optional)
-x|--excludeutr       (optional)
-h|--help"
EOF
	exit 0;
}

(Getopt::Long::Parser->new)->getoptions(
	'g|outgtf:s' => \$filegtf,
	'd|outdexseq=s' => \$filedex,
	'x|excludeutr' => \$excludeutr,
	'b|biotype:s' => sub {
		$biotype = qr/gene_biotype\s+\"$_[1]\"/;
	},
	'h|help' => \&help
);
&help if $args < 2;

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
		$genesdex{"$tid\@$gid"} = [$l[0],"gtf4dexseq.pl","aggregate_gene",$l[3],$l[3],".",$l[6],".","gene_id \"$tid\@$gid\""];
		$genesgtf{"$tid\@$gid"} = [$l[0],"gtf4dexseq.pl","gene",$l[3],$l[3],".",$l[6],".","gene_id \"$tid\@$gid\""];
		push @ids, "$tid\@$gid";
	}
	$genesdex{"$tid\@$gid"}[4] = $l[4];
	$exonnumber{"$tid\@$gid"}++;

	my $n = join("",(0)x(3-length($exonnumber{"$tid\@$gid"})),$exonnumber{"$tid\@$gid"});
	push @{$exonsdex{"$tid\@$gid"}}, [$l[0],"gtf4dexseq.pl","exonic_part",$l[3],$l[4],".",$l[6],".","transcripts \"$tid\"; exonic_part_number \"$n\"; gene_id \"$tid\@$gid\""];
	push @{$exonsgtf{"$tid\@$gid"}}, [$l[0],"gtf4dexseq.pl","exon",$l[3],$l[4],".",$l[6],".","transcript_id \"$tid\"; exon_number \"$n\"; gene_id \"$tid\@$gid\"; exon_id \"$tid\@$gid:$n\""];
}

open GTF, ">$filegtf" or die $! if $filegtf;
open DEX, ">$filedex" or die $!;

for my $id (@ids){
	if ($#{$exonsdex{$id}} > 1){
		say DEX join("\t",@{$genesdex{$id}});
		say DEX join("\t",@$_) for @{$exonsdex{$id}};
		if ($filegtf){
			say GTF join("\t",@{$genesgtf{$id}});
            say GTF join("\t",@$_) for @{$exonsgtf{$id}};
		}
	}
}

close DEX;
close GTF if $filegtf;
