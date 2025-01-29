#! /usr/bin/env perl
# (c) Konstantin Riege

# http://www.ensembl.org/Help/Glossary?id=346
# The canonical transcript is used in the gene tree analysis in Ensembl
# and does not necessarily reflect the most biologically relevant transcript
# of a gene. For human, the canonical transcript for a gene is set according
# to the following hierarchy: 1. Longest CCDS translation with no stop
# codons. 2. If no (1), choose the longest Ensembl/Havana merged translation
# with no stop codons. 3. If no (2), choose the longest translation with no
# stop codons. 4. If no translation, choose the longest non-protein-coding
# transcript.

# 1. grep for 'tag "CCDS"' -> multi transcripts? grep for $2==ensembl_havana -> multi transcripts? get longest
# 2. grep for $2==ensembl_havana -> multi transcripts? get longest
# 3. multi transcripts? get longest

use strict;
use warnings;
use feature ":5.10";
use List::Util qw(min max);

if ($#ARGV < 3){
	say "usage: tpm.pl <gtf> <feature-level> <feature-tag> <htseq.counts>";
	say 'extracts canonical transcripts searching for flags CCDS/*canonical, fallbacks to the longest meta-feature composed of sub-transcript features given the associated <feature-level> e.g. "exon"';
	say 'for each canonical transcript, the parental feature of interest is looked-up by its <feature-tag> e.g. "gene_id" or "gene_name" within the first column of supplied <htseq.counts> file';
	say 'the count value stored in the second column gets normalized by total count (library size) and canonical or longest transcript length';
	exit 1;
}

my $lvl=$ARGV[1];
my $ft=$ARGV[2];
my %m;
open F, "<".$ARGV[0] or die $!;
while(<F>){
	chomp;
	next if /^(\s*$|#)/;
	my @l = split /\t+/;
	next unless $l[2] eq $lvl;
	$l[-1]=~/$ft\s+(\S+)/;
	my $g = $1 ? $1 : $l[-1];
	$g=~s/("|;)//g;
	$l[-1]=~/transcript_id\s+(\S+)/;
	my $t = $1 ? $1=~s/("|;)//gr : 'transcript';
	if ($l[-1]=~/tag\s+\"(CCDS|[^"]*canonical)\"/) {
		$m{$g}{1}{$t} += $l[4]-$l[3]+1;
	} elsif ($l[2] eq 'ensembl_havana' || $l[2] eq 'HAVANA') {
		$m{$g}{2}{$t} += $l[4]-$l[3]+1;
	} else{
		$m{$g}{3}{$t} += $l[4]-$l[3]+1;
	}
}
close F;
my %id2len;
for my $g (keys %m){
	my $k = min(keys %{$m{$g}});
	my @ls;
	while (my ($t, $l) = each %{$m{$g}{$k}}) {
		push @ls,$l;
	}
	$id2len{$g} = max(@ls);
}

my @ids;
my %id2count;
my $f;
open F, "<".$ARGV[3] or die $!;
while(<F>){
	chomp;
	my @l = split /\s+/;
	unless(exists $id2len{$l[0]}){
		warn "$ARGV[0] does not contain $l[0]";
		next;
	}
	push @ids,$l[0];
	my $c = $l[1]/($id2len{$l[0]}/1000);
	$id2count{$l[0]} = $c;
	$f += $c;
}
close F;
$f /= 1000000;

say $_."\t".($id2count{$_}/$f) for @ids;

exit 0;
