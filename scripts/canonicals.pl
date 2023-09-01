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

if ($#ARGV < 2){
	say "usage: canonicals.pl file.gtf feature-level feature-tag";
	say "gtf needs to contain feature level e.g. exon";
	say "and feature id tag e.g. gene_id with transcript_id";
	exit 1;
}

my $lvl=$ARGV[1];
my $ft=$ARGV[2];
my %m;
my @ids;
open F, "<".$ARGV[0] or die $!;
while(<F>){
	chomp;
	next if /^(\s*$|#)/;
	my @l = split /\t+/;
	next unless $l[2] eq $lvl;
	$l[-1]=~/$ft\s+(\S+)/;
	my $g = $1 ? $1 : $l[-1];
	$g=~s/("|;)//g;
	push @ids,$g unless exists $m{$g};
	$l[-1]=~/transcript_id\s+(\S+)/;
	my $t = $1 ? $1 : 'transcript';
	$t=~s/("|;)//g;
	if ($l[-1]=~/tag\s+\"CCDS\"/) {
		$m{$g}{1}{$t} += $l[4]-$l[3]+1;
	} elsif ($l[2] eq 'ensembl_havana') {
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
	my @ts;
	while (my ($t, $l) = each %{$m{$g}{$k}}) {
		push @ls,$l;
		push @ts,$t;
	}
	my $max = max(@ls);
	say $g."\t".$ts[(grep {$ls[$_]==$max} 0..$#ls)[0]];
}

exit 0;
