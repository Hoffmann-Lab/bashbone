#! /usr/bin/env perl
# (c) Konstantin Riege
use strict;
use warnings;
use feature ":5.10";

if ($#ARGV < 2){
	say "usage: annotate.pl [gtf.descr|0] [gtf|0] [bed|deseq.tsv|heatmap.ps| ..]";
	say "0 or gtf.descr => 4 tab seperated columns: geneID geneName biotype description";
	say "0 or gtf => parsed for gene_id gene_name gene_biotype";
	exit 1;	
}

my %m;
my %mps;
if ($ARGV[0]){
	open I,"<$ARGV[0]" or die $!;
	while(<I>){
		chomp;
		my @l=split/\t/;
		push @l,"" while $#l<3;
		$l[3]=~s/\s+\[Source.+//;
		$m{$l[0]}=join"\t",@l[1..3]; # name, biotype, description
		$mps{$l[0]}=$l[1];
	}
	close I;
}

if ($ARGV[1]){
	open I,"<$ARGV[1]" or die $!;
		while(<I>){
			chomp;
			my @l=split/\t/;
			next unless $l[2] eq 'gene';
			my ($g,$n,$b) = ('','','');
			if ($l[-1]=~/gene_id\s+\"([^\"]+)/){
				$g = $1;
				next if exists $mps{$g};
				if ($l[-1]=~/gene_name\s+\"([^\"]+)/){
					$n = $1;
					$mps{$g} = $n;
				}
				if ($l[-1]=~/gene_biotype\s+\"([^\"]+)/){
					$b = $1;
				}
				$m{$g} = join"\t",($n,$b,$n);
			}
		}
	close I;
}

for (2..$#ARGV){
	my $f=$ARGV[$_];
	my @o=split/\./,$f;
	my $e=$o[-1];
	$o[-1]="annotated";
	push @o,$e;
	my $o=join".",@o;
	open I,"<$f" or die $!;
	open O,">$o" or die $!;
	my $ps;
	my $deseq;
	while(<I>){
		chomp;
		my @l=split/\t/;
		if ($.==1){
			if (/PS-Adobe/){
				$ps=1;
				say O;
				next;
			} elsif(/log2FoldChange/){
				$deseq=1;
				say O join"\t",(@l,"geneName","biotype","description");
				next;
			} # else: not next!
		}
		if ($ps) {
			if (/\((\S+)(\s\+|\s-)*\)\s+(0|0\.25)\s+0\s+t$/) {
				my $n = $mps{$1};
				$n = $mps{(split/\@/,$1)[-1]} unless $n;
				s/$1/$n/ if $n;
			}
			say O;
		} elsif($deseq) {
			my $n = $m{$l[0]};
			$n = $m{(split/\@/,$l[0])[-1]} unless $n;
			$n = join"\t",('','','') unless $n;
			say O join"\t",(@l,$n);
		} else {
			my $n = $m{$l[3]};
			$n = join"\t",('','','') unless $n;
			say O join"\t",(@l,$n);
		}
	}
	close I;
	close O;
}
