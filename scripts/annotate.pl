#! /usr/bin/env perl
# (c) Konstantin Riege
use strict;
use warnings;
use feature ":5.10";

if ($#ARGV < 3){
	say "usage: annotate.pl <gtf.descr|0> <gtf|0> <feature> <bed|deseq.tsv|heatmap.ps> [<bed|deseq.tsv|heatmap.ps> ..] [tpm]";
	say "0 or gtf.descr => 4 tab seperated columns: geneID geneName biotype description";
	say "0 or gtf => parsed for gene_id gene_name gene_biotype";
	say "feature => e.g. exon -> gtf parsed for exon_id exon_name exon_biotype";
	exit 1;
}
my $ft=$ARGV[2];
my %m;
my %mps;
if ($ARGV[0]){
	open I,"<$ARGV[0]" or die $!;
	while(<I>){
		chomp;
		my @l=split/\t/;
		push @l,"" while $#l<3;
		$l[3]=~s/\s+\[Source.+//;
		$l[3]="NA" if $l[3]=~/^\s*$/;
		$m{$l[0]}=join"\t",@l[1..3]; # name, biotype, description
		$mps{$l[0]}=$l[1];
	}
	close I;
}

if ($ARGV[1]){
	open I,"<$ARGV[1]" or die $!;
		while(<I>){
			chomp;
			next if /^(\s*$|#)/;
			my @l=split/\t/;
			next unless $l[2] eq $ft;
			my ($g,$n,$b) = ('','','');
			if ($l[-1]=~/${ft}_id\s+\"([^\"]+)/){
				$g = $1;
				next if exists $mps{$g};
				if ($l[-1]=~/(${ft}_n|N)ame\s+\"([^\"]+)/){
					$n = $+;
					$mps{$g} = $n;
				}
				if ($l[-1]=~/${ft}_(bio)?type\s+\"([^\"]+)/){
					$b = $+;
				}
				$m{$g} = join"\t",($n,$b,'NA');
			}
		}
	close I;
}

my $files = $#ARGV == 3 ? 3 : $#ARGV-1;
my %tpm;
if ($files < $#ARGV){
	open I,"<$ARGV[-1]" or die $!;
	while(<I>){
		chomp;
		my @l=split/\t/;
		$tpm{$l[0]} = join"\t",@l[1..$#l];
	}
	close I;
	$tpm{NA}=$tpm{id}=~s/\S+/NA/gr;
}

for (3..$files){
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
				if ($files < $#ARGV){
					say O join"\t",(@l,"name","biotype","description",$tpm{id});
				} else {
					say O join"\t",(@l,"name","biotype","description");
				}
				next;
			} # else: not next!
		}
		# do splits by '@' in case of merged gtf with ids of type feature@subfeature@subfeature...
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
			$n = join"\t",('NA','NA','NA') unless $n;
			if ($files < $#ARGV){
				my $t = $tpm{$l[0]};
				$t = $tpm{NA} unless $t;
				$n .= "\t$t";
			}
			say O join"\t",(@l,$n);
		} else {
			my $n = $m{$l[3]};
			$n = $m{(split/\@/,$l[3])[-1]} unless $n;
			$n = join"\t",('NA','NA','NA') unless $n;
			say O join"\t",(@l,$n);
		}
	}
	close I;
	close O;
}
