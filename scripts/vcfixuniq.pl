#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use List::MoreUtils qw(indexes);
use List::Util qw(sum);

my @mem;
while(<>){
    my $l = $_;
    chomp $l;
    my @l = split /\s+/,$l;
    if ($l[0]=~/^#/){
        say $l;
        next;
    }
    if ($#mem == -1){
        push @mem,$l;
    } else {
        my @m = split /\s+/,$mem[-1];
        unless ($l[0] eq $m[0] && $l[1] == $m[1]){
            &select;
            @mem=();
        }
        push @mem,$l;
    }
}
&select;

sub select (){
    my %gq;
    for my $m (@mem){
        my @m = split /\s+/ , $m;
        my ($i) = indexes { $_ =~ /^GQ$/i } split /:/,$m[8];
        my ($j) = indexes { $_ =~ /^MAF$/i } split /:/,$m[8];
        my @v = split /:/,$m[-1];
        $m[-1] = join ":",@v;
        $m = join "\t",@m;
        $gq{$v[$i]}{$v[$j]} = $m;
    }

    for my $i (reverse sort {$a <=> $b} keys %gq){
        for my $j (reverse sort {$a <=> $b} keys %{$gq{$i}}){
            say $gq{$i}{$j};
            last;
        }
        last;
    }
}
