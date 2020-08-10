#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use List::MoreUtils qw(indexes uniq);
use List::Util qw(min max sum);
use Getopt::Long;
die "option -i|--in is missing" if $#ARGV == -1;

# Possible Types for INFO fields are: Integer, Float, Flag, Character, and String
# If the field has <X> values value should be <X>.
# If the field has one value per alternate allele then this value should be ‘A’.
# If the field has one value for each possible allele (including the reference), then this value should be ‘R’.
# If the field has one value for each possible genotype then this value should be ‘G’.
# If the number of possible values varies, is unknown, or is unbounded, then this value should be ‘.’.
#!!! Number=G dont work for bcftools < v1.9 - better use . instead?
my %format = (
    GQ => '##FORMAT=<ID=GQ,Number=A,Type=Float,Description="Phred-scaled genotype quality">',
    PL => '##FORMAT=<ID=PL,Number=G,Type=Float,Description="List of Phred-scaled genotype likelihoods">', 
    GL => '##FORMAT=<ID=GL,Number=G,Type=Float,Description="List of genotype likelihoods">',
    MAF => '##FORMAT=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequence">',
    COV => '##FORMAT=<ID=COV,Number=1,Type=Integer,Description="Position read coverage">',
    AD => '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
    DP => '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth at the locus">',
    DP4 => '##FORMAT=<ID=DP4,Number=4,Type=Float,Description="Strand specific ref and alt read counts: ref-fwd, ref-rev, alt-fwd, alt-rev">',
    ASF => '##FORMAT=<ID=ASF,Number=1,Type=Float,Description="Alt strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
    RSF => '##FORMAT=<ID=RSF,Number=1,Type=Float,Description="Ref strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
);
my $caller='';

(Getopt::Long::Parser->new)->getoptions(
    'c|caller=s' => \$caller,
    'i|in' => sub {
        while(<>){
            my $l = $_;
            chomp $l;
            my @l = split /\s+/,$l;
            if ($l[0]=~/^#/){
                if ($l[0]=~/^##FORMAT=<ID=([^,]+)/) {
                    say exists $format{$1} ? $format{$1} : $l;
                    delete $format{$1};
                } elsif ($l[0]=~/^#CHROM/) {
                    delete $format{PL};
                    delete $format{GL};
                    say $_ for values %format;
                    say $l;
                } else {
                    say $l;
                }
                next;
            }

            # skip undefined GT
            next if $l[-1] =~ /^\./;

            my @t;
            for my $s (9..$#l){
                @t = split /:/,$l[8];
                my @v = split /:/,$l[$s];
                my $i;

                #bcftools norm adapts GT, DP, AD, PL, GL during slipt of multiallelic sites
                ($i) = indexes { /^GT$/ } @t;
                my @gt = split/[\/\|]/,$v[$i];
                @gt = uniq(sort {$a <=> $b} @gt); # bcftools fix - simple sites 0/2 splitted into 1/0 and 0/1, multiallelic sites 0/1/2 splitet to 0/1/0 and 0/0/1
                push @gt, $gt[0] if $#gt == 0; # take care of 0/0, 1/1, ...
                $v[$i] = join '/', @gt;

                my ($dp4) = indexes { /^DP4$/ } @t;
                if(! defined $dp4){
                    my @dp4;
                    if(($i) = indexes { /^SB$/ } @t){ #gatk4 haplotypecaller fix (requieres -A StrandBiasBySample)
                        @dp4 = split/,/,$v[$i];
                    } elsif($l[7]=~/[\s;]SRF=([^\s;]+)/){ #freebayes fix
                        @dp4 = ($1);
                        $l[7]=~/[\s;]SRR=([^\s;]+)/;
                        push @dp4 , $1;
                        $l[7]=~/[\s;]SAF=([^\s;]+)/;
                        push @dp4 , sum(split/,/,$1);
                        $l[7]=~/[\s;]SAR=([^\s;]+)/;
                        push @dp4 , sum(split/,/,$1);
                    } elsif($l[7]=~/[\s;]TCF=([^\s;]+)/){ #platypus fix
                        my $allf = $1;
                        $l[7]=~/[\s;]TCR=([^\s;]+)/;
                        my $allr=$1;
                        $l[7]=~/[\s;]NF=([^\s;]+)/;
                        my $altf=sum(split/,/,$1);
                        $l[7]=~/[\s;]NR=([^\s;]+)/;
                        my $altr=sum(split/,/,$1);
                        @dp4 = ($allf-$altf,$allr-$altr,$altf,$altr);
                    } elsif(($i) = indexes { /^ALD$/ } @t){ #vardict fix
                        @dp4 = ($v[$i]);
                        ($i) = indexes { /^RD$/ } @t;
                        push @dp4 , $v[$i];
                    } else {
                        next;
                    }

                    if(defined $dp4){
                        $v[$dp4] = join(",",@dp4);
                    } else {
                        splice @t, 1, 0, 'DP4';
                        splice @v, 1, 0, join(",",@dp4);
                    }
                }
                
                ($i) = indexes { /^AD$/ } @t;                
                if (defined $i && $v[$i]!~/,/){ #varscanfix (only correct for 1 or 2 alleles)
                    if (my ($j) = indexes { /^RD$/ } @t){
                        my @alleles = split /,/,$l[4];
                        pop @alleles;
                        my @ad;
                        push @ad, $v[$j];
                        push @ad, $v[$j]-$v[$i] for @alleles;
                        push @ad, $v[$i];
                        $v[$i] = join ",",@ad;
                    }
                }

                unless (grep { $_ =~ /^DP$/ } @t){
                    if (($i) = indexes { /^NR$/ } @t) { #platypus fix
                        my ($j) = indexes { /^NV$/ } @t;
                        my $dp = max(split /,/,$v[$i]);
                        my @ad = split /,/,$v[$j];
                        unshift @ad, $dp-sum(@ad);
                        splice @t, 1, 0, 'DP';
                        splice @t, 1, 0, 'AD';
                        splice @v, 1, 0, $dp;
                        splice @v, 1, 0, join(",",@ad);
                    } elsif (($i) = indexes { /^DP4$/ } @t) {
                        my $dp = sum split /,/,$v[$i];
                        splice @v, 1, 0, $dp;
                        splice @t, 1, 0, 'DP';
                    } elsif (($i) = indexes { /^AD$/ } @t) {
                        my $dp = sum split /,/,$v[$i];
                        splice @v, 1, 0, $dp;
                        splice @t, 1, 0, 'DP';
                    } else {
                        next;
                    }
                }

                # add or recalculate GQ as difference from best and second best observed (pred scaled)-genotype likelihood
                # NOT CORRECT FOR NORMAL SAMPLES IN CASE OF SOMATIC INPUT
                my $gq=0;
                if(($i) = indexes { $_ =~ /^PL$/ } @t){
                    my @q = sort {$a <=> $b} split /,/,$v[$i];
                    $gq = sprintf("%.4f",$q[1] - $q[0]); 
                } elsif(($i) = indexes { $_ =~ /^GL$/ } @t){
                    my @pl = split /,/,$v[$i];
                    $_ = ($_*-1)/10 for @pl;
                    my @q = sort {$a <=> $b} @pl;
                } elsif($l[7]=~/[\s;]TLOD=([^\s;]+)/){ #gatk mutect2 fix - INFO splitted by bcftools 
                    $gq = $1;
                } elsif (($i) = indexes { $_ =~ /^GQ$/ } @t){
                    $gq = $v[$i];
                }
                my @gq = split /,/,$gq; #in case of multiallelic sites
                for(@gq){
                    $_ = sprintf("%.4f",$_); #fix for vt normalize which cannot handle E-10 like representation
                    $_ =~ s/0+$//g; # trim trailing zeros
                    $_ =~ s/\.$//g;    
                }
                $gq = join ',' , @gq;
                ($i) = indexes { $_ =~ /^GQ$/ } @t; # replace after recalculated from bcftools splittet PL fields
                if (defined $i){
                    $v[$i] = $gq;
                } else {
                    splice @t, 1, 0, 'GQ';
                    splice @v, 1, 0, $gq;
                }

                # finally add COV, MAV and strand usage fraction info tags
                ($i) = indexes { /^DP$/ } @t;
                my ($j) = indexes { /^DP4$/ } @t;
                my @dp4 = split/,/,$v[$j];
                ($j) = indexes { /^AD$/ } @t;
                my @ad = split /,/,$v[$j]; 

                $v[$i] = min(sum(@ad), $v[$i]); #varscan, vardict and freebayes fix DP >= sum(AD|DP4) - (== total COV, not filtered depth), whereas AD values are based on filtered reads
                my $cov = max(sum(@ad), sum(@dp4), $v[$i]); #try to find real COV
                shift @ad;
                my @maf;
                push @maf, $cov == 0 ? 0 : sprintf("%.4f",$_/$cov) for @ad;
                $_=~s/(\.[^0]*)0+$/$1/ for @maf;
                $_=~s/\.$// for @maf;
                my $rsf = max($dp4[0],$dp4[1]) == 0 ? 0 : sprintf("%.4f",min($dp4[0],$dp4[1])/max($dp4[0],$dp4[1]));
                $rsf=~s/(\.[^0]*)0+$/$1/;
                $rsf=~s/\.$//g;
                my $asf = max($dp4[2],$dp4[3]) == 0 ? 0 : sprintf("%.4f",min($dp4[2],$dp4[3])/max($dp4[2],$dp4[3]));
                $asf=~s/(\.[^0]*)0+$/$1/;
                $asf=~s/\.$//;

                ($i) = indexes { /^COV$/ } @t;
                if(defined $i){
                    $v[$i]=$cov
                } else {
                    splice @t, 1, 0, 'COV';
                    splice @v, 1, 0, $cov;
                }
                ($i) = indexes { /^MAF$/ } @t;
                if(defined $i){
                    $v[$i]=join(",",@maf);
                } else {
                    splice @t, 1, 0, 'MAF';
                    splice @v, 1, 0, join(",",@maf);
                }
                ($i) = indexes { /^ASF$/ } @t;
                if(defined $i){
                    $v[$i]=$asf;
                } else {
                    splice @t, 1, 0, 'ASF';
                    splice @v, 1, 0, $asf;
                }
                ($i) = indexes { /^RSF$/ } @t;
                if(defined $i){
                    $v[$i]=$rsf;
                } else {
                    splice @t, 1, 0, 'RSF';
                    splice @v, 1, 0, $rsf;
                }

                # print
                $l[$s] = join ':' , @v;
            }
            $l[8] = join ':' , @t;
            $l = join "\t" , @l;
            say $l;
        }
    }
) or die "Unkown option";
