#! /usr/bin/env perl
use strict;
use warnings;
use feature ":5.10";
use List::MoreUtils qw(indexes uniq);
use List::Util qw(min max sum);
use Getopt::Long;

sub usage {
    say "usage";
    say "-g|--germline (optional. use if multi sample vcf from germline calling)";
    say "-n|--normals <number> (optional. use if multi sample vcf from somatic calling) [1]";
    say "-i|--in <vcf> (required. needs to be last option!)";
    exit 1;
}
&usage if $#ARGV == -1;


# bcftools norm adapts info or format number=G|R i.e. e.g. TLOD info, GT, AD, PL, GL format during split of multiallelic sites but no other number=\d like DP4
# -> depends on present G and A keys in vcf header section

# Possible Types for INFO fields are: Integer, Float, Flag, Character, and String
# If the field has <X> values value should be <X>.
# If the field has one value per alternate allele then this value should be ‘A’.
# If the field has one value for each possible allele (including the reference), then this value should be ‘R’.
# If the field has one value for each possible genotype then this value should be ‘G’.
# If the number of possible values varies, is unknown, or is unbounded, then this value should be ‘.’.

my %header = (
    PL => '##FORMAT=<ID=PL,Number=G,Type=Float,Description="List of Phred-scaled genotype likelihoods">',
    GL => '##FORMAT=<ID=GL,Number=G,Type=Float,Description="List of genotype likelihoods">',
    NV => '##FORMAT=<ID=NV,Number=A,Type=Integer,Description="Number of reads containing variant in this sample">',
    NR => '##FORMAT=<ID=NR,Number=A,Type=Integer,Description="Number of reads covering variant location in this sample">',

    GQ => '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Phred-scaled genotype quality">',
    MAF => '##FORMAT=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequence">',
    COV => '##FORMAT=<ID=COV,Number=1,Type=Integer,Description="Position read coverage">',
    AD => '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
    DP => '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered read depth at the locus">',
    DP4 => '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Strand specific ref and alt read counts: ref-fwd, ref-rev, alt-fwd, alt-rev">',
    ASF => '##FORMAT=<ID=ASF,Number=1,Type=Float,Description="Alt strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
    RSF => '##FORMAT=<ID=RSF,Number=1,Type=Float,Description="Ref strand specific reads fraction: min(fwd,rev)/max(fwd/rev)">',
    vcfix_germlinerisk => '##FILTER=<ID=vcfix_germlinerisk,Description="if paired vcf, not marked as vcfix_somatic and vcfsamplediff reported germline">',
    vcfix_strandbias => '##FILTER=<ID=vcfix_strandbias,Description="ASF < 0.3">',
    vcfix_lowcov => '##FILTER=<ID=vcfix_lowcov,Description="COV < 10">',
    vcfix_lowqual => '##FILTER=<ID=vcfix_lowqual,Description="QUAL < 5">',
    vcfix_somatic => '##FILTER=<ID=vcfix_somatic,Description="if paired vcf, site is somatic according to the caller or vcfix QUAL and MAF fraction filter">',
    vcfix_dp4bias => '##FILTER=<ID=vcfix_dp4bias,Description="if source is platypus and normal samples significantly support variant, tumor alt dp4 values and thus ASF is biased">',
    vcfix_ok => '##FILTER=<ID=vcfix_ok,Description="site passed all naive vcfix filters">',
);

my $normals=1;
my $germline;

(Getopt::Long::Parser->new)->getoptions(
    'n|normals=i' => \$normals,
    'g|germline' => \$germline,
    'i|in' => sub {
        $normals+=8; # last normal index

        my $warnedgt;
        my $warneddp;
        my $warneddp4;
        my $warnedgq;
        while(<>){
            my $l = $_;
            chomp $l;
            my @l = split /\s+/,$l;
            if ($l[0]=~/^#/){
                if ($l[0]=~/^##FILTER=<ID=([^\s,]+)/) {
                    say exists $header{$1} ? $header{$1} : $l;
                    delete $header{$1};
                } elsif ($l[0]=~/^##FORMAT=<ID=([^\s,]+)/) {
                    say exists $header{$1} ? $header{$1} : $l;
                    delete $header{$1};
                } elsif ($l[0]=~/^#CHROM/) {
                    delete $header{PL};
                    delete $header{GL};
                    delete $header{NV};
                    delete $header{NR};
                    say $_ for values %header;
                    say $l;
                } else {
                    say $l;
                }
                next;
            }

            next if $l[5] ne "." && $l[5] < 1; # simple quality filter

            my @smaf;
            my @sgq;
            my %filter;
            my $somatic;
            my $vsupn = 0;
            my @t;
            $germline=1 if $#l==9;
            for my $s (9..$#l){
                @t = split /:/,$l[8]; # redo due to splice
                my @v = split /:/,$l[$s];
                my $i;


                # simplify gt
                ($i) = indexes { /^GT$/ } @t;
                if ($v[$i] =~ /\./){ # freebayes undefined genotype - sometimes happens if normal is poorly covered
                    warn "warning: unable to infer GT. output will be truncated." unless $warnedgt;
                    $warnedgt=1;
                    $filter{"vcfix_formaterror"}=1;
                    last;
                }
                my @gt = split/[\/\|]/,$v[$i];
                @gt = uniq(sort {$a <=> $b} @gt); # bcftools fix - simple sites 0/2 splitted into 1/0 and 0/1, multiallelic sites 0/1/2 splittet to 0/1/0 and 0/0/1
                push @gt, $gt[0] if $#gt == 0; # take care of 0/0, 1/1, ...
                $v[$i] = join '/', @gt;


                # insert or recalculate dp4
                my ($dp4) = indexes { /^DP4$/ } @t;
                my @dp4;
                if(($i) = indexes { /^SB$/ } @t){ #gatk4 haplotypecaller fix (requieres -A StrandBiasBySample)
                    @dp4 = split/,/,$v[$i];
                } elsif(($i) = indexes { /^RDF$/ } @t){ #varscan (germline) fix - ok, no multiallelic sites reported. attention: bcftools also reports ADF and ADR but as comma-sep. ref,alt
                    push @dp4 , $v[$i];
                    ($i) = indexes { /^RDR$/ } @t;
                    push @dp4 , $v[$i];
                    ($i) = indexes { /^ADF$/ } @t;
                    push @dp4 , $v[$i];
                    ($i) = indexes { /^ADR$/ } @t;
                    push @dp4 , $v[$i];
                } elsif(($i) = indexes { /^ADF$/ } @t){ #bcftools fix (below varscan fix)- do not use $l[7]=~/[\s;]DP4=([^\s;]+)/ which is the sum of all input samples
                    my @adf = split/,/,$v[$i];
                    ($i) = indexes { /^ADR$/ } @t;
                    my @adr = split/,/,$v[$i];
                    $dp4[0]=shift @adf;
                    $dp4[1]=shift @adr;
                    $dp4[2]=sum(@adf);
                    $dp4[3]=sum(@adr);
                } elsif($l[7]=~/[\s;]SRF=([^\s;]+)/){ #freebayes fix (info properly defined for bcftools norm). unlike platypus, reads are always variant exclusive (no read sharing within multiallelic cases)
                    if (! $germline && $s <= $normals){ # normal
                        @dp4=(-1,-1,-1,-1);
                    } else {
                        @dp4 = ($1);
                        push @dp4 , $1;
                        $l[7]=~/[\s;]SAF=([^\s;]+)/;
                        push @dp4 , sum(split/,/,$1);
                        $l[7]=~/[\s;]SAR=([^\s;]+)/;
                        push @dp4 , sum(split/,/,$1);
                    }
                } elsif($l[7]=~/[\s;]TCF=([^\s;]+)/){ #platypus fix
                    # if reads are shared among rare multiallelic cases, sum(NF) and/or sum(NR) can be larger than TCF and or TCR
                    # CG -> C,TG -----T----| |----  read (internally left-aligned)
                    #            NNNNNCGGGGGGGNNNN  ref
                    # note also that TC* and N* are total values across all input samples and thus is just works germline or single tumor normal pair, if normal does not support any variant
                    if (! $germline && $s <= $normals){ # normal
                        @dp4=(-1,-1,-1,-1);
                    } else {
                        my $allf = $1;
                        $l[7]=~/[\s;]TCR=([^\s;]+)/;
                        my $allr=$1;
                        $l[7]=~/[\s;]NF=([^\s;]+)/;
                        my $altf=sum(split/,/,$1);
                        $l[7]=~/[\s;]NR=([^\s;]+)/;
                        my $altr=sum(split/,/,$1);
                        if($germline) {
                            @dp4 = (max(-1,$allf-$altf),max(-1,$allr-$altr),$altf,$altr);
                        } else {
                            @dp4 = (-1,-1,$altf,$altr); # ref-fw ref-rv cannot be calculated, since unclear how to split and assign TCF/TCR to normal and tumor ref
                            $filter{"vcfix_dp4bias"}=1 if $s > $normals && $vsupn/($altf+$altr+1) >= 0.05;
                        }
                    }
                } elsif(($i) = indexes { /^ALD$/ } @t){ #vardict fix - ok, no multiallelic sites called
                    @dp4 = ($v[$i]);
                    ($i) = indexes { /^RD$/ } @t;
                    push @dp4 , $v[$i];
                } elsif (defined $dp4) { # varscan (somatic) fix and fallback
                     @dp4 = split/,/,$v[$dp4];
                } else {
                    warn "warning: unable to infer DP4. output will be truncated." unless $warneddp4;
                    $warneddp4 = 1;
                    $filter{"vcfix_formaterror"}=1;
                    last;
                }
                if(defined $dp4){ # update values, e.g. after splitting multiallelic hits
                    $v[$dp4] = join(",",@dp4);
                } else {
                    splice @t, 1, 0, 'DP4';
                    splice @v, 1, 0, join(",",@dp4);
                }


                # insert or correct ad
                ($i) = indexes { /^AD$/ } @t;
                if( ! defined $i){ # platypus fix
                    if (($i) = indexes { /^NR$/ } @t) {
                        my @nr = split /,/,$v[$i];
                        ($i) = indexes { /^NV$/ } @t;
                        my @ad = split /,/,$v[$i];

                        $vsupn += max(@ad) if $s <= $normals; # need for warning

                        $nr[$_]-=$ad[$_] for 0..$#nr;
                        unshift @ad, max(@nr);
                        splice @t, 1, 0, 'AD';
                        splice @v, 1, 0, join(",",@ad);
                    }
                } elsif (defined $i && $v[$i]!~/,/){ # varscanfix missing alt support (only correct for 1 or 2 alleles)
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


                # insert or recalculate dp
                my ($dpi) = indexes { /^DP$/ } @t;
                my $dp;
                if (($i) = indexes { /^AD$/ } @t) {
                    $dp = sum split /,/,$v[$i];
                } elsif (defined $dpi) { #fallback
                    $dp = $v[$dpi];
                } else {
                    warn "warning: unable to infer DP. output will be truncated." unless $warneddp;
                    $warneddp = 1;
                    $filter{"vcfix_formaterror"}=1;
                    last;
                }
                if(defined $dpi){ # update values, e.g. after splitting multiallelic hits
                    $v[$dpi] = $dp;
                } else {
                    splice @t, 1, 0, 'DP';
                    splice @v, 1, 0, $dp;
                }



                # add COV, MAV and strand usage fraction info
                ($i) = indexes { /^DP4$/ } @t;
                @dp4 = split/,/,$v[$i];
                ($i) = indexes { /^AD$/ } @t;
                my @ad = split /,/,$v[$i];
                ($i) = indexes { /^DP$/ } @t;

                # get cov
                my $cov = max(sum(@ad), $v[$i]); #try to find real COV, do not use dp4 since not present in freebayes, platypus tumor normal paired vcf
                $v[$i] = min(sum(@ad), $v[$i]); #varscan, vardict and freebayes fix DP >= sum(AD) i.e. DP ~ total COV, not filtered depth, whereas AD values are wrongly based on filtered reads

                # get maf
                shift @ad;
                my @maf;
                push @maf, $cov == 0 ? 0 : sprintf("%.4f",$_/$cov) for @ad;
                $_=~s/(\.[^0]*)0+$/$1/ for @maf;
                $_=~s/\.$// for @maf;

                # get asf and rsf
                my $rsf = max($dp4[0],$dp4[1]) == 0 ? 0 : max($dp4[0],$dp4[1]) == -1 ? -1 : sprintf("%.4f",min($dp4[0],$dp4[1])/max($dp4[0],$dp4[1]));
                $rsf=~s/(\.[^0]*)0+$/$1/;
                $rsf=~s/\.$//g;
                my $asf = max($dp4[2],$dp4[3]) == 0 ? 0 : max($dp4[2],$dp4[3]) == -1 ? -1 : sprintf("%.4f",min($dp4[2],$dp4[3])/max($dp4[2],$dp4[3]));
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


                # add or recalculate GQ (pred scaled)-genotype likelihoods

                # haplotypecaller info
                # https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
                # Genotype    A/A A/T T/T
                # Raw PL  -10 * log(0.000001) = 60  |  -10 * log(0.000100) = 40  |  -10 * log(0.010000) = 20
                # Normalized PL   60 - 20 = 40  |  40 - 20 = 20  |  20 - 20 = 0    (subtract 20 to get a 0)
                # genotype GQ is simply the difference between the second lowest PL and the lowest PL
                # -> gatk and freebayes cap at 99, bcftools at 127
                # Raw GL  0.000001  |  0.0001  |  0.01
                # log10-GL  -6  |  -4  |  -2
                # => GQ = sort(PL)[1 minus 0] or better min(PL) - PL[0]

                # handle germline call similar to somatic tumor call
                # PL: min(alts) - ref * -1 # var with highest gt prob (==1 aka 10^-0) minus prob of ref means to be 10^(diff/10) times more sure to be actually alt than ref
                # GL: max(alt) - ref * 10
                # handle normal of somatic call that way
                # PL: ref - min(alt) * -1 # von der erwarteten referenz den naechsten genotyp mit der howchsten wsk abziehen
                # GL: ref - max(alt) * 10

                # example
                # A->T  0/1  60,0,20
                # A->T  1/1  60,20,0
                # PL=> 0 - 60 * -1
                # GL=> 0 - -6 * 10
                # A->T  0/0  0,30,10  0/1  60,0,20
                # => 0 - 10 * -1 = NLOD
                # => 0 - 60 * -1 = TLOD
                # A->T,G  0/0 NLOD=10 0/1/2 TLOD=50,40 <- keep max(60,10,40 ; 60,20,40)

                # above is identical to tlod calc @ https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py
                # TumorLOD = max( map {$_ - $gl[0] } @gl[1..$#gl]) # in tumor, best likelihood not to be reference
                # -> max(@gl[1..$#gl]) - $gl[0]
                # NormalLOD = min( map {$gl[0] - $_ } @gl[1..$#gl]) # in normal, best likelihood to be reference
                # -> gl[0] - max(@gl[1..$#gl])

                # above also in agreement with
                # github.com/broadgsa/gatk-protected/blob/master/protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2/MuTect2.java
                # normalLod = normalGLs[0] - normalGLs[1]; tumorLod = tumorGLs[1] - tumorGLs[0];

                my $gq=999999; # vardict/varscan fix
                my @pl;
                if(($i) = indexes { $_ =~ /^PL$/ } @t){ # bcftools, haplotypecaller
                    @pl = split /,/,$v[$i];
                    if($germline){
                        $gq =  min(@pl[1..$#pl]) - $pl[0] * -1;
                    } elsif($s <= $normals){ # normal
                        $gq = $pl[0] - min(@pl[1..$#pl]) * -1;
                    } else { # tumor
                        $gq =  min(@pl[1..$#pl]) - $pl[0] * -1;
                    }
                } elsif(($i) = indexes { $_ =~ /^GL$/ } @t){ # freebayes, platypus
                    if ($v[$i] eq "-1,-1,-1"){ # platypus at multiallelic sites, should be rr,ra,rb,aa,ab,bb i.e. (n*[n-1])/2 , n=ref+#alt
                        my @alleles = split /,/,$l[4];
                        $v[$i] = join",",((-1)x( (($#alleles+3)*($#alleles+2))/2 ));
                        # one may calculate them like mutect does @ calcGenotypeLikelihoodsOfRefVsAny @ https://github.com/broadgsa/gatk-protected/blob/master/protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2/MuTect2.java
                        # for p in pileup: next unless p.isDeletion() || p.getQual > minBaseQual; x=1-10^-3; if isNonRef(refBase, p) : gl[ab] += log10( MAF*x + (1-MAF)*x/3 ); ql[aa] += log10( (1-x)/3 ) else gl[ab] = log10( MAF*(1-x)/3 + (1-MAF)*x); gl[aa] += log10(x)
                        # isNonRef : return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip()
                        # -> note for normal gatk uses MAF=0.5
                    } else {
                        @pl = split /,/,$v[$i];
                        $_ = 10*-$_ for @pl; # despite being defined as -log10 likelihoods (10**-$_), GL are like PL but not multiplied by 10. ( -10 * log(10^-$_)/log(10) )
                        # kinda proof: cap is ~ -300 which is the maximum of educible digits by double precision data type (2^-([2048-2]/2) = 2^-1023 or 10^-308 )
                        if($germline){
                            $gq =  min(@pl[1..$#pl]) - $pl[0] * -1;
                        } elsif($s <= $normals){ # normal
                            $gq = $pl[0] - min(@pl[1..$#pl]) * -1;
                        } else { # tumor
                            $gq =  min(@pl[1..$#pl]) - $pl[0] * -1;
                        }
                    }
                } elsif($l[7]=~/[\s;]TLOD=([^\s;]+)/){ # mutect
                    $gq = max(split /,/,$1);
                    if($s <= $normals){ # normal - TODO check if mutect can be run on multiple bams and how NLOD | TLOD changes
                        $l[7]=~/[\s;]NLOD=([^\s;]+)/;
                        $gq = max(split /,/,$1);
                    }
                } elsif (($i) = indexes { $_ =~ /^GQ$/ } @t){ # fallback
                    $gq = $v[$i] unless $v[$i] eq ".";
                } else {
                    warn "warning: unable to infer GQ. define GQ=999999." unless $warnedgq;
                    $warnedgq=1;
                }
                ($i) = indexes { $_ =~ /^GQ$/ } @t;
                if (defined $i){
                    $v[$i] = $gq;
                } else {
                    splice @t, 1, 0, 'GQ';
                    splice @v, 1, 0, $gq;
                }


                # add kinda decision helper
                push @smaf,max(@maf);
                push @sgq,max($gq);
                if( ! $germline && $s > $normals ){
                    # pcbio suggest to reject freebayes somatic calls if TumorLOD < 3.5 || NormalLOD < 3.5 AND Normal-MAF <0.001 || Normal-MAF < Normal-MAF/2.7
                    # -> values optimized according to DREAM synthetic benchmark dataset
                    # https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py
                    # -> somatic check for callers with paired input which do not report a status (bcftools, platypus, freebayes)
                    #if ( $sgq[0] >= 3.5 && $sgq[-1] >= 3.5 && ($smaf[0] <= 0.001 || $smaf[0] <= $smaf[-1]/2.7)){
                    if ( min(@sgq[0..$normals-9]) >= 3.5 && $sgq[-1] >= 3.5 && ( max(@smaf[0..$normals-9]) <= 0.001 || max(@smaf[0..$normals-9]) <= $smaf[-1]/2.7)){
                        $somatic=1 if max(@gt)>0; # otherwise could be LOH
                        $filter{"vcfix_germlinerisk"}=1 if $l[7]=~/=[gG]ermline/; # vcfsamplediff, vardict
                    }
                }
                if($germline || $s > $normals){
                    $filter{"vcfix_strandbias"}=1 if $asf < 0.3;
                }
                $filter{"vcfix_lowcov"}=1 if $cov < 10;


                $l[$s] = join ':' , @v;
            }

            next if exists $filter{"vcfix_formaterror"};

            $l[8] = join ':' , @t;

            $filter{"vcfix_lowqual"}=1 if $l[5] ne "." && $l[5] < 5;

            unless($germline){
                if($l[7]=~/[\s;]SS=(\d)[\s;]/){ # varscan
                    if ($1 == 2){
                        $somatic=1;
                    } elsif ($1 == 1){
                        $filter{"vcfix_germlinerisk"}=1;
                    }
                } elsif ($l[7]=~/[\s;]STATUS=[^;]*[sS]omatic/ && $l[6]=~/PASS/){ # vardict
                     $somatic=1;
                } elsif ($l[7]=~/[\s;]TLOD=/ && $l[6]=~/PASS/){ # mutect
                    $somatic=1;
                } elsif ($l[6]=~/germline/){ # mutect
                     $filter{"vcfix_germlinerisk"}=1;
                }
            }

            $filter{"vcfix_ok"}=1 if scalar(keys %filter) == 0 || (scalar(keys %filter) == 1 && exists $filter{"vcfix_dp4bias"});
            $filter{"vcfix_somatic"}=1 if $somatic;

            $filter{$_}=1 for split/;/,$l[6];
            delete $filter{"."}; # vcfix_dp4bias
            my @filter = sort {$a cmp $b} keys %filter;

            $l[6] = join";",@filter;
            $l = join "\t" , @l;
            say $l;
        }
    }
) or &usage;
