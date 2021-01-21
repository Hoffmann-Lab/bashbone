#! /usr/bin/env perl
# (c) Konstantin Riege
use feature ":5.10";

while(<>){
	next if /^\s*(#|$)/;
	@l=split/\t/;
	if ($l[2] eq "gene"){
		$l[-1]=~/gene_id\s+\"([^\"]+)/;
		@{$gpos{$1}}=($l[3],$l[4]);
	}
	next unless $l[2] eq "exon";
	$l[-1]=~/gene_id\s+\"([^\"]+)/;
	$g=$1;
	$l[-1]=~/transcript_id\s+\"([^\"]+)/;
	$t=$1;
	$l[-1]=~/exon_number\s+\"([^\"]+)/;
	$e=$1;
	$l[-1]="$t:$e";
	push @{$m{$g}}, join"\t",@l;
}

$all=scalar(keys %m);
for my $g (sort {$a cmp $b} keys %m){
	$c++;
	say STDERR "merged exons for $c of $all features" if $c%1000==0;
	$l = join"\n", sort {@a=split/\t/,$a; @b=split/\t/,$b; $a[3] <=> $b[3] || $b[4] <=> $b[4]} @{$m{$g}};
	$e=0;
	$i=0;
	$gsta=0;
	@o = split /\n/,`echo -e '$l' | bedtools merge -s -i - -c 6,7,8,9 -o distinct -delim '+'`;
	unless(@o){ #in case cmd becomes too long
		$f=$ENV{tmp}.'.exons.MERGERTMP';
		open O,">$f" or die $!;
		say O $l;
		close O;
		@o = split /\n/,`bedtools merge -s -i $f -c 6,7,8,9 -o distinct -delim '+'`;
	}
	for (@o){
		@l=split/\t/;
		$l[1]++; #bedtools outcome is 0-based
		# now handle sloped genes, exons respectively
		if ($e == 0){
			$gsta = $gpos{$g}[0];
			$l[1] = $gsta if $gsta < $l[1];
			$gsta = $l[1];
			$t = "first";
		}
		if ($e == $#o){
			$gsto = $gpos{$g}[1];
			$l[2] = $gsto if $gsto > $l[2];
			$gsto = $l[2];
			$t = "last";
			$t = "single" if $e == 0;
		}
		if ($e > 0 && $e < $#o){
			$t = "inner";
		}
		$e++;

		$l[-1]='feature_id "exon:'.$e.'@'.$g.'"; gene_id "'.$g.'"; full_id "'.$l[-1].'"; feature_tag "'.$t.'";';
		say join"\t",($l[0],"merger","exon",@l[1..$#l]);

		if ($i > 0){
			$t = "first"; #if $i == 1;
			$t = "inner" if $i > 1;
			$t = "last" if $i > 1 && $e-1 == $#o;
			$t = "single" if $i==1 && $e-1 == $#o;
			# from previous esto to current sta
			say join"\t",($l[0],"merger","intron",++$esto,--$l[1],@l[3..5],'feature_id "intron:'.$i.'@'.$g.'"; gene_id "'.$g.'"; feature_tag "'.$t.'";');
		}
		$i++;
		$esto=$l[2];
	}
	if($gsta){
		say join"\t",($l[0],"merger","gene",$gsta,$gsto,@l[3..5],'feature_id "gene:'.$g.'"; gene_id "'.$g.'";');
	} else {
		say STDERR "ERROR shit happened for feature $g";
	}
}
say STDERR "merged exons for $c of $all features" if $c%1000==0;
