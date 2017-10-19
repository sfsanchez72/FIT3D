#!/usr/bin/perl
#
#
#
use PGPLOT;


use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;

# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fix one parameter to do
# a force fitting latter
if ($#ARGV<6) {
    print "USE: poly_force.pl PARAM MODEL POLY_ORDER N_MIN N_MAX INPUT_FILE OUT_FILE\n";        
    exit;
}

$n_param=$ARGV[0];
$n_model=$ARGV[1];
$npoly=$ARGV[2];
$n_min=$ARGV[3];
$n_max=$ARGV[4];
$input_file=$ARGV[5];
$out_file=$ARGV[6];
$dev='?';
if ($#ARGV==7) {
    $dev=$ARGV[7];
}

#
# Reading the input file
#
$y_min=1000000;
$y_max=-1000000;
$n=0;
open(INPUT,"<$input_file");
while($input=<INPUT>) {
    chop($input);
    $n_mod=$input;
    for ($i=0;$i<$n_mod;$i++) {
	$input=<INPUT>; chop($input);
	@data=split(" ",$input);
	$model[$i]=$data[0];
	for ($j=0;$j<9;$j++) {
	    $a[$i][$j][$n]=$data[1+2*$j];
	    $ia[$i][$j][$n]=$data[1+2*$j+1];
	    if (($j==$n_param)&&($i==$n_model)) {
		$a_fix[$n]=$a[$i][$j][$n];
		if ($y_min>$a_fix[$n]) {
		    $y_min=$a_fix[$n];
		}
		if ($y_max<$a_fix[$n]) {
		    $y_max=$a_fix[$n];
		}
	    }
	}
    }
    $x[$n]=$n;
    $n++;
};
close(INPUT);

#
# We fit to a polynomical function
# of the order npoly (using PDL);


$y_pdl = pdl(@a_fix);
($s_y,$coeff) = fitpoly1d $y_pdl,$npoly;
my @s_x_out;
for ($j=0;$j<$n;$j++) {
    $s_y_out[$j] = $s_y->slice($j)->sclr;
    $sigma_y[$j]=sqrt(($a_fix[$j]-$s_y_out[$j])**2);
}

$nn=0;

open(OUT,">junk.txt");
my @xx; my @yy; my @w;
for ($j=$n_min;$j<$n_max;$j++) {
    $xx[$nn] = $j;
    $yy[$nn] = $a_fix[$j];
    $w[$nn]=0.1;
    print OUT "$nn $xx[$nn] $yy[$nn]\n";
    $nn++;
}
close(OUT);

#print "NN=$nn\n";
$xx_pdl = pdl(@xx);
$yy_pdl = pdl(@yy);
#($s_y2,$coeff2) = fitpoly1d($xx_pdl,$yy_pdl,$npoly);#, {Weights => pdl(@w)});
#my $c;

($s_y2,$coeff2) = fitpoly1d(pdl(@xx),pdl(@yy),$npoly);#, {Weights => pdl(@w)});

for ($j=0;$j<$npoly;$j++) {
    $c[$j]=$coeff2->slice($j)->sclr;
#    print "$c[$j]\n";
}
my @s_x_out2;
for ($j=0;$j<$n;$j++) {
    $s_y_out2[$j]=0;
    for ($k=0;$k<$npoly;$k++) {
	$s_y_out2[$j]=$s_y_out2[$j]+$c[$k]*($x[$j]**$k);
    }
}



pgbeg(0,$dev,1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgsci(7);
pgline($n,\@x,\@a_fix);
pgsci(1);
pgpoint($n,\@x,\@a_fix,2);
pgsci(3);
pgpoint($nn,\@xx,\@yy,3);
#pgsci(15);
pgsci(2);
pgline($n,\@x,\@s_y_out2);

#pgline($n,\@x,\@s_y_out2);
pgsci(1);
pgclos();
pgend();


open(OUTPUT,">$out_file");
for ($k=0;$k<$n;$k++) {
    print OUTPUT "$n_mod\n";
    for ($i=0;$i<$n_mod;$i++) {
	print OUTPUT "$model[$i] ";
	for ($j=0;$j<9;$j++) {
	    if (($j==$n_param)&&($i==$n_model)) {
		$a[$i][$j][$k]=$s_y_out2[$k];
#		$a[$i][$j][$k]=$s_y_out[$k];
	    }
	    print OUTPUT "$a[$i][$j][$k] $ia[$i][$j][$k] ";
	}
	print OUTPUT "\n";
    }
};
close(OUTPUT);

pgbeg(0,'poly_force.ps/CPS',1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgsci(7);
pgline($n,\@x,\@a_fix);
pgsci(1);
pgpoint($n,\@x,\@a_fix,2);
#pgsci(15);
pgsci(2);
pgline($n,\@x,\@s_y_out);

#pgline($n,\@x,\@s_y_out2);
pgsci(1);
pgclos();
pgend();

exit;
#




