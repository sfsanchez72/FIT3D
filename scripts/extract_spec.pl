#!/usr/bin/perl
#
#
#
use PGPLOT;


use Math::SpecFun::Gamma qw(gamma gammaln gammainc gammaincc);
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;



# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fix one parameter to do
# a force fitting latter
if ($#ARGV<3) {
    print "USE: extract_spec.pl start_w delta_w INPUT_FILE OUT_FILE\n";        
    exit;
}

$start_w=$ARGV[0];
$delta_w=$ARGV[1];
$input_file=$ARGV[2];
$out_file=$ARGV[3];


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
	}
	#
	# We detemine the flux depending of the model
	#
	$norm=0;
	if ($model[$i] eq "gauss") {
	    $norm=1;
	}
	if ($model[$i] eq "moffat") {
	    $norm=1;
	}
	if ($model[$i] eq "sersic") {
	    $norm=gamma((2/$a[$i][4][$n]))/$a[$i][4][$n];
	}
	if ($model[$i] eq "gsersic") {
	    $norm=gamma((2/$a[$i][4][$n]))/$a[$i][4][$n];
	    $norm=$norm/10;
	}
	if ($model[$i] eq "disk") {
	    $norm=1;
	}
	if ($model[$i] eq "vauc") {
	    $norm=43200;
	}
	
	$flux[$i][$n]=$a[$i][2][$n]*2*3.14159*($a[$i][3][$n]**2)*$norm;
	$e_flux[$i][$n]=($ia[$i][2][$n]*($a[$i][3][$n]**2)+2*$a[$i][2][$n]*$a[$i][3][$n])*2*3.14159**$norm;
	if ($y_max<$flux[$i][$n]) {
	    $y_max=$flux[$i][$n];
	}
	if ($y_min>$flux[$i][$n]) {
	    $y_min=$flux[$i][$n];
	}
	
    }
    $x[$n]=$start_w+$delta_w*$n;


    $n++;
};
close(INPUT);

#
# We fit to a polynomical function
# of the order npoly (using PDL);
#


open(OUTPUT,">$out_file");
print OUTPUT "# $n_mod\n";
print OUTPUT "# (1) wavelength \n";
for ($i=0;$i<$n_mod;$i++) {
    $j=$i+2;
    print OUTPUT "# ($j) $model[$i]  flux\n";
    $j=$i+3;
    print OUTPUT "# ($j) $model[$i] Errorflux\n";
}
for ($k=0;$k<$n;$k++) {
    print OUTPUT "$x[$k] ";
    for ($i=0;$i<$n_mod;$i++) {
	print OUTPUT "$flux[$i][$k] $e_flux[$i][$k] ";	
    }
    print OUTPUT "\n";
}
close(OUTPUT);

pgbeg(0,'/XS',1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgsci(7);
for ($i=0;$i<$n_mod;$i++) {
    pgsci($i+1);
    for ($k=0;$k<$n;$k++) {
	$a_fix[$k]=$flux[$i][$k];
    }
    pgline($n,\@x,\@a_fix);
}

pgsci(1);
pgclos();
pgend();


pgbeg(0,'extract_spec.ps/CPS',1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgsci(7);
for ($i=0;$i<$n_mod;$i++) {
    pgsci($i+1);
    for ($k=0;$k<$n;$k++) {
	$a_fix[$k]=$flux[$i][$k];
    }
    pgline($n,\@x,\@a_fix);
}

pgsci(1);
pgclos();
pgend();

exit;
#




