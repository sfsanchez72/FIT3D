#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


#
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<5) {
    print "USE: imcntr_cube.pl input_cube.fits X Y NP_X NP_Y centers_list.txt [N_START] [N_END]\n";
    exit;
}

$input_cube=$ARGV[0];
$x0=$ARGV[1];
$y0=$ARGV[2];
$npoly_x=$ARGV[3];
$npoly_y=$ARGV[4];
$centers_list=$ARGV[5];
$ps="cn_".$input_cube;
$ps =~ s/.fits/.ps/;
$dev=$ps."/CPS";
$n_start=0;
$n_end=0;
if ($#ARGV==7) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
}

$na=read_naxes($input_cube);   
@naxis=@$na;

open (TMP,">find_centers.cl");
for ($k=0;$k<$naxis[2];$k++) {
    $kk=$k+1;
    $line="imcntr ".$input_cube."[*,*,".$kk."] ".$x0." ".$y0;
    print TMP "$line\n";
}
print TMP "logout\n";
close(TMP);

system("/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh < find_centers.cl > centers.tmp");

$n=0;
if ($x0<$y0) {
    $min=$x0-3;
    $max=$y0+3;
} else {
    $min=$y0-3;
    $max=$x0+3;
}
print "$x0 $y0 $min $max $n_start $n_end\n";

open (TMP,"<centers.tmp");
while($line=<TMP>) {
    chop($line);
    $sec=substr($line,0,3);
#    print "$sec\n";
    if (substr($line,0,3) eq "cl>" ) {
	@data=split(" ",$line);
	$id[$n]=$n+1;
	$x[$n]=$data[3];
	$y[$n]=$data[5];
	print "$n $x[$n] $y[$n]\n";
	$n++;
    }
}
close(TMP);

print "$n lines found";

if ($n_start>0) {
    for ($j=$n_start;$j>-1;$j--) {
	$x[$j]=$x[$n_start];
	$y[$j]=$y[$n_start];
    }
}
if ($n_end<$n) {
    for ($j=$n_end;$j<$n;$j++) {
	$x[$j]=$x[$n_end];
	$y[$j]=$y[$n_end];
    }
}


$y_pdl = pdl(@y);
($s_y,$coeff) = fitpoly1d $y_pdl,$npoly_y;
for ($j=0;$j<$n;$j++) {
    $y_out[$j] = $s_y->slice($j)->sclr;
}


$x_pdl = pdl(@x);
($s_x,$coeff) = fitpoly1d $x_pdl,$npoly_x;
for ($j=0;$j<$n;$j++) {
    $x_out[$j] = $s_x->slice($j)->sclr;
}

pgbegin(0,"/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
pgenv(1,$n,$min-1,$max+1,0,0);
pgsci(2);
pgpoint($n,\@id,\@x,3);
pgsci(4);
pgpoint($n,\@id,\@y,3);
pgsci(1);
pgline($n,\@id,\@x_out);
pgline($n,\@id,\@y_out);
pgclose;
pgend;

open(FH,">$centers_list");
for ($j=0;$j<$n;$j++) {
    print FH "$j $x_out[$j] $y_out[$j]\n";
}
close(FH);

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
pgenv(1,$n,$min-1,$max+1,0,0);
pgsci(2);
pgpoint($n,\@id,\@x,3);
pgsci(4);
pgpoint($n,\@id,\@y,3);
pgsci(1);
pgline($n,\@id,\@x_out);
pgline($n,\@id,\@y_out);
pgclose;
pgend;


exit;
