#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

#use Statistics::OLS;
#use Math::FFT;
#use Math::Stat;
#use Math::Spline qw(spline linsearch binsearch);
#use Math::Derivative qw(Derivative2);
#use Math::Approx;
#use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<6) {
    print "USE: extract_spec_g3d.pl INPUT_LOGFILE QUALITY.txt OUTPUT_FILE model_number mag_zeropoint CRVAL CDELT\n";
    exit;
}

$input_log=$ARGV[0];
$qfile=$ARGV[1];
$output_file=$ARGV[2];
$model_number=$ARGV[3];
$mag_zero=$ARGV[4];
$crval=$ARGV[5];
$cdelt=$ARGV[6];

if ($qfile ne "null") {
$nq=0;
open(FH,"<$qfile");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $q_id[$nq]=$data[0];
    $q_val[$nq]=$data[1];
    $nq++;
}
close(FH);
} else {
    $nq=10000;
    for ($i=0;$i<$nq;$i++) {
    	$q_val[$i]=1;
    }
}


$n=0;
open(FH,"<$input_log");
open(OUT,">$output_file");
#while ($line=<FH>) {
for ($i=0;$i<$nq;$i++) {
    if ($q_val[$i]==1) {
	do {
	    $line=<FH>;
	    if (eof(FH)) {
		$nq=$n;
	    }
	} while (($line !~ "Output")&&($n != $nq));
#while ($line !~ "Output");

#	if ($line =~ "Output") {
	    $line=<FH>; #blank line
	    for ($j=0;$j<$model_number;$j++) {
		$results=<FH>;
		$errors=<FH>;
	    }
	    $results =~ s/\://;
	    $errors =~ s/\://;
	    $results =~ s/\,//;
	    $errors =~ s/\,//;
	    $results =~ s/\(//;
	    $errors =~ s/\(//;
	    $results =~ s/\)//;
	    $errors =~ s/\)//;
	    @data=split(" ",$results);
	    $model_name[$i]=$data[0];
	    $x[$i]=$data[1];
	    $y[$i]=$data[2];
	    $mag[$i]=$data[3];
	    $re[$i]=$data[4];
	    $nsersic[$i]=$data[5];
	    $ab[$i]=$data[6];
	    $pa[$i]=$data[7];
	    @data=split(" ",$errors);	
	    $e_x[$i]=$data[0];
	    $e_y[$i]=$data[1];
	    $e_mag[$i]=$data[2];
	    $e_re[$i]=$data[3];
	    $e_nsersic[$i]=$data[4];
	    $e_ab[$i]=$data[5];
	    $e_pa[$i]=$data[6];
	    $flux[$i]=10**(0.4*($mag_zero-$mag[$i]));
#	    $e_flux[$i]=0.4*$e_mag[$i]/(0.1+$mag[$i]);
	    $wave=$crval+$cdelt*$i;
	    $k=$n+1;
	print OUT "$i $wave $flux[$i]\n";
#	print "$k $wave $flux[$i] $mag[$i]\n";
	$n++;
#	}
    } else {
	$wave=$crval+$cdelt*$i;
	print OUT "$i $wave 0\n";
    }
}
close(OUT);
close(FH);

print "DONE\n";
exit;
