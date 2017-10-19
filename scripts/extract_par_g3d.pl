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

if ($#ARGV<5) {
    print "USE: extract_par_g3d.pl INPUT_LOGFILE QUALITY.txt OUTPUT_FILE model_number par_number DEVICE [min max]\n";
    exit;
}

$input_log=$ARGV[0];
$qfile=$ARGV[1];
$output_file=$ARGV[2];
$model_number=$ARGV[3];
$par_number=$ARGV[4];
$dev=$ARGV[5];
$def=0;
if ($#ARGV==7) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $def=1;
}

$y_min=1e12;
$y_max=-1e12;

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
#else {
#    $nq=    
#}

$n=0;
$i=0;
open(FH,"<$input_log");
open(OUT,">$output_file");
for ($i=0;$i<$nq;$i++) {
#    print "$i/$nq\n";
#    if ($qfile eq "null") {
#	$q_val[$i]=1;
#    }
    if ($q_val[$i]>0) {
#	$line=<FH>;
#	print "$line\n";
#while ($line=<FH>) {
#	 ($line =~ "Output") {
	do {
	    $line=<FH>;	    
	    if (eof(FH)) {
		$nq=$n;
	    }
#	    print "$i $n $nq\n";
	} while (($line !~ "Output")&&($n != $nq));
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
	    $wave[$i]=$i;
	    $par[$i]=$data[$par_number];
	    $k=$n+1;
#	    print "PAR=$i $par[$i]\n";
	    print OUT "$i $i $par[$i]\n";
	    if ($par[$i]>$y_max) {
		$y_max=$par[$i];
	    }
	    if ($par[$i]<$y_min) {
		$y_min=$par[$i];
	    }
#	print "$k $wave $flux[$n] $mag[$n]\n";
	    $n++;
#	}
    } else {
	$model_name[$i]="none";
	$wave[$i]=$i;
	$par[$i]=0;
	$k=$n+1;
	print OUT "$i $i $par[$i]\n";
    }
}
    close(OUT);
close(FH);

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
} else {
    $y_min=$y_min-0.1;    
    $y_max=$y_max+0.1;
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wave[0],$wave[$n-1],$y_min,$y_max,0,0);
pglabel("Index","Par $par_number","");
pgline($n,\@wave,\@par);
pgclose;
pgend;



print "DONE\n";
exit;
