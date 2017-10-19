#!/usr/bin/perl
#
use Statistics::OLS;
use Math::FFT;
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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<1) {
    print "USE: imcombine.pl images.list output.fits\n";
    exit;
}

$config=$ARGV[0];
$outfile=$ARGV[1];

$n=0;
open(FH,"<$config");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $file[$n]=$data[0];
    $n++;
}

print "$n files to combine\n";

for ($i=0;$i<$n;$i++) {
    my $pdl_in=rfits($file[$i]);
    ($nx,$ny)=$pdl_in->dims();
    $h=$pdl_in->gethdr;
    if ($i==0) {
	$pdl_tmp=zeroes($nx,$ny,$n);	
       
    }
    my $t = $pdl_tmp->slice(":,:,($i)");
    $t .= $pdl_in;    
}
print "Files readed\n";

$pdl_out=medover($pdl_tmp->xchg(0,2));#$xchg(0,2);#zeroes($nx,$ny);
$pdl_out=$pdl_out->xchg(0,1);
$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);


exit;


for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	my $t = $pdl_tmp->slice("($i),($j),:");
	($mean,$rms,$median,$min,$max) = stats($t);
	$median=$med;
	set($pdl_out,$i,$j,$med);
    }
}

	my @a=0;
	my $na=0;
	for ($k=0;$k<$n;$k++) {
	    $val=$t->at($k);
	    if (abs($val-$median)<0.7*$median) {
		$a[$na]=$val;
		$na++;
	    } 
	}
	if ($na<2) {
	    $med=$median;
	} else {
	    $med=mean(@a);
	}
