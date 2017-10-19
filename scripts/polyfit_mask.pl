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


if ($#ARGV<3) {
    print "USE: poly.pl INPUT1.FITS NPOLY OUTPUT.FITS MASK_FILE [NX1 NX2] [DEV]\n";
    print "mask_list: Start_wave End_wave [to remove]\n";
    exit;
}

$infile=$ARGV[0];
$npoly=$ARGV[1];
$outfile=$ARGV[2];
$mask_file=$ARGV[3];
$dev="/xs";
$n_mask=0;
open(FH,"<$mask_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $start_mask[$n_mask]=$data[0];
    $end_mask[$n_mask]=$data[1];
    $n_mask++;
}
close(FH);
print "$n_mask spectral regions to mask\n";


$a_in=rfits($infile);
$h=$a_in->gethdr;
$crval=$h->{CRVAL1};
$cdelt=$h->{CDELT1};
if ($cdelt==0) {
    $cdelt=1;
}
($nx,$ny)=$a_in->dims;
$nx1=0;
$nx2=0;
if ($#ARGV>3) {
    $nx1=$ARGV[4];
    $nx2=$ARGV[5];
    $dev=$ARGV[6];
} else {
    $nx2=$nx;
}
 
$a_out=zeroes($nx,$ny);
for ($j=0;$j<$ny;$j++) {
    my @a_spec;
    my $w_out;
    $nt=0;
    $min=1e12;
    $max=-1e12;

    for ($i=0;$i<$nx;$i++) {
	$w_out[$i]=$crval+$cdelt*$i;
	if (($i>$nx1)&&($i<$nx2)) {
	    $copy=1;
	    for ($k=0;$k<$n_mask;$k++) {
		if (($w_out[$i]>$start_mask[$k])&&($w_out[$i]<$end_mask[$k])) {
		    $copy=0;
		}
	    }
	    if ($copy==1) {
		$w_in[$nt]=$crval+$cdelt*$i;
		$a_spec[$nt]=$a_in->at($i,$j);
		if ($min>$a_spec[$nt]) {
		    $min=$a_spec[$nt];
		}
		if ($max<$a_spec[$nt]) {
		    $max=$a_spec[$nt];
		}	    
#		print "$w_in[$nt] $a_spec[$nt]\n";
		$nt++;
	    }
	}
    }
    ($s_x,$coeff) = fitpoly1d(pdl(@w_in),pdl(@a_spec),$npoly);
    @c=list($coeff);
    for ($i=0;$i<$nx;$i++) {
	$out_spec[$i]=0;
	for ($k=0;$k<$npoly;$k++) {
	    $out_spec[$i]=$out_spec[$i]+$c[$k]*(($crval+$cdelt*$i)**$k);
#	    print "$k $c[$k] $out_spec[$i]\n";
	}
#	if ($min>$out_spec[$i]) {
#       $min=$out_spec[$i];
#	}
#	if ($max<$out_spec[$i]) {
#	    $max=$out_spec[$i];
#	}

	set($a_out,$i,$j,$out_spec[$i]);
    }

    pgbegin(0,$dev,1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($crval,$crval+$cdelt*$nx,$min,$max,0,0);
    pglabel("X-axis","Y-Axis","$j/$ny");
    pgsci(1);
#    pgpoint($nt,\@w_in,\@a_spec,5);
    pgslw(2);
    pgline($nt,\@w_in,\@a_spec);
    pgsci(3);
    pgslw(1);
    pgline($nx,\@w_out,\@out_spec);
    pgsci(1);
    pgsci(1);
    pgclose;
    
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
