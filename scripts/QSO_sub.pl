#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#



use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
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

if ($#ARGV<8) {
    print "USE: QSO_sub.pl input.fits X Y output.fits model.fits We1 We2 Wc1 Wc2\n";
    exit;
}
$infile=$ARGV[0];
$xc=$ARGV[1];
$yc=$ARGV[2];
$output=$ARGV[3];
$model=$ARGV[4];
$we1=$ARGV[5];
$we2=$ARGV[6];
$wc1=$ARGV[7];
$wc2=$ARGV[8];



$in_pdl=rfits($infile);
$h=$in_pdl->gethdr;
($nx,$ny,$nz)=$in_pdl->dims;
$crpix1=$h->{CRPIX1};
$crval1=$h->{CRVAL1};
$cdelt1=$h->{CDELT1};
$crpix2=$h->{CRPIX2};
$crval2=$h->{CRVAL2};
$cdelt2=$h->{CDELT2};
$crpix3=$h->{CRPIX3};
$crval3=$h->{CRVAL3};
$cdelt3=$h->{CDELT3};

$slice_pdl=zeroes($nx,$ny);
$k=500;
$ip=0;
$jp=0;
$peak=-1e12;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$in_pdl->at($i,$j,$k);
#	set($slice_pdl,$i,$j,$val);
	if ($val>$peak) {
	    $peak=$val;
	    $ip=$i;
	    $jp=$j;	    
	}
    }
}
print "PEAK at $ip,$jp\n";

$ip=int($xc);
$jp=int($yc);


$spec2D = sumover $in_pdl;
$spectrum = sumover $spec2D;
#$spectrum = $spectrum/($nx*$ny);
@spec_sum=list($spectrum);
#print "@spec_sum\n";

$slice_peak=$in_pdl->slice("$ip,$jp,:");
@spec_peak=list($slice_peak);
($min,$max)=minmax(@spec_sum);

$c_sum=0;
$e_sum=0;
$n_c=0;
$n_e=0;
$c_peak=0;
$e_peak=0;
for ($k=0;$k<$nz;$k++) {
    $wave[$k]=$crval3+($k-($crpix3-1))*$cdelt3;
#    print "$wave[$k] $spec_peak[$k] $spec_sum[$k]\n";
    if (($wave[$k]>$we1)&&($wave[$k]<$we2)) {
	$e_sum=$e_sum+$spec_sum[$k];
	$e_peak=$e_peak+$spec_peak[$k];
	$n_e++;
    }
    if (($wave[$k]>$wc1)&&($wave[$k]<$wc2)) {
	$c_sum=$c_sum+$spec_sum[$k];
	$c_peak=$c_peak+$spec_peak[$k];
	$n_c++;
    }
}
$e_sum=$e_sum/$n_e;
$e_peak=$e_peak/$n_e;
$c_sum=$c_sum/$n_c;
$c_peak=$c_peak/$n_c;
$ratio=($e_sum-$c_sum)/($e_peak-$c_peak);
for ($k=0;$k<$nz;$k++) {
    $spec_scaled[$k]=$spec_peak[$k]*$ratio;
}

pgbegin(0,"/xs",1,1);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wave[0],$wave[$nz-1],0,$max,0,0);
pglabel("Wavelength","Flux","QSO 0-spec");
pgsci(1);
pgline($nz,\@wave,\@spec_peak);
pgsci(2);
pgline($nz,\@wave,\@spec_sum);
pgsci(3);
pgline($nz,\@wave,\@spec_scaled);
pgclos();
pgend();

#
# Scale everything!
#
$scaled_pdl=zeroes($nx,$ny,$nz);
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$slice_ij=$in_pdl->slice("$i,$j,:");
	@spec_ij=list($slice_ij);
	$c_ij=0;
	$e_ij=0;
	$n_c=0;
	$n_e=0;
	$c_peak=0;
	$e_peak=0;
	for ($k=0;$k<$nz;$k++) {
#	    $wave[$k]=$crval3+($k-($crpix3-1))*$cdelt3;
	    if (($wave[$k]>$we1)&&($wave[$k]<$we2)) {
		$e_ij=$e_ij+$spec_ij[$k];
		$e_peak=$e_peak+$spec_peak[$k];
		$n_e++;
	    }
	    if (($wave[$k]>$wc1)&&($wave[$k]<$wc2)) {
		$c_ij=$c_ij+$spec_ij[$k];
		$c_peak=$c_peak+$spec_peak[$k];
		$n_c++;
	    }
	}
	$e_ij=$e_ij/$n_e;
	$e_peak=$e_peak/$n_e;
	$c_ij=$c_ij/$n_c;
	$c_peak=$c_peak/$n_c;
	$ratio=($e_ij-$c_ij)/($e_peak-$c_peak);
	for ($k=0;$k<$nz;$k++) {
	    $spec_scaled[$k]=$spec_peak[$k]*$ratio;
	    set($scaled_pdl,$i,$j,$k,$spec_scaled[$k]);
	}
    }
}
$res_pdl=$in_pdl-$scaled_pdl;
$scaled_pdl->sethdr($h);
$res_pdl->sethdr($h);
$scaled_pdl->wfits($model);
$res_pdl->wfits($output);

($min,$max)=minmax(@spec_peak);
pgbegin(0,"/xs",1,1);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wave[0],$wave[$nz-1],0,$max,0,0);
pglabel("Wavelength","Flux","QSO 0-spec");
$color=1;
#
# We iterate
#
$niter=0;
$delta_peak=1000;
while ($delta_peak>0.0005) {
    #
    # We create the EELR of the center
    # 
    pgsci($color);
    pgline($nz,\@wave,\@spec_peak);
    ($min_old,$max_old)=minmax(@spec_peak);
   $color++;
    $slice_ring=zeroes($nz);
    @DIM=$slice_ring->dims();
#    print "DIM=@DIM\n";
    $nr=0;
    for ($i=$ip-1;$i<$ip+2;$i++) {
	for ($j=$jp-1;$j<$jp+2;$j++) {
#	    $slice_peak=$in_pdl->slice("$ip,$jp,:");    
	    if (($i!=$ip)||($j!=$jp)) {
		for ($k=0;$k<$nz;$k++) {
		    $val=$slice_ring->at($k);
		    $val2=$res_pdl->at($i,$j,$k);    
#		    print "$val $val2\n";
		    $val=$val+$val2;
		    set($slice_ring,$k,$val);
		}
		$nr++;
	    }
	}
    }
#    print "NR=$nr\n";
    $slice_ring=$slice_ring/$nr;
    @DIM=$slice_ring->dims();
#    print "DIM=@DIM\n";
    $slice_peak=$in_pdl->slice("$ip,$jp,:");
#    $slice_peak=$slice_peak-$slice_ring;
    @spec_peak=list($slice_peak);
    for ($k=0;$k<$nz;$k++) {
	$val=$spec_peak[$k];
	$val2=$slice_ring->at($k);
	$spec_peak[$k]=$spec_peak[$k]-$val2;
#	print "$k $wave[$k] $val $val2\n";
    }
    ($min_new,$max_new)=minmax(@spec_peak);
    $delta_peak=abs($max_new-$max_old)/$max;
#    print "@spec_peak\n";
#
# Scale everything!
#
    $scaled_pdl=zeroes($nx,$ny,$nz);
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
	    $slice_ij=$in_pdl->slice("$i,$j,:");
	    @spec_ij=list($slice_ij);
	    $c_ij=0;
	    $e_ij=0;
	    $n_c=0;
	    $n_e=0;
	    $c_peak=0;
	    $e_peak=0;
	    for ($k=0;$k<$nz;$k++) {
#	    $wave[$k]=$crval3+($k-($crpix3-1))*$cdelt3;
		if (($wave[$k]>$we1)&&($wave[$k]<$we2)) {
		    $e_ij=$e_ij+$spec_ij[$k];
		    $e_peak=$e_peak+$spec_peak[$k];
		    $n_e++;
		}
		if (($wave[$k]>$wc1)&&($wave[$k]<$wc2)) {
		    $c_ij=$c_ij+$spec_ij[$k];
		    $c_peak=$c_peak+$spec_peak[$k];
		    $n_c++;
		}
	    }
	    $e_ij=$e_ij/$n_e;
	    $e_peak=$e_peak/$n_e;
	    $c_ij=$c_ij/$n_c;
	    $c_peak=$c_peak/$n_c;
	    if ($e_peak!=$c_peak) {
		$ratio=($e_ij-$c_ij)/($e_peak-$c_peak);
	    } else {
		$ratio=0;
	    }
#	    print "$niter ($i,$j) $n_e $n_c $e_ij $c_ij $e_peak $c_peak $ratio\n";

	    for ($k=0;$k<$nz;$k++) {
		$spec_scaled[$k]=$spec_peak[$k]*$ratio;
		set($scaled_pdl,$i,$j,$k,$spec_scaled[$k]);
	    }
	}
    }
    $res_pdl=$in_pdl-$scaled_pdl;
    $scaled_pdl->sethdr($h);
    $res_pdl->sethdr($h);
    $scaled_pdl->wfits("scaled_pdl_I.fits");
    $res_pdl->wfits("res_pdl_I.fits");



    print "$niter $delta_peak\n";
    $niter++;
}
pgclos();
pgend();

exit;

