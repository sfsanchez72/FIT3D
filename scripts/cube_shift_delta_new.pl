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
use PDL::Transform;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<3) {
    print "USE: cube_shift_delta_new.pl input_cube.fits output_cube.fits delta_x delta_y\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$dx=$ARGV[2];
$dy=$ARGV[3];


$pdl_in=rfits($input_cube);
($nx,$ny,$nz)=$pdl_in->dims;
$H=$pdl_in->gethdr;
$scale=1;
$nx_0=int($nx*$scale);
$ny_0=int($ny*$scale);

print "$scale $nx_0 $ny_0 $nz\n";
$pdl_out=zeroes($nx_0,$ny_0,$nz);
$x_c=zeroes($nz);
$y_c=zeroes($nz);
$flux=zeroes($nz);
$n_start=0;
$n_end=$nz;

for ($k=0;$k<$nz-1;$k++) {
    $DX=$dx;
    $DY=$dy;
    $pdl_in_sec=$pdl_in->slice(",,($k)");
    $shift=pdl($DX,$DY);
    $tr=t_offset($shift);
    my $a=$pdl_in_sec->map($tr,{pix=>1});
    my $t = $pdl_out->slice(",,($k)");
    $t .= $a;
}

$nz_med=int($nz/2);

for ($i=0;$i<$nx_0;$i++) {
    for ($j=0;$j<$ny_0;$j++) {
	$val=$pdl_out->at($i,$j,$nz_med);
	if (abs($val)>1e100) {
	    my $t = $pdl_out->slice("($i),($j),");
	    my $zero = zeroes($nz);
	    $t .= $zero;
	}
    }
}


$CDELT1=$H->{CDELT1};
$CDELT2=$H->{CDELT2};
$CDELT3=$H->{CDELT3};
$CRPIX1=$H->{CRPIX1};
$CRPIX2=$H->{CRPIX2};
$CRPIX3=$H->{CRPIX3};
$CRVAL1=$H->{CRVAL1};
$CRVAL2=$H->{CRVAL2};
$CRVAL3=$H->{CRVAL3};

$pdl_out=$pdl_out/($scale*$scale);
$pdl_out->hdr->{CDELT1}=$CDELT1/$scale;
$pdl_out->hdr->{CDELT2}=$CDELT2/$scale;
$pdl_out->hdr->{CDELT3}=$CDELT3;
$pdl_out->hdr->{CRPIX1}=$CRPIX1;
$pdl_out->hdr->{CRPIX2}=$CRPIX2;
$pdl_out->hdr->{CRPIX3}=$CRPIX3;
$pdl_out->hdr->{CRVAL1}=$CRVAL1;
$pdl_out->hdr->{CRVAL2}=$CRVAL2;
$pdl_out->hdr->{CRVAL3}=$CRVAL3;


$pdl_out->wfits($output_cube);

exit;

