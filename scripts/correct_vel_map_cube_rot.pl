#!/usr/bin/perl
#
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
#use PDL::Graphics::TriD;
#use PDL::Graphics::TriD::Image;
use PDL::Fit::Gaussian;
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;

use PDL::Func;
use PDL::Math;
use PDL::Primitive;
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


$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<2) {
    print "USE: correct_vel_map_cube_rot.pl INCUBE.fits VEL_MAP.fits OUTCUBE.fits [FINAL_VEL]\n";
    print "It will use images coordinates: 0,nx 0,ny\n";
    exit;
}

$input=$ARGV[0];
$velmap=$ARGV[1];
$output=$ARGV[2];
$final_vel=0;
if ($#ARGV==3) {
    $final_vel=$ARGV[3];
}
$pdl_cube=rfits("$input");
($NX,$NY,$NZ)=$pdl_cube->dims();
$out_pdl=$pdl_cube;
$h=$pdl_cube->gethdr;
$crval=$pdl_cube->hdr->{"CRVAL3"};
$cdelt=$pdl_cube->hdr->{"CDELT3"};
$crpix=$pdl_cube->hdr->{"CRPIX3"};
for ($i=0;$i<$NZ;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
#    print "$i $wave[$i]\n";
}
#exit;
$pdl_wave=pdl(@wave);
$wmin=$wave[0];
$wmax=$wave[$NZ-1];

$pdl_vel_map=rfits($velmap);
($nx,$ny)=$pdl_vel_map->dims;

if (($nx!=$NX)||($ny!=$NY)) {
    print "Cube Dimensions ($NX,$NY,$NZ), do not match with vel-map dimensions ($nx,$ny)\n";
    exit;
}

for ($i=0;$i<$NX;$i++) {
    for ($j=0;$j<$NY;$j++) {
	my $pdl_spec=$pdl_cube->slice("($i),($j),");
	my $pdl_spec_out=$out_pdl->slice("($i),($j),");	
	my $pdl_new_wave=zeroes($NZ);;
	$vel=$pdl_vel_map->at($i,$j);
	$vel=$vel-$final_vel;
	for ($k=0;$k<$NZ;$k++) {	    
	    $val=$wave[$k]*(1+$vel/300000);
	    set($pdl_new_wave,$k,$val);
	} 
	$pdl_spec->inplace->setbadtoval(0); 
	$out_spec_pdl=interpol($pdl_new_wave,$pdl_wave, $pdl_spec);
	$pdl_spec_out .= $out_spec_pdl;
#	print "$i/$NX; $j/$NY\n";
    }
	print ".";
}

$out_pdl->wfits("$output");

exit;
