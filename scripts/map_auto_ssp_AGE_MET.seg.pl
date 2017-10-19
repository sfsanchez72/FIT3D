#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
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
use PDL::Image2D;
use  PDL::Fit::Linfit;



use PDL::Core;
use PDL::Basic;
use PDL::Exporter;
@ISA    = qw( PDL::Exporter );
use PDL::Options ':Func';
use PDL::Slatec; # For matinv()


$vel_light=299792.458;


if ($#ARGV<1) {
    print "USE: map_auto_ssp.pl AUTO.OUTseg_file.fits  PREFIX_OUT\n";
    exit;
}

$infile=$ARGV[0];
$seg_file=$ARGV[1];
$prefix=$ARGV[2];

$pdl_seg=rfits($seg_file);
($nx,$ny)=$pdl_seg->dims;

$pdl_chi=zeroes($nx,$ny);
$pdl_age=zeroes($nx,$ny);
$pdl_met=zeroes($nx,$ny);
$pdl_Av=zeroes($nx,$ny);
$pdl_disp=zeroes($nx,$ny);
$pdl_vel=zeroes($nx,$ny);
$pdl_flux=zeroes($nx,$ny);


print "Reading input files\n";
$nc=-1;
$nc_max=-1;
$coeffs_infile="coeffs_".$infile;
open(FH,"<$coeffs_infile");
while ($line=<FH>) {
    chop($line);
    if ($line !~ "ID") {
	@data=split(" ",$line);
	$i=$data[0];
	if ($i==0) {
	    $nc++;
	}
	$AGE[$i][$nc]=$data[1];
	$MET[$i][$nc]=$data[2];
	$COEFF[$i][$nc]=$data[3];
	$NORM[$i][$nc]=$data[4];
	$ML[$i][$nc]=$data[5];
	$AV[$i][$nc]=$data[6];
	if ($nc_max<=$i) {
	    $nc_max=$i+1;
	}
    }
}
close(FH);
#print "$i $nc $nc_max $nx,$ny\n";
$pdl_COEFF=zeroes($nx,$ny,$nc_max);
$pdl_NORM=zeroes($nx,$ny,$nc_max);	
$pdl_AV=zeroes($nx,$ny,$nc_max);	
print "DONE\n";

$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$chi[$n]=$data[0];
	$age[$n]=$data[1];
	$met[$n]=$data[2];
	$Av[$n]=$data[3];
	$vel[$n]=$data[4];
	$disp[$n]=$data[5];
	$flux[$n]=$data[8];
	$n++;
    }
}
close(FH);

print "Feeding the arrays\n";
for ($ix=0;$ix<$nx;$ix++) {
    for ($iy=0;$iy<$ny;$iy++) {
	$iseg=$pdl_seg->at($ix,$iy);
	if ($iseg>0) {
	    $iseg=$iseg-1;
	    set($pdl_chi,$ix,$iy,$chi[$iseg]);
	    set($pdl_age,$ix,$iy,$age[$iseg]);
	    set($pdl_met,$ix,$iy,$met[$iseg]);
	    set($pdl_Av,$ix,$iy,$Av[$iseg]);
	    set($pdl_vel,$ix,$iy,$vel[$iseg]);
	    set($pdl_disp,$ix,$iy,$disp[$iseg]);
	    set($pdl_flux,$ix,$iy,$flux[$iseg]);
	    for ($i=0;$i<$nc_max;$i++) {
	#	print "$ix,$iy [$nx,$ny] $i [$nc_max]\n";
		set($pdl_NORM,$ix,$iy,$i,$NORM[$i][$iseg]);
	    }
	}
    }
}
print "DONE\n";
print "Writting the output files\n";


$chi_name=$prefix."_chi_ssp.fits";
$pdl_chi->wfits($chi_name);

$age_name=$prefix."_age_ssp.fits";
$pdl_age->wfits($age_name);

$met_name=$prefix."_met_ssp.fits";
$pdl_met->wfits($met_name);

$Av_name=$prefix."_Av_ssp.fits";
$pdl_Av->wfits($Av_name);

$disp_name=$prefix."_disp_ssp.fits";
$pdl_disp->wfits($disp_name);

$pdl_vel=$pdl_vel*300000;
$vel_name=$prefix."_vel_ssp.fits";
$pdl_vel->wfits($vel_name);

$flux_name=$prefix."_flux_ssp.fits";
$pdl_flux->wfits($flux_name);

for ($j=0;$j<$nc_max;$j++) {
    $file_name=$prefix."_NORM_".$j."_ssp.fits";
    $sec=$pdl_NORM->slice(":,:,$j");
    $sec->wfits($file_name);
    $call="write_img_header.pl ".$file_name." AGE ".$AGE[$j][0];
    system($call);
    $call="write_img_header.pl ".$file_name." MET ".$MET[$j][0];
    system($call);

}

exit;





