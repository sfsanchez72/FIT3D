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


if ($#ARGV<2) {
    print "USE: map_auto_ssp.seg.pl elines_OUT seg_file.fits PREFIX_OUT\n";
    exit;
}

$infile=$ARGV[0];
$sspfile=$ARGV[0];
$cut="elines_";
$sspfile =~ s/$cut//;
$seg_file=$ARGV[1];
$prefix=$ARGV[2];

$pdl_seg=rfits($seg_file);
($nx,$ny)=$pdl_seg->dims;



$call="map_auto_ssp_AGE_MET.seg.pl ".$sspfile." ".$seg_file." ".$prefix;
system($call);
print "Reading fit output file\n";
$n=0;
$l=0;
$nmod=0;
$N=0;
$first=0.0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line =~ "eline") {
	@data=split(" ",$line);
	if (($first>0)&&($first==$data[1])) {
	    $n++;
	    $l=0;
	}
	if (($n==0)&&($l==0)) {
	    $first=$data[1];
	}
	$wave[$l]=$data[1];
	($w,$dw)=split(/\./,$wave[$l]);
	$wave_name[$l]=$w;
	$flux[$n][$l]=$data[3];
	$eflux[$n][$l]=$data[4];
	$sig[$n][$l]=$data[5];
	$esig[$n][$l]=$data[6];
	$vel[$n][$l]=$data[7];
	$evel[$n][$l]=$data[8];
	if ($l>$L) {
	    $nmod=$l;
	}
	$l++;
    }
}
close(FH);

print "DONE\n";
print "Feeding arrays\n";
$pdl_flux=zeroes($nx,$ny,$nmod);
$pdl_eflux=zeroes($nx,$ny,$nmod);
$pdl_vel=zeroes($nx,$ny,$nmod);
$pdl_disp=zeroes($nx,$ny,$nmod);
$N=$n+1;
for ($ix=0;$ix<$nx;$ix++) {
    for ($iy=0;$iy<$ny;$iy++) {
	$iseg=$pdl_seg->at($ix,$iy);
	if ($iseg>0) {
	    $iseg=$iseg-1;
	    for ($i=0;$i<$nmod;$i++) {
		set($pdl_flux,$ix,$iy,$i,$flux[$iseg][$i]);
		set($pdl_eflux,$ix,$iy,$i,$eflux[$iseg][$i]);
		set($pdl_disp,$ix,$iy,$i,$sig[$iseg][$i]);
		set($pdl_vel,$ix,$iy,$i,$vel[$iseg][$i]);
	    }
	}
    }
}
print "Done\n";
print "Writing output\n";

for ($i=0;$i<$nmod;$i++) {
    $ii=$i;
    if ($ii<10) {
	$ii="0".$i;
    }
    $flux_name=$prefix."_flux_".$wave_name[$i].".fits";
    $eflux_name=$prefix."_eflux_".$wave_name[$i].".fits";
    $disp_name=$prefix."_disp_".$wave_name[$i].".fits";
    $vel_name=$prefix."_vel_".$wave_name[$i].".fits";
    my $map_flux=$pdl_flux->slice(",,($i)");
    my $map_eflux=$pdl_eflux->slice(",,($i)");
    my $map_vel=$pdl_vel->slice(",,($i)");
    my $map_disp=$pdl_disp->slice(",,($i)");
    #print "$vel_name\n";
    $map_flux->wfits($flux_name);
    $map_eflux->wfits($eflux_name);
    $map_vel->wfits($vel_name);
    $map_disp->wfits($disp_name);
    $call="write_img_header.pl ".$flux_name." CRVAL1 ".$wave[$i];
    system($call);
    $call="write_img_header.pl ".$eflux_name." CRVAL1 ".$wave[$i];
    system($call);
    $call="write_img_header.pl ".$vel_name." CRVAL1 ".$wave[$i];
    system($call);
    $call="write_img_header.pl ".$disp_name." CRVAL1 ".$wave[$i];
    system($call);
    print "$wave_now[$i] DONE\n";
}
print "DONE\n";

exit;





