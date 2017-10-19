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
    print "USE: map_auto_ssp.pl elines_OUT PREFIX_OUT\n";
    exit;
}

$infile=$ARGV[0];
$sspfile=$ARGV[0];
$cut="elines_";
$sspfile =~ s/$cut//;
$prefix=$ARGV[1];

$call="map_auto_ssp_AGE_MET.pl ".$sspfile." ".$prefix;
system($call);

$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    ($nx,$ny,$ix,$iy)=split(" ",$line);
    $line=<FH>;
    chop($line);
    $nmod=$line;
#    print "$nmod\n";
    if ($n==0) {
	$pdl_flux=zeroes($nx,$ny,$nmod);
	$pdl_eflux=zeroes($nx,$ny,$nmod);
	$pdl_vel=zeroes($nx,$ny,$nmod);
	$pdl_disp=zeroes($nx,$ny,$nmod);
	$nmod_max=$nmod;
    }
    for ($i=0;$i<$nmod;$i++) {
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	#print "$ix $iy $nx,$ny $i,$nmod,$nmod_max $data[1] | $wave_name[$i] $wave[$i]\n";
	set($pdl_flux,$ix,$iy,$i,$data[3]);
	set($pdl_eflux,$ix,$iy,$i,$data[4]);
	set($pdl_disp,$ix,$iy,$i,$data[5]);
	set($pdl_vel,$ix,$iy,$i,$data[7]);	
	if ($nmod==$nmod_max) {
	    $wave[$i]=$data[1];
	    ($w,$dw)=split(/\./,$wave[$i]);
	    $wave_name[$i]=$w;
	}

    }
    $n++;
}
close(FH);
#exit;
#print "NMOD=$nmod_max\n";
for ($i=0;$i<$nmod_max;$i++) {
    $ii=$i;
    if ($ii<10) {
	$ii="0".$i;
    }
    if ($wave_name[$i]>0) {
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
    }
}

exit;





