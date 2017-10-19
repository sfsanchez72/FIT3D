#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#f
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
$red_elines=0.0;

   

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<4) {
    print "USE: montecarlo.pl org_spec.txt res_joint_spec.txt Nsim FWHM_res OUTPUT.RSS.fits [PHOT_FILE]\n";
    print "'org_spec.txt' and 'res_joing_spec.txt' should be\n";
    print "the output of:\n";
    print "auto_ssp_elines_several_Av_log.pl\n";
    print "PHOT_FILE:\n";
    print "ID WAVE FLUX ERROR\n";
    print "FULLPATH!\n";
    print "\n";
    exit;
}

$org_file=$ARGV[0];
$res_file=$ARGV[1];
$nsim=$ARGV[2];
$fwhm=$ARGV[3];
$output_file=$ARGV[4];
$nx=0;
open(FH,"<$org_file");
open(FH2,"<$res_file");
$dpix[0]=1e12;
$cdelt=1e12;
while($line=<FH>) {
    chop($line);
    my @data=split(" ",$line);
    $id[$nx]=$data[0];
    $wave[$nx]=$data[1];
    $flux[$nx]=$data[2];
    $line2=<FH2>;
    my @data=split(" ",$line2);
    $res[$nx]=$data[2];
    $nx++;
}
close(FH2);
close(FH);
$crval=$wave[0];
$crpix=1;
$cdelt=$wave[1]-$wave[0];

#print "$crval $cdelt\n";

$pdl_flux=pdl(@flux);
@med_res=median_filter(3*$fwhm/$cdelt,\@res);
$pdl_res=pdl(@med_res);
#$pdl_res=pdl(@res);

my $pdl_output=zeroes($nx,$nsim);
$h = {NAXIS=>2, NAXIS1=>$nx, NAXIS2=>$nsim , COMMENT=>"output MONTECARLO",
      CRVAL1 => $crval,
      CRPIX1 => $crpix,
      CDELT1 => $cdelt
      }; 
$$h{FILENAME} =$output_file;
if ($#ARGV==5) {
    $out_phot="phot_".$output_file;
    $in_phot=$ARGV[5];
    $$h{PHOTFILE}=$out_phot;
}

$pdl_output->sethdr( $h );


for ($j=0;$j<$nsim;$j++) {
    my $pdl_noise=grandom($nx);
    $pdl_noise=$pdl_noise*$pdl_res;
    my $pdl_now=$pdl_flux+$pdl_noise;
    $t=$pdl_output->slice(",($j)");
    $t.=$pdl_now;
}
$pdl_output->wfits($output_file);

if ($#ARGV==5) {
    $np=0;
    open(FH,"<$in_phot");
    while($line=<FH>) {
	chop($line);
	my @data=split(" ",$line);
	$wp[$np]=$data[1];
	$fp[$np]=$data[2];
	$efp[$np]=$data[3];
	$np++;
    }
    close(FH);
    $pdl_wp=pdl(@wp);
    $pdl_fp=pdl(@fp);
    $pdl_efp=pdl(@efp);
#    print "$np $nsim $in_phot\n";
    
    $pdl_phot=zeroes($np,$nsim+1);
    for ($j=0;$j<$nsim;$j++) {
	my $pdl_noise=grandom($np);
	$pdl_noise=$pdl_noise*$pdl_efp;
	my $pdl_now=$pdl_fp+$pdl_noise;
	$t=$pdl_phot->slice(",($j)");
	$t.=$pdl_now;
#	print "$pdl_now $pdl_fp\n";
    }
    
    $t=$pdl_phot->slice(",($nsim)");
    $t.=$pdl_wp;
    $pdl_phot->wfits($out_phot);
}
	
exit;


sub add_file {
    my $J=$_[0];
    my $K=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
	my @data=split(" ",$line);
	my $val=$data[2];
	set($pdl_output,$i,$J,$K,$val);
	$i++;
    }
    close(FH);
    return;
}
