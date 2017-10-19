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




$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<10) {
    print "USE: create_sim_ssp_dust.pl BACK_LIST.fits OUT_SPEC.txt REDSHIFT SIGMA_DISP Av RMS.fits N.MODELS CRVAL CDELT NPIX FLUX_MEAN\n";
    exit;
}

$back_list=$ARGV[0];
$outfile=$ARGV[1];
$redshift=$ARGV[2];
$sigma=$ARGV[3];
$Av=$ARGV[4];
$rms=$ARGV[5];
$nmod=$ARGV[6];
$crval0=$ARGV[7];
$cdelt0=$ARGV[8];
$npix0=$ARGV[9];
$flux_mean=$ARGV[10];
$nmask=0;



$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};

for ($i=0;$i<$nf;$i++) {
    $header="NORM".$i;
    $ml[$i]=$pdl_flux_c_ini->hdr->{$header};
}

$pdl_rms=rfits($rms);
($nx_rms,$ny_rms)=$pdl_rms->dims;

$ny_now=int(rand($ny_rms));
$pdl_rms_now=$pdl_rms->slice(":,($ny_now)");


for ($j=0;$j<$npix0;$j++) {
    $wave_unc[$j]=($crval0+$cdelt0*$j);
    $wave_rest[$j]=$wave_unc[$j]/(1+$redshift);
    $dust_rat[$j]=A_l(3.1,$wave_rest[$j]);
    $dust[$j]=10**(-0.4*$Av*$dust_rat[$j]);  
}

for ($j=0;$j<$n_c;$j++) {
    $wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
#    print "$wave_c[$j] $dust[$j] $dust_rat[$j]\n";
}
$dpix_c=$wave_c[1]-$wave_c[0];
$rsigma=$sigma/$dpix_c;



#
# We create a kernel
#
$box=int(3*$rsigma);

#print "$sigma $rsigma $dpix_c $box\n";
if ($box==0) {
    $box=1;
}
$kernel=zeroes(2*$box);
$norm=0;
$flux_c[$i]=0;
for ($j=0;$j<2*$box;$j++) {
    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
    set($kernel,$j,$gaus);
    $norm=$norm+$gaus;
}
$kernel=$kernel/$norm;

$pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;

#
#
#

$sum=0;
$sum_mass=0;
for ($i=0;$i<$nmod;$i++) {
    $models[$i]=int(rand($nf));
#    $weight[$i]=1/(($i+1)**2);
    $weight[$i]=1/(($i+1));
    $w_mass[$i]=$weight[$i]*$ml[$models[$i]];
    $sum=$sum+$weight[$i];
    $sum_mass=$sum_mass+$w_mass[$i];

}


for ($i=0;$i<$nmod;$i++) {
    $weight[$i]=$weight[$i]/$sum;
    if ($sum_mass>0) {
	$w_mass[$i]=$w_mass[$i]/$sum_mass;
    }
}


$pdl_new=zeroes($npix0);

$age_now=0;
$met_now=0;
$age_mass_now=0;
$met_mass_now=0;
for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;
    ($AGE,$MET)=split("_",$name_min);
    if ($AGE =~ "Myr") {
	$age=$AGE;
	$age =~ s/Myr//;
	$age=$age/1000;
    } else {
	$age=$AGE;
	$age =~ s/Gyr//;
    }
    $met=$MET;
    $met =~ s/z/0\./;
    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
    for ($j=0;$j<$nmod;$j++) {
	if ($iii==$models[$j]) {
	    my $pdl_slice=$pdl_flux_c_conv->slice(",($iii)");
	    my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_slice);
	    $pdl_new=$pdl_new+$weight[$j]*$out_spec_pdl;
	    $age_now=$age_now+$weight[$j]*$age_mod[$iii];
	    $met_now=$met_now+$weight[$j]*$met_mod[$iii];
	    $age_mass_now=$age_mass_now+$w_mass[$j]*$age_mod[$iii];
	    $met_mass_now=$met_mass_now+$w_mass[$j]*$met_mod[$iii];
#	    print "$i $nmod $age_now $age_mod[$iii]\n";
	}
    }
}

$outfile_clean=$outfile.".clean";
open(OUT,">$outfile");
open(OUT_CLEAN,">$outfile_clean");
for ($j=0;$j<$npix0-1;$j++) {
    $k=$j+1;
    $val=$pdl_new->at($j);
    $val=$val*$dust[$j]*$flux_mean;
    $noise=$pdl_rms_now->at($j);
    $val_clean=$val;
    if ($noise>10*$val) {
	$noise=1;
    }
    $val=$val+$noise;
#    print "$k $wave_unc[$j] $val $dust[$j]\n";
    print OUT "$k $wave_unc[$j] $val\n";
    print OUT_CLEAN "$k $wave_unc[$k] $val_clean\n";

}
close(OUT_clean);
close(OUT);

print "$age_now $met_now $age_mass_now $met_mass_now\n";

exit;

