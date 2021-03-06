#!/usr/bin/perl

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

if ($#ARGV<7) {
    print "USE: kin_back_rss.pl rss.fits config_file back_file nback start_wavelength end_wavelength out_file DEVICE [MASK_FILE]\n";
    exit;
}

$input_cube=$ARGV[0];
$config_file=$ARGV[1];
$back_file=$ARGV[2];
$nback=$ARGV[3];
$start_w=$ARGV[4];
$end_w=$ARGV[5];
$out_file=$ARGV[6];
$device=$ARGV[7];
if ($#ARGV==8) {
    $mask_list=$ARGV[8];
}

if ($mask_list eq "none") {
    $nmask=0;
} else {
    open(FH,"<$mask_list");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$start_mask[$nmask]=$data[0];
	$end_mask[$nmask]=$data[1];
	$nmask++;
    }
    close(FH);
}

print "$nmask regions\n";


#
# We create a temporary BACK fits file
#
$nf=0;
open(FH,"<$back_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$file[$nf]=$line;
	$nf++;
    }
}
close(FH);
print "NF=$nf\n";

for ($j=0;$j<$nf;$j++) {
    open(FH,"<$file[$j]");
    $nx=0;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($j==0) {
	    if ($nx==0) {
		$crval=$data[1];
	    }
	    if ($nx==1) {
		$cdelt=$data[1]-$crval;
	    }
	}
	$id[$nx]=$data[0];
	$w[$nx]=$data[1];
	$f[$nx]=$data[2];
	$nx++;
    }
    close(FH);
    if ($j==0) {
	$pdl_out=zeroes($nx,$nf);
    }
    for ($i=0;$i<$nx;$i++) {
	set($pdl_out,$i,$j,$f[$i]);
    }
}
$$h{"CRPIX1"}=1;
$$h{"CRVAL1"}=$crval;
$$h{"CDELT1"}=$cdelt;

if ($nf>0) {
    $pdl_out->sethdr($h);
    $pdl_out->wfits("BACK_JUNK.fits");
}
#
# DONE
#

$pdl_input_cube=rfits($input_cube);
$h=$pdl_input_cube->gethdr;
$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$ny=$pdl_input_cube->hdr->{"NAXIS2"};
$crpix=$pdl_input_cube->hdr->{"CRPIX1"};
$crval=$pdl_input_cube->hdr->{"CRVAL1"};
$cdelt=$pdl_input_cube->hdr->{"CDELT1"};


$call="cp ".$config_file." tmp.config";
system($call);
for ($j=0;$j<$ny;$j++) {
    open(FHOUT,">fit_spectra.input");
    $nx_out=0;
    for ($i=0;$i<$nx;$i++) {
	$w=$cdelt*($i-$crpix+1)+$crval;
	$f=$pdl_input_cube->at($i,$j);
	if (($w>$start_w)&&($w<$end_w)) {
	    $masked=0;
	    for ($jm=0;$jm<$nmask;$jm++) {
		if (($w>$start_mask[$jm])&&($w<$end_mask[$jm])) {
		    $masked=1;
		    
		}
	    }
	    if ($masked==0) {
		printf FHOUT "$i $w $f\n";
		if ($nx_out==0) {
		    $crpix_out=$i+1;
		}
		$nx_out++;
	    }
	}
    }
    close(FHOUT);
    if ($j==0) {
	system("emacs tmp.config &");
    }
    
    if ($auto !~ "y") {
	$command="r";
	$n_fit=1;
	while ($command !~ "s") {
	    if ($command =~ "r") {
		print "Fitting Num. $n_fit of the spectrum ($i,$j)\n";
#		$call="fit_spec_back fit_spectra.input tmp.config ".$back_file." ".$nback." ".$device;
		if ($nf>0) {
		    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
		} else {
		    $call="fit_spec_back fit_spectra.input tmp.config none 0 ".$device;
		}
		print "CALL=$call\n";
		system($call);
		$n_fit++;
		}
	    print "Options:\n";
		print "[s] save the results\n";
	    print "[r] repeat the fitting\n";
	    print "[m] copy the output of the last fit to the new config\n";
	    print "[a] Set the automatic fitting\n\n";
	    $command=<stdin>;
	    chop($command);
	    
	    if ($command =~ "a") {
		$command="s";
		$auto="y";
	    }
	    
	    if ($command =~ "m") {
		system("cp out_config.fit_spectra tmp.config");
		system("emacs tmp.config &");
	    }	    
	}
    } else {
	print "Fitting the spectrum ($i,$j)\n";
	#$call="fit_spec_back fit_spectra.input tmp.config ".$back_file." ".$nback." ".$device;
	if ($nf>0) {
	    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
	} else {
	    $call="fit_spec_back fit_spectra.input tmp.config none 0 ".$device;
	}
	system($call);	    
    }
#
# We copy the output
#
    open(FH,"<out.fit_spectra");
    if ($j==0) {
	open(FHOUT,">$out_file");
    } else {
	open(FHOUT,">>$out_file");
    }
    while($line=<FH>) {
	print FHOUT "$line";
    }
    close(FHOUT);
    close(FH);
#
# We copy the residuals
#	
    print "$nx_out\n";
    open(FH,"<out_mod_res.fit_spectra");
    if ($j==0) {
	$org=zeroes($nx_out,$ny);
	$mod=zeroes($nx_out,$ny);
	$res=zeroes($nx_out,$ny);
    }
    $n_line=0;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($n_line==0) {
	    $start_out_w=$data[0];
	}
	if ($n_line==1) {
	    $delta_out_w=$data[0];
	}
	set($org,($n_line,$j),$data[1]);
	set($mod,($n_line,$j),$data[2]);
	set($res,($n_line,$j),$data[3]);
	$n_line++;
    }
    close(FH);
	
}
$h->{CRPIX1}=$crpix_out;
$h->{NAXIS1}=$nx_out;
#$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$org->sethdr($h);
$org->wfits("kin_back_rss_org.fits");
$res->sethdr($h);
$res->wfits("kin_back_rss_res.fits");
$mod->sethdr($h);
$mod->wfits("kin_back_rss_mod.fits");



exit;
