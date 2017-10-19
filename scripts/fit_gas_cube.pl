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
    print "USE: kin_back_cube.pl cube.fits config_file back_file nback start_wavelength end_wavelength out_file DEVICE [MASK_FILE] [ASK]\n";
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
$ask=0;
if ($#ARGV==9) {
    $mask_list=$ARGV[8];
    $ask=$ARGV[9];
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
    $pdl_out->sethdr( $h );
    $pdl_out->wfits("BACK_JUNK.fits");
}

$pdl_input_cube=rfits($input_cube);
$h=$pdl_input_cube->gethdr;
$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$ny=$pdl_input_cube->hdr->{"NAXIS2"};
$nz=$pdl_input_cube->hdr->{"NAXIS3"};
$crpix=$pdl_input_cube->hdr->{"CRPIX3"};
$crval=$pdl_input_cube->hdr->{"CRVAL3"};
$cdelt=$pdl_input_cube->hdr->{"CDELT3"};


$call="cp ".$config_file." tmp.config";
system($call);
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	open(FHOUT,">fit_spectra.input");
	$nz_out=0;
	for ($k=0;$k<$nz;$k++) {
	    $w=$cdelt*($k-$crpix+1)+$crval;
	    $f=$pdl_input_cube->at($i,$j,$k);
	    if (($w>$start_w)&&($w<$end_w)) {
		$masked=0;
		for ($jm=0;$jm<$nmask;$jm++) {
		    if (($w>$start_mask[$jm])&&($w<$end_mask[$jm])) {
			$masked=1;			
		    }
		}
		if ($masked==0) {
		    if ($f eq "nan") {
			$f=0;
		    }
		    if ($f eq "BAD") {
			$f=0;
		    }
		    if ($f != 1*$f) {
			$f=0;
		    }
		    printf FHOUT "$k $w $f\n";
		    if ($nz_out==0) {
			$crpix_out=$k+1;
		    }
		    $nz_out++;
		}
	    }
	}
	close(FHOUT);
	if (($i==0)&&($j==0)&&($ask==0)) {
	    system("emacs tmp.config &");
	   
	} else {
	    if ($ask==1) {
		$auto="y";
	    }
	}

	if ($auto !~ "y") {
	    $command="r";
	    $n_fit=1;
	    while ($command !~ "s") {
		if ($command =~ "r") {
		    print "Fitting Num. $n_fit of the spectrum ($i,$j)\n";
#		    $call="fit_spec_back fit_spectra.input tmp.config ".$back_file." ".$nback." ".$device;
		    if ($nf>0) {
		    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
		} else {
		    $call="fit_spec_back fit_spectra.input tmp.config none 0 ".$device;
		}
#		    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
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
#	    $call="fit_spec_back fit_spectra.input tmp.config ".$back_file." ".$nback." ".$device;
#	    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
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
	if (($i==0)&($j==0)) {
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
	open(FH,"out_mod_res.fit_spectra");
	if (($i==0)&&($j==0)) {
	    $org=zeroes($nx,$ny,$nz_out);
	    $mod=zeroes($nx,$ny,$nz_out);
	    $res=zeroes($nx,$ny,$nz_out);
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
	    set($org,($i,$j,$n_line),$data[1]);
	    set($mod,($i,$j,$n_line),$data[2]);
	    set($res,($i,$j,$n_line),$data[3]);
	    $n_line++;
	}
	close(FH);
	
    }
}
$h->{CRPIX3}=$crpix_out;
$h->{NAXIS3}=$nz_out;
$org->sethdr($h);
$org->wfits("kin_back_cube_org.fits");
$res->sethdr($h);
$res->wfits("kin_back_cube_res.fits");
$mod->sethdr($h);
$mod->wfits("kin_back_cube_mod.fits");



exit;
