#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

#use Statistics::OLS;
#use Math::FFT;
#use Math::Stat;
#use Math::Spline qw(spline linsearch binsearch);
#use Math::Derivative qw(Derivative2);
#use Math::Approx;
#use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";


if ($#ARGV<5) {
    print "USE: galfit3d.pl INPUT_CUBE.fits MOD_CUBE.fits RES_CUBE.fits GALFIT_CONF LOGFILE QUALITY.txt\n";
    print "GALFIT_CONF: A) object.fits, B) output.fits\n";
    exit;
}

$input_cube=$ARGV[0];
$mod_cube=$ARGV[1];
$res_cube=$ARGV[2];
$galfit_conf=$ARGV[3];
$logfile=$ARGV[4];
$qfile=$ARGV[5];

#$pdl=rfits($input_cube,{data=>0});

system("rm output.fits");
system("rm fit.log");
#
#We create the output files
#

$pdl_input_cube=rfits($input_cube);
$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$ny=$pdl_input_cube->hdr->{"NAXIS2"};
$nz=$pdl_input_cube->hdr->{"NAXIS3"};
$pdl_mod_cube=rfits($input_cube);
$pdl_res_cube=rfits($input_cube);

#system("rm $mod_cube");
#system("rm $res_cube");
#open(TMP,">tmp.cl");
#print TMP "\n";
#print TMP "\n";
#print TMP "\n";
#print TMP "imcopy $input_cube $mod_cube\n";
#print TMP "imcopy $input_cube $res_cube\n";
#print TMP "\n";
#print TMP "logout\n";
#print TMP "\n";
#close(TMP);
#$call=$cl." < tmp.cl > junk.tmp";
#system($call);
    



#
# Loop over the Z direction:
#
for ($k=0;$k<$nz;$k++) {
    $kk=$k+1;
    print "$k/$nz...";
    my $slice=$pdl_input_cube->slice(",,($k)");
#    $slice->wfits("object.fits");
    system("rm object.fits");  
    open(TMP,">tmp.cl");
    print TMP "\n";
    print TMP "imcopy $input_cube\[*,*,$kk\] object.fits\n";
    print TMP "\n";
    print TMP "logout\n";
    print TMP "\n";
    close(TMP);
    $call=$cl." < tmp.cl > junk.tmp";
    system($call);
   print "copied...";
    $call=$galfit." ".$galfit_conf." > galfit.junk";
#    $call=$galfit." ".$galfit_conf." ";
    system($call);
    system ("rm galfit.??");
    system ("rm galfit.junk");
#    $call=$galfit." ".$galfit_conf;
#    system($call);
    print "modelled...";

#    my $pdl_output=rfits("output.fits");
    my $mod;
    open(DIR,"ls output.fits |");
    $dir=<DIR>;
    close(DIR);
    chop($dir);
#    print "\n'$dir'\n";
    $fit=0.0;
#    $test=rfits("output.fits");
#    $btest=$test->dims;
#    print "BTEST=$btest\n";
    if ($dir !~ "output.fits" ) {
	$mod=zeroes($nx,$ny);
    } else {
	$mod=rfits("output.fits[2]");
	$fit=1.0;
    }    
    ($nx_c,$ny_c)=$mod->dims;

#    print "PASO\n";
#    $mod->wfits("test_mod.fits");
    for ($i=0;$i<$nx_c;$i++) {
	for ($j=0;$j<$ny_c;$j++) {	
	    $val=$mod->at($i,$j);
	    set($pdl_mod_cube,$i,$j,$k,$val);
	}
    }
    my $res;
    if ($dir !~ "output.fits" ) {
	$res=$slice;
    } else {
	$res=rfits("output.fits[3]");
    }    

#=rfits("output.fits[3]");
#    my $res=$pdl_output->slice(",,(2)");
#    $res->wfits("test_res.fits");
    for ($i=0;$i<$nx_c;$i++) {
	for ($j=0;$j<$ny_c;$j++) {	
	    $val=$res->at($i,$j);
	    set($pdl_res_cube,$i,$j,$k,$val);
	}
    }
#    open(TMP,">tmp.cl");
#    print TMP "\n";
#    print TMP "imcopy output.fits\[2\] $mod_cube\[*,*,$kk\]\n";
#    print TMP "\n";
#    print TMP "imcopy output.fits\[3\] $res_cube\[*,*,$kk\]\n";
#    print TMP "\n";
#    print TMP "logout\n";
#    print TMP "\n";
#    print TMP "display output.fits\[1\] 1 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[2\] 2 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[3\] 3 zs- zr- z1=-0.001 z2=0.02\n";
#    close(TMP);
#    $call=$cl." < tmp.cl > junk.tmp";
#    system($call);
    system("cp output.fits out_put.fits");
    print "DONE\n";
    open(TMP,">display.cl");
    print TMP "\n";
    print TMP "display out_put.fits\[1\] 1 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display out_put.fits\[2\] 2 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display out_put.fits\[3\] 3 zs- zr- z1=-0.001 z2=0.02\n";
    close(TMP);
    system ("rm galfit.??");  
    system("rm output.fits");
    system("cp fit.log $logfile");
    if ($k==0) {
	open(Q,">$qfile");
    } else {
	open(Q,">>$qfile");
    }
    print Q "$k $fit\n";
    close(Q);
    print "coping cubes...";
    $pdl_mod_cube->wfits($mod_cube);
    $pdl_res_cube->wfits($res_cube);
    print "done\n";
}


print "Model done\n";
$pdl_res_cube->wfits($res_cube);
print "$mod_cube and $res_cube saved\n";
exit;




