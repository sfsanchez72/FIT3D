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
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<6) {
    print "USE: median_clean.pl INPUT.fits MEDIAN_CLEAN.fits med_X_box med_Y_box reject_frac x_box_limit y_box_limit [niter] [llimit] [hlimit]\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$xbox=$ARGV[2];
$ybox=$ARGV[3];
$r_frac=$ARGV[4];
$xbox2=$ARGV[5];
$ybox2=$ARGV[6];
$llimit=-1e12;
$hlimit=1e12;
$niter=1;
if ($#ARGV==7) {
    $niter=$ARGV[7];
}
if ($#ARGV==8) {
    $niter=$ARGV[7];
    $llimit=$ARGV[8];
}
if ($#ARGV==9) {
    $niter=$ARGV[7];
    $llimit=$ARGV[8];
    $hlimit=$ARGV[9];
}


#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
#$nax=read_naxes($infile);   
#@naxis=@$nax;
#$nx=$naxis[0];
#$ny=$naxis[1];s
#@in_array=read_img($infile);

$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
#@out_array=list($a_in);
#$a_in=0;
#print "$nx $ny\n";
for ($n=0;$n<$niter;$n++) {
    $a_smooth=med2df($a_in,$xbox,$ybox,{Boundary => Reflect});
    $a_med=med2df($a_in,$xbox2,$ybox2,{Boundary => Reflect});
    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	    $val=$a_in->at($i,$j);
	    $median=$a_smooth->at($i,$j);
	    $max=$median+$r_frac*$median;
	    $min=$median-$r_frac*$median;
	    $smooth=$a_smooth->at($i,$j);
	    $out=$a_med->at($i,$j);
	    if (($val<$max)&&($val>$min)) {
		$out_array[$j][$i]=$val;
	    } else {	    
		$out_array[$j][$i]=$smooth;
	    }
	    if ($val<$llimit) {
#	    print "$val $llimit\n";
		$out_array[$j][$i]=$out;
	    }
	    if ($val>$hlimit) {
		$out_array[$j][$i]=$out;
	    }
	}
	if ($i==(50*int($i/50))) {
	    print "$i/$nx\n";
	}
    }
    $a_in=pdl(@out_array);
}

print "Writting the median cleaned file $out_file\n";
system("rm $out_file");
#write_img($out_file,$nx,$new_ny,\@out_array);
$pdl_out=pdl(@out_array);
$pdl_out->sethdr($h);
$pdl_out->wfits($out_file);



print "DONE\n";

exit;
