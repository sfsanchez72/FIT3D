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


if ($#ARGV<2) {
    print "USE: create_sky.pl RSS.fits SKY.fits NSIGMA [NY] [MASK_LIST.TXT]\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$nsigma=$ARGV[2]; 
#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";
$new_ny=$ny;
if ($#ARGV==3) {
    $new_ny=$ARGV[3];
}

if ($#ARGV==4) {
    $new_ny=$ARGV[3];
    $mask_list=$ARGV[4];
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


my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    for ($j=0;$j<$ny;$j++) {
	$cut[$j]=$in_array[$j][$i];
    }
    $median=median(@cut);
    $mean=mean(@cut);
    $sigma=sigma(@cut);
    my @new_spec;
    $n_new=0;
    for ($j=0;$j<$ny;$j++) {
	if ($cut[$j]<($median+$nsigma*$sigma)) {
	    $new_spec[$n_new]=$cut[$j];
	    $n_new++;
	}
    }
#    print "$i $n_new\n";

    $masked=0;
    for ($k=0;$k<$nmask;$k++) {
	if (($i>$start_mask[$k])&&($i<$end_mask[$k])) {
	    $masked=1;	    
	}
    }

    if ($masked==0) {
	if ($n_new==0) {
	    $spec[$i]=$median;
	} else {
	    $median2=median(@new_spec);
	    $mean2=mean(@new_spec);
#	$spec[$i]=(3*$median2-2*$mean2)*0.5+0.5*$median2;	
	    $spec[$i]=(2*$median2-$mean2);
	    #$spec[$i]=$median2;
	}
    } else {
	$spec[$i]=$spec[$i-1];
    }
	
    for ($j=0;$j<$new_ny;$j++) {
	$out_array[$j][$i]=$spec[$i];
    }
}

print "Writting the median spectrum $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$new_ny,\@out_array);
print "DONE\n";

exit;
