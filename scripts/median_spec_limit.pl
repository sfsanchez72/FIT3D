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
    print "USE: median_spec.pl RSS.fits MEDIAN_SPEC.fits LIMIT [NY]\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$limit=$ARGV[2];
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

#
# We define the minimul range on the basis of a possible vignetting
#

$nx1=$nx;
$nx2=0;
for ($i=0;$i<$nx;$i++) {
    $val=$in_array[0][$i];
    if (($i>0)&&($i<$nx)) {
	if (($val>$limit)&&($old_val<=$limit)&&($nx2==0)) {
	    $nx1=$i;
	}
	if (($old_val>$limit)&&($val<=$limit)&&($nx1<$nx)) {
	    $nx2=$i;
	}
    }
#    print "$limit $nx1 $nx2 $val\n";
    $old_val=$val;
}
$nx1=$nx1+5;
$nx2=$nx2-5;
print "Range ($nx1,$nx2) (0,$nx)\n";

my @mask;
my @mean_vals;
$Jmin=0;
$min_val=1e12;
for ($j=0;$j<$ny-1;$j++) {
    $mask[$j]=1;
    my @tmp;
    my $k=0;
    for ($i=$nx1;$i<$nx2;$i++) {
	$tmp[$k]=$tmp[$k]+$in_array[$j][$i];
	$k++;
    }
    $mean_vals[$j]=median(@tmp);
    if (($min_val>$mean_vals[$j])&&($mean_vals[$j]>0)) {
	$Jmin=$j;
	$min_val=$mean_vals[$j];
    }
}
$med_mean_val=mean(@mean_vals);
$sig_mean_val=sigma(@mean_vals);

#($MIN,$MAX)=minmax(@mean_vals);
#print "$med_mean_val $sig_mean_val $MIN $MAX\n";
#if ((($MAX-$MIN)/$MIN)>1) {
for ($j=0;$j<$ny;$j++) {
    if (($mean_vals[$j]-$mean_vals[$Jmin])>1.0*$mean_vals[$Jmin]) {
	$mask[$j]=0;
    } 

    if (abs($mean_vals[$j]-$med_mean_val)>2*$sig_mean_val) {
	$mask[$j]=0;
    } 
    
#    if (abs($mean_vals[$j]-$med_mean_val)>0.7*$sig_mean_val) {
#	$mask[$j]=0;
#    } 

    if ($mean_vals[$j]==0) {
	$mask[$j]=0;
    }
#    print "$j $Jmin $mean_vals[$j]-$mean_vals[$Jmin])>1.5*$mean_vals[$Jmin] $mask[$j]\n";
}

#} else {
#    for ($j=0;$j<$ny;$j++) {

#	if ($mean_vals[$j]==0) {
#	    $mask[$j]=0;
#	}
#
#    }
#}




my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    $kk=0;
    for ($j=0;$j<$ny;$j++) {
	$val=$in_array[$j][$i];
	if (($val>$limit)&&($mask[$j]==1)) {
	    $cut[$kk]=$val;
	    $kk++;
	}
    }
    $spec[$i]=median(@cut);
    for ($j=0;$j<$new_ny;$j++) {
	$out_array[$j][$i]=$spec[$i];
    }
}

print "Writting the median spectrum $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$new_ny,\@out_array);
print "DONE\n";

exit;
