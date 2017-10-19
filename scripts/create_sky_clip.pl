#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
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




$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<4) {
    print "USE: create_sky.pl INPUT SKY WIDTH_BOX NSIGMA PLOT\n";
    exit;
}

$inputfile=$ARGV[0];
$skyfile=$ARGV[1];
$npoly=$ARGV[2];
$nsigma=$ARGV[3];
$spy=$ARGV[4];

print "Reading file\n";
$pdl=rfits($inputfile);
($n1,$n2)=$pdl->dims;
$pdl_sky=zeroes($n1,$n2);
print "Done\n";
print "$n1 $n2\n";



for ($j=0;$j<$n2;$j++) {
    $new_pos[$j]=$j;
    $pos[$j]=$j;
}

my @sky_331;
for ($i=0;$i<$n1;$i++) {
    my @array;
    $y_min=1e12;
    $y_max=-1e12;
    $pdl_sec=$pdl->slice("$i,:");
#    print "$pdl_sec\n";
    ($y_min,$y_max)=$pdl_sec->minmax;
    @nn=$pdl_sec->dims;
    my @array=list($pdl_sec);
    my @marray=median_clip1d(\@array,$n2,$npoly,$nsigma);
#($pdl_array,$coeff) = fitpoly1d(pdl(@pos),pdl(@marray),$npoly);
#($pdl_array,$coeff) = fitpoly1d(pdl(@pos),pdl(@array),$npoly);
#my $spline=new Math::Spline(\@x_new,\@x_old);
    #my $pdl_array=pdl(@array);
    
#my @new_array=list($pdl_array);
    if ($spy==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.6);           # Set character height
	pgenv(0,$n2,$y_min,$y_max,0,0);
	pgsch(1.6);           # Set character height
	pglabel("Spec.Id","Counts","Row $i");
	pgsch(2.2);           # Set character height
	pgsci(1);
	pgpoint($n2,\@pos,\@array,1);
	pgsci(3);
	pgpoint($n2,\@pos,\@marray,3);
    #pgsci(4);
	#pgpoint($n2,\@pos,\@array,2);
#
	pgsci(2);
#    pgline($n2,\@new_pos,\@new_array);    
	pgsci(1);
	pgclose;
	pgend;
#    <stdin>;
    }
    $ji=0;
    
    my @mmarray=median_filter(5,\@marray);
    for ($j=0;$j<$n2;$j++) {
	$val=$mmarray[$j];
	set($pdl_sky,$i,$j,$val);
    }
#    my $t = $pdl_sky->slice("$i,:");
#    $t .= $pdl_array; 
    $pdl_sky->wfits($skyfile);
    print "$i/$n1\n";
}

#system("rm $skyfile");



exit;

