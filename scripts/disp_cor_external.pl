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
use PDL::Func;
use PDL::Math;

#



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<5) {
    print "USE: disp_cor_external.pl EXTRACTED.fits CRVAL CDELT OUTPUT.fits DISPERSION.txt plot\n";
    exit;
}

$infile=$ARGV[0];
$crval=$ARGV[1];
$cdelt=$ARGV[2];
$out_file=$ARGV[3];
$disp_file=$ARGV[4];
$plot=$ARGV[5];

print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";

#
# We load the dispersion correction
#
print "Loading the dispersion solution from $disp_file\n";
open(FH,"<$disp_file");
$n_cut=0;
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$a_i_old[$n_cut]=$data[2];
	$n_cut++;
    }
}
close(FH);
print "DONE\n";

#exit;
#
# We correct for the new dispersion
#
my @out_array;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$spec[$i]=$in_array[$j][$i];
	$wave[$i]=$i;
    }

    my $pdl_wave=pdl(@wave);
    my $pdl_spec=pdl(@spec);
    my $pdl_a_i_old=pdl(@a_i_old);
#    my $INTER = init PDL::Func( Interpolate => "Hermite");
#    $INTER->set( x => $pdl_wave, y => $pdl_spec, bc => "simple" );              
#    my $out_spec_pdl = $INTER->interpolate($pdl_a_i_old);             


    my $out_spec_pdl = interpol(pdl(@a_i_old), pdl(@wave), pdl(@spec));
    $min=1e12;
    $max=-1e12;
    for ($i=0;$i<$n_cut;$i++) {
	$wave_new[$i]=$crval+$cdelt*$i;
	$out_spec[$i]=$out_spec_pdl->at($i);
	$out_array[$j][$i]=$out_spec[$i];
	if ($min>$spec[$i]) {
	    $min=$spec[$i];
	}
	if ($max<$spec[$i]) {
	    $max=$spec[$i];
	}
    }

    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($crval,$crval+$cdelt*$n_cut,$min,$max,0,0);
	pglabel("Wavelength","Counts","$j/$ny");
	pgsci(1);
	pgline($n_cut,\@wave_new,\@out_spec);
	pgclos();
	pgend();
	print "Press Enter";
	$command=<stdin>;
	chop($command);
	if ($command eq "A") {
	    $plot=0;
	}
    }    
}


print "Writting the $out_file\n";
system("rm $out_file");
write_rss($out_file,$n_cut,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;
