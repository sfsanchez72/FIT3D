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
use PDL::Func;
use PDL::Math;

#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<4) {
    print "USE: mdist_cor_external.pl EXTRACTED.fits DISTORTION.txt DISTORSION.fits OUTPUT.fits plot [offset]\n";
    exit;
}

$infile=$ARGV[0];
$dist_txt=$ARGV[1];
$dist_file=$ARGV[2];
$out_file=$ARGV[3];
$plot=$ARGV[4];
$offset=0;
if ($#ARGV==5) {
    $offset=$ARGV[5];
}



print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";

if ($dist_txt ne "none") {
    print "Reading the distortion correction\n";
    open(FH,"<$dist_txt");
    $nd=0;
    $shift_min=1e12;
    $shift_max=-1e12;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$index_all[$nd]=$data[0];
	$shift_all[$nd]=$data[1]+$offset;
	if ($shift_min>$shift_all[$nd]) {
	    $shift_min=$shift_all[$nd];
	}
	if ($shift_max<$shift_all[$nd]) {
	    $shift_max=$shift_all[$nd];
	}
	$nd++;
    }
    close(FH);
    print "DONE\n";
} else {
    for ($nd=0;$nd<$ny;$nd++) {
	$shift_all[$nd]=0;
    }
}

#$ny_cen=int($ny/2);
#$shift_cen=$shift_all[$ny_cen];
#for ($j=0;$j<$nd;$j++) {
#    $shift_all[$j]=$shift_all[$j]-$shift_cen;
#    print "$j $shift_all[$j]\n";
#}


print "Reading 2nd order distortion $dist_file ";
@a_dist=read_img($dist_file);
print "DONE\n";
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$a_dist[$j][$i]=$a_dist[$j][$i]+$shift_all[$j];
    }
}

#
# We correct for the new dispersion
#
my @out_array;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$a_i_old[$i]=$a_dist[$j][$i];
	$spec[$i]=$in_array[$j][$i];
	$wave[$i]=$i;
    }
    my $pdl_wave=pdl(@wave);
    my $pdl_spec=pdl(@spec);
    my $pdl_a_i_old=pdl(@a_i_old);
    my $INTER = init PDL::Func( Interpolate => "Hermite");
    $INTER->set( x => $pdl_wave, y => $pdl_spec, bc => "simple" );
    my $out_spec_pdl = $INTER->interpolate($pdl_a_i_old);

#    my $out_spec_pdl = interpol(pdl(@a_i_old), pdl(@wave), pdl(@spec));
    $min=1e12;
    $max=-1e12;
    for ($i=0;$i<$nx;$i++) {
	$wave_new[$i]=$i;
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
	pgenv(0,$nx,$min,$max,0,0);
	pglabel("Wavelength","Counts","$j/$ny");
	pgsci(2);
	pgline($nx,\@wave_new,\@spec);
	pgsci(1);
	pgline($nx,\@wave_new,\@out_spec);
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
write_img($out_file,$nx,$ny,\@out_array);
#write_rss($out_file,$n_cut,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;

