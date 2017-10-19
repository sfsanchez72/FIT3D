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


if ($#ARGV<4) {
    print "USE: spots_center.pl INPUT.FITS SPOTS_LIST(X Y) NX_BOX NY_BOX OUTPUT_SPOTS [DX DX]\n";
    exit;
}

$infile=$ARGV[0];
$spots_file=$ARGV[1];
$bx=$ARGV[2];
$by=$ARGV[3];
$spots_out=$ARGV[4];
$reg_out=$spots_out.".reg";
$DX=$ARGV[5];
$DY=$ARGV[6];

open(HDR,"read_img_header.pl $infile HIERARCH |");
while($hdr=<HDR>) {
    chop($hdr);
    @data=split(" ",$hdr);
    if ($hdr =~ "AZ_END") {
	$az=$data[5];
    }
    if ($hdr =~ "EL_END") {
	$el=$data[5];
    }
}
close(HDR);


#print "AZ,EL = $az,$el\n";

$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;

open(FH,"<$spots_file");
open(OUT,">$spots_out");
print OUT "# $az $el\n";
open(REG,">$reg_out");
print REG "# Region file format: DS9 version 4.0\n";
print REG "# Filename: $infile\n";
print REG "global color=red font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n";
print REG "physical\n";

while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$x=$data[0]+$DX;
	$y=$data[1]+$DY;
	$ii=int($x);
	$jj=int($y);
	$i_min=$ii-$bx;
	$i_max=$ii+$bx;
	$j_min=$jj-$by;
	$j_max=$jj+$by;
	if ($i_min<0) {
	    $i_min=0;
	}
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($i_max>=$nx) {
	    $i_max=$nx-1;
	}
	if ($j_max>=$ny) {
	    $j_max=$ny-1;
	}

	$D_lim=1e12;
	$val_lim=-1e12;
	for ($i=$i_min;$i<$i_max;$i++) {
	    for ($j=$j_min;$j<$j_max;$j++) {
		$val=$a_in->at($i,$j);
		$D=sqrt(($i-$x)**2+($j-$y)**2);
#		if (($val>$val_lim)&&($D<$D_lim)) {
		if ($val>$val_lim) {
		    $val_lim=$val;
		    $D_lim=$D;
		    $x_cen_peak=$i;
		    $y_cen_peak=$j;
		}

	    }
	}



	$ii=int($x_cen_peak);
	$jj=int($y_cen_peak);
	$i_min=$ii-2;
	$i_max=$ii+2;
	$j_min=$jj-2;
	$j_max=$jj+2;
	if ($i_min<0) {
	    $i_min=0;
	}
	if ($j_min<0) {
	    $j_min=0;
	}
	if ($i_max>=$nx) {
	    $i_max=$nx-1;
	}
	if ($j_max>=$ny) {
	    $j_max=$ny-1;
	}



	$x_cen=0;
	$y_cen=0;
	$I=0;
	for ($i=$i_min;$i<$i_max;$i++) {
	    for ($j=$j_min;$j<$j_max;$j++) {
		$val=$a_in->at($i,$j);
		$I=$I+$val;
		$x_cen=$x_cen+$i*$val;
		$y_cen=$y_cen+$j*$val;

	    }
	}
	if ($I>0) {
	    $x_cen=$x_cen/$I+1;
	    $y_cen=$y_cen/$I+1;
	}


	$x_cen_peak=$x_cen_peak+1;
	$y_cen_peak=$y_cen_peak+1;

#	print "$x $y $x_cen $y_cen $x_cen_peak $y_cen_peak\n";
	print OUT "$x $y $x_cen $y_cen\n";
	print REG "circle($x_cen,$y_cen,2)\n";
    }

}
close(REG);
close(OUT);
close(FH);


exit;
