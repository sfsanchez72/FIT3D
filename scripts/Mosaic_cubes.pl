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


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<1) {
    print "USE: Mosaic.pl CONFIG.txt OUTPUT [CUT] [NBOX]\n";
    exit;
}

$config=$ARGV[0];
$outfile=$ARGV[1];
$cut=0;
$box=0;
if ($#ARGV==2) {
    $cut=$ARGV[2];
}
if ($#ARGV==3) {
    $cut=$ARGV[2];
    $box=$ARGV[3];
}
$nf=0;
$dxmin=1e12;
$dxmax=-1e12;
$dymin=1e12;
$dymax=-1e12;
open(FH,"<$config");
while($line=<FH>) {
#    chop($line);
    @data=split(" ",$line);
    $infile[$nf]=$data[0];
    $dx[$nf]=$data[1];
    $dy[$nf]=$data[2];
    if ($dxmin>$dx[$nf]) {
	$dxmin=$dx[$nf];
    }
    if ($dxmax<$dx[$nf]) {
	$dxmax=$dx[$nf];
    }
    if ($dymin>$dy[$nf]) {
	$dymin=$dy[$nf];
    }
    if ($dymax<$dy[$nf]) {
	$dymax=$dy[$nf];
    }
    $nf++;
}
close(FH);
print "$nf files\n";
#$nax=read_naxes($infile[0]);   
#@naxis=@$nax;
$pdl_in=rfits($infile[0]);
@naxis=$pdl_in->dims();
$nx=$naxis[0];
$ny=$naxis[1];
$nz=$naxis[2];
print "[$nx,$ny,$nz]\n";
$h=$pdl_in->gethdr;



#@tr=read_img_headers($infile[0],["CRPIX1","CRVAL1","CDELT1","CRPIX2","CRVAL2","CDELT2","CRPIX3","CRVAL3","CDELT3"]);

$xmin=int($dxmin);
$xmax=int($nx+$dxmax);
$ymin=int($dymin);
$ymax=int($ny+$dymax);

$nx_new=$xmax-$xmin;
$ny_new=$ymax-$ymin;

$h->{NAXIS1}=$nx_new;
$h->{NAXIS2}=$ny_new;


print "$xmin $xmax $ymin $ymax $nx_new $ny_new\n";

$pdl_out=zeroes($nx_new,$ny_new,$nz);
$pdl_sum=zeroes($nx_new,$ny_new);
for ($f=0;$f<$nf;$f++) {
    print "Reading file $infile[$f], ";    
    my $pdl_in=rfits($infile[$f]);
    $nx_start=$dx[$f]-$dxmin;
    $ny_start=$dy[$f]-$dymin;
    print "added at $nx_start,$ny_start\n";
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
	    $nz_med=int($nz/2);
#	    print "$i,$j,$nz_med\n";
	    my $val=$pdl_in->at($i,$j,$nz_med);
	    $clean=0;
	    if ($box>0) {
		$dist=sqrt($i**2+$j**2);
		if ($dist<$box) {
		    $clean=1;
		}
		$dist=sqrt(($nx-$i)**2+$j**2);
		if ($dist<$box) {
		    $clean=1;
		}
		$dist=sqrt($i**2+($ny-$j)**2);
		if ($dist<$box) {
		    $clean=1;
		}
		$dist=sqrt(($nx-$i)**2+($ny-$j)**2);
		if ($dist<$box) {
		    $clean=1;
		}
		$j_min=$j-$box;
		$j_max=$j+$box;
		$i_min=$i-$box;
		$i_max=$i+$box;
		if ($j_min<0) {
		    $j_min=0;
		}
		if ($j_max>=$ny) {
		    $j_max=$ny-1;
		}
		if ($i_min<0) {
		    $i_min=0;
		}
		if ($i_max>=$nx) {
		    $i_max=$nx-1;
		}
		for ($ib=$i_min;$ib<$i_max;$ib++) {
		    for ($jb=$j_min;$jb<$j_max;$jb++) {
			my $val2=$pdl_in->at($ib,$jb,$nz_med);
			if ($val2==0) {
			    $clean=1;
			}
		    }
		}
	    }
	    if (($val>$cut)&&($clean==0)) {
		$ii=$i+$nx_start;
		$jj=$j+$ny_start;

		my $sec=$pdl_out->slice("($ii),($jj),:");
		my $pdl=$pdl_in->slice("($i),($j),:");
		$sec .= $sec+$pdl;
		my $sum=$pdl_sum->at($ii,$jj);
		$sum=$sum+1;
		set($pdl_sum,$ii,$jj,$sum);
	    }
	}
    }    
}
#print "Done 1\n";

for ($j=0;$j<$ny_new;$j++) {
    for ($i=0;$i<$nx_new;$i++) {
	my $sum=$pdl_sum->at($i,$j);
	if ($sum>0) {
	    my $sec=$pdl_out->slice("($i),($j),:");
	    my $pdl=$pdl_in->slice("($i),($j),:");
	    $sec .= $sec/$sum;
	}
    }
}


system("rm $outfile");
$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);
exit;


