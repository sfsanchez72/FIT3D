#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#



use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<3) {
    print "USE: rss2cube.pl input_RSS.fits pos_table.txt output_cube.fits final_dpix\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$output_cube=$ARGV[2];
$dpix=$ARGV[3];

$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$dpix_ini=$data[1]/$dpix;
$shape=$data[0];
if ($shape eq "R") {
    $area=$dpix_ini*$dpix_ini;
}
if ($shape eq "C") {
    $area=3.1416*($dpix_ini*$dpix_ini);
}
if ($shape eq "H") {
    $area=3.1416*($dpix_ini*$dpix_ini);
}

$nx_min=1e12;
$ny_min=1e12;
$nx_max=-1e12;
$ny_max=-1e12;

while($line=<PT>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $xx[$n]=$data[1]/$dpix;
    $yy[$n]=$data[2]/$dpix;
    if ($xx[$n]>$nx_max) {
	$nx_max=$xx[$n];
    }
    if ($yy[$n]>$ny_max) {
	$ny_max=$yy[$n];
    }
    if ($xx[$n]<$nx_min) {
	$nx_min=$xx[$n];
    }
    if ($yy[$n]<$ny_min) {
	$ny_min=$yy[$n];
    }
#    print "$xx[$n] $yy[$n] $nx_min,$nx_max $ny_min,$ny_max\n";
    $n++;
}
close(PT);


$nx=intval($nx_max-$nx_min+1+$dpix_ini);
$ny=intval($ny_max-$ny_min+1+$dpix_ini);
for ($j=0;$j<$n;$j++) {
    $xx[$j]=$xx[$j]-$nx_min;
    $yy[$j]=$yy[$j]-$ny_min;
    $ii=intval($xx[$j]);
    $jj=intval($yy[$j]);
#    print "$id[$j] $xx[$j] $yy[$j] $ii,$jj\n";
}

#exit;
#$nx=$ARGV[2];
#$ny=$ARGV[3];


$nax=read_naxes($input_rss);   
@naxis=@$nax;
@a_rss=read_img($input_rss);
$check=$nx*$ny;
($start_w,$delta_w)=read_img_headers($input_rss,["CRVAL1","CDELT1"]);
if ($delta_w==0) {
    $delta_w=1;
}

$pdl=zeroes($nx,$ny,$naxis[0]);
#$pdl=pdl(@a_cube);
#    for ($j=0;$j<$ny;$j++) {
#	for ($i=0;$i<$nx;$i++) {
#	    $a_cube[$k][$j][$i]=$a_rss[$j+($nx-1-$i)*$ny][$k];	    
for ($nn=0;$nn<$n;$nn++) {
    $xc=$xx[$nn];
    $yc=$yy[$nn];
#    print "$nn *********************\n";
    for ($i=($xc-1.5*$dpix_ini);$i<($xc+1.5*$dpix_ini);$i++) {
	for ($j=($yc-1.5*$dpix_ini);$j<($yc+1.5*$dpix_ini);$j++) {
	    $d=sqrt(($i-$xc)**2+($j-$yc)**2);
	    $ii=intval($i);
	    $jj=intval($j);

	    if (($d<$dpix_ini)&&($ii>=0)&&($ii<$nx)&&($jj>=0)&&($jj<$ny)) {
#	    if (($d<ceil($dpix_ini))&&($i>=0)&&($i<$nx)&&($j>=0)&&($j<$ny)) {
#		print "$nn ($xx[$nn],$yy[$nn]) ($i,$j), [$nx,$ny] \n";
		for ($k=0;$k<$naxis[0];$k++) {   
		    set($pdl,$ii,$jj,$k,$a_rss[$nn][$k]/$area);
#		    $a_cube[$k][$j][$i]=$a_rss[$nn][$k];	    	    
#		    $a_done[$j][$i]=1;
		}
	    }
	}
    }
    print "$nn/$n\n";
}
#print "PASO\n";
$nax=[$nx,$ny,$naxis[0]];
#$pdl=pdl(@a_cube);
$pdl->hdr->{CRPIX1}=1;
$pdl->hdr->{CRVAL1}=$nx_min;
$pdl->hdr->{CDELT1}=$dpix;
$pdl->hdr->{CRPIX2}=1;
$pdl->hdr->{CRVAL2}=$ny_min;
$pdl->hdr->{CDELT2}=$dpix;
$pdl->hdr->{CRPIX3}=1;
$pdl->hdr->{CRVAL3}=$start_w;
$pdl->hdr->{CDELT3}=$delta_w;


$pdl->wfits($output_cube);
#write_fits($output_cube,$nax,3,\@a_cube);
#write_crval_cube($output_cube,[1,0,1,1,0,1,1,$start_w,$delta_w]);
exit;


sub intval {
    my $a=@_[0];
    my $ia=int($a);
    my $d=$a-$ia;
    if ($d>0.5) {
	$ia=$ia+1;
    }
    return $ia;
}
