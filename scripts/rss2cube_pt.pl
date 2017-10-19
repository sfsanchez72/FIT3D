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

if ($#ARGV<2) {
    print "USE: rss2cube.pl input_RSS.fits pos_table.txt output_cube.fits\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$output_cube=$ARGV[2];

$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$dpix=$data[1];
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
    $n++;
}
close(PT);

#
# Error???
#
$nx=intval($nx_max-$nx_min+1);
$ny=intval($ny_max-$ny_min+1);
for ($j=0;$j<$n;$j++) {
    $xx[$j]=$xx[$j]-$nx_min;
    $yy[$j]=$yy[$j]-$ny_min;
    $ii=intval($xx[$j]);
    $jj=intval($yy[$j]);
#    print "$id[$j] $xx[$j] $yy[$j] $ii,$jj $nx $ny\n";
}

#exit;
#$nx=$ARGV[2];
#$ny=$ARGV[3];

$pdl_rss=rfits($input_rss);
@naxis=$pdl_rss->dims();
#$nax=read_naxes($input_rss);   
#@naxis=@$nax;
#@a_rss=read_img($input_rss);
$check=$nx*$ny;

$start_w=$pdl_rss->hdr->{CRVAL1};
$delta_w=$pdl_rss->hdr->{CDELT1};

#($start_w,$delta_w)=read_img_headers($input_rss,["CRVAL1","CDELT1"]);
if ($delta_w==0) {
    $delta_w=1;
}
if ($check!=$naxis[1]) {
    print "Error: The Y-axis of the RSS has $naxis[1] values\n";
    print "It does not correspond with $nx X $ny\n";
    exit;
}

$pdl=zeroes($nx,$ny,$naxis[0]);

for ($nn=0;$nn<$n;$nn++) {
    $i=intval($xx[$nn]);
    $j=intval($yy[$nn]);
    my $t = $pdl->slice("($i),($j),:");
    $t .= $pdl_rss->slice(":,($nn)");
    print "$i $j $xx[$nn] $yy[$nn] $nx_min,$nx_max $ny_min,$ny_max\n";

#    for ($k=0;$k<$naxis[0];$k++) {    
#	$a_cube[$k][$j][$i]=$a_rss[$nn][$k];	    	    
#	$a_done[$j][$i]=1;
#    }
}
#$nax=[$nx,$ny,$naxis[0]];
$nax=[$nx,$ny,$naxis[0]];
#system("rm $output_cube");
#print "PASO\n";
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
