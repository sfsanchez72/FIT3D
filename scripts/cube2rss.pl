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

if ($#ARGV<1) {
    print "USE: cube2rss.pl input_cube.fits output_RSS.fits [POS_TABLE.txt]\n";
    exit;
}

$input_cube=$ARGV[0];
$output_rss=$ARGV[1];
$pos_table="cube2rss.pt.txt";
if ($#ARGV==2) {
    $pos_table=$ARGV[2];
}


$in_pdl=rfits($input_cube);
$h=$in_pdl->gethdr;

open(OUT,">$pos_table");
$pix=$h->{CDELT1};
$crpix1=$h->{CRPIX1};
$crval1=$h->{CRVAL1};
$cdelt1=$h->{CDELT1};
$crpix2=$h->{CRPIX2};
$crval2=$h->{CRVAL2};
$cdelt2=$h->{CDELT2};
$crpix3=$h->{CRPIX3};
$crval3=$h->{CRVAL3};
$cdelt3=$h->{CDELT3};


print OUT "S $pix $pix 1\n";
($nx,$ny,$nz)=$in_pdl->dims;
$npx=$nx*$ny;
$out_pdl=zeroes($nz,$npx);
#for ($k=0;$k<$nz;$k++) {    
print "$nx,$ny,$nz\n";
$k=0;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
#	$px=$j+($nx-1-$i)*$ny;
	$px=$i+($ny-1-$j)*$nx;
        $t = $out_pdl->slice(":,($px)");
	$t .= $in_pdl->slice("($i),($j),:");
#	$val=$in_pdl->at($i,$j,$k);
#	my $t = $in_pdl->slice("($i)

#	set($out_pdl,$k,$px,$val);
	if ($k==0) {
	    $x=$crval1+$cdelt1*(($crpix1-1)-$j);
	    $y=$crval2+$cdelt2*($i-($crpix2-1));
	    $id=$i+$j*$nx+1;
	    $x=$j;
	    $y=$i;

	    print OUT "$id $y $x 1\n";

	}
	
    }
}
#}
close(OUT);
$$h{NAXIS}=2;
#$out_pdl->sethdr($h);
$out_pdl->hdr->{CRPIX1}=$crpix3; 
$out_pdl->hdr->{CRVAL1}=$crval3; 
$out_pdl->hdr->{CDELT1}=$cdelt3; 
$out_pdl->wfits($output_rss);

exit;
