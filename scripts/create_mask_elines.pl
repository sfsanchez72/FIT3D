#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#


if ($#ARGV<3) {
    print "USE: create_mas_eline.pl elines_file redshift width mask_file\n";
    exit;
}

$efile=$ARGV[0];
$redshift=$ARGV[1];
$width=$ARGV[2];
$mask_file=$ARGV[3];

open(FH,"<$efile");
open(OUT,">$mask_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$w1=$data[0]*(1+$redshift)-0.5*$width;
	$w2=$data[0]*(1+$redshift)+0.5*$width;
	print OUT "$w1 $w2\n";
    }
}
close(OUT);
close(FH);
exit;

