#!/usr/bin/perl

if ($#ARGV<3) {
    print "USE: write_ws.pl FITSFILE CRVAL CDELT CRPIX\n";
    exit;
}
$file=$ARGV[0];
$crval=$ARGV[1];
$cdelt=$ARGV[2];
$crpix=$ARGV[3];


$call="write_img_header.pl ".$file." CRVAL1 ".$crval;
system($call);
$call="write_img_header.pl ".$file." CRPIX1 ".$crpix;
system($call);
$call="write_img_header.pl ".$file." CDELT1 ".$cdelt;
system($call);

exit;
