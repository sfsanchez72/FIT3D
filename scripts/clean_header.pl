#!/usr/bin/perl
#
#
#


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");


if ($#ARGV<1) {
    print "clean_header.pl input output\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
print "Start Reading\n";
my $naxes=read_naxes($infile);
($nx,$ny) = @$naxes;
@array=read_img($infile);
write_img($outfile,$nx,$ny,\@array); 

exit;
