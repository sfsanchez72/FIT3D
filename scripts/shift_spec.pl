#!/usr/bin/perl
#
# This program reads a 1d ascii spectrum
# and interpolated it to change the wavelength solution
#
#
#

if ($#ARGV<4) {
    print "USE: shift_spec.pl SPEC.txt OUT_SPEC.txt CRVAL CDELT NPIX\n";
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$crval=$ARGV[2];
$cdelt=$ARGV[3];
$npix=$ARGV[4];


$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $w[$n]=$data[1];
    $f[$n]=$data[2];
    $n++;
}
close(FH);
