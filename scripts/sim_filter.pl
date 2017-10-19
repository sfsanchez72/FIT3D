#!/usr/bin/perl


if ($#ARGV<2) {
    print "USE: sim_filter.pl CENTRAL_WAVE   FWHM trans.txt\n";
    exit;
}

$wcen=$ARGV[0];
$fwhm=$ARGV[1];
$output=$ARGV[2];


$wmin=$wcen-2*$fwhm;
$wmax=$wcen+2*$fwhm;

$w=$wmin;
$s=$fwhm/2.45;
open(FH,">$output");
while($w<$wmax) {
    $f=exp(-0.5*(($w-$wcen)/$s)**2);
print FH "$w $f\n";

$w=$w+30;
}

close(FH);

exit;
