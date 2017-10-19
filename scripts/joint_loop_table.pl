#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<2) {
    print "USE: joint_loop_table.pl auto_ssp.out coeffcs.out OUTPOUT.txt\n";
    exit;
}


$n=0;
open(FH,"<$ARGV[0]");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$chi[$n]=$data[0];
	$age[$n]=$data[1];
	$met[$n]=$data[2];
	$av[$n]=$data[3];
	$n++;
    }
}
close(FH);

$nt=0;
$nn=0;
open(FH,"<$ARGV[1]");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$LINE[$nt]=$line;
	if ($data[0]>$nn) {
	    $nn=$data[0];
	}
	$nt++;
    }
}
close(FH);
$nn=$nn+1;


open(OUT,">$ARGV[2]");
print OUT "# 1 Chi-sq\n";
print OUT "# 2 Age_log\n";
print OUT "# 3 Met_log\n";
print OUT "# 4 Av_log\n";
print OUT "# 5 Age_mean\n";
print OUT "# 6 Met_mean\n";
print OUT "# 7 Av_mean\n";
print OUT "# 8 N_SSP\n";
for ($j=0;$j<$nn;$j++) {
    $k1=9+$j;
    $k2=9+$j+1;
    $k3=9+$j+2;
    $k4=9+$j+3;
    $J=$j+1;
    print OUT "# $k1 Age_$J\n"; 
    print OUT "# $k2 Met_$J\n"; 
    print OUT "# $k3 Av_$J\n"; 
    print OUT "# $k4 Coeff_$J\n"; 
}
for ($j=0;$j<$n;$j++) {
    printf(OUT "%5.3f %5.3f %5.3f %5.3f ",$chi[$j],$age[$j],$met[$j],$av[$j]);
    $age_now=0;
    $met_now=0;
    $av_now=0;
    for ($i=1;$i<($nn+1);$i++) {
	$val=$LINE[($nn+1)*$j+$i];
	@data=split(" ",$val);
	$f_now=$data[4];
	$age_now=$age_now+$f_now*$data[1];
	$met_now=$met_now+$f_now*$data[2];
	$av_now=$av_now+$f_now*$data[6];
    }
    printf(OUT "%5.3f %5.3f %5.3f %d ",$age_now,$met_now,$av_now,$nn);
    for ($i=1;$i<($nn+1);$i++) {
	$val=$LINE[($nn+1)*$j+$i];
	@data=split(" ",$val);
	$age_now=$data[1];
	$met_now=$data[2];
	$f_now=$data[4];
	$av_now=$data[6];
#	print "$val\n";
	printf(OUT "%5.3f %5.3f %5.3f %5.3f ",$age_now,$met_now,$av_now,$f_now);
    }
    print OUT "\n";
#    print "Enter\n"; <stdin>;
}
close(OUT);



exit;
