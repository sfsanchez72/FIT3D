#!/usr/bin/perl
#
#
#

use DBI;
use PGPLOT;  # Load PGPLOT module
use Statistics::OLS;
use Math::Stat;


$file=$ARGV[0];
$out=$file.".out";
if ($#ARGV<0) {
    $ps=0;
} else {
    $ps=1;
}

$ENV{'PGPLOT_ENVOPT'}="V";

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");
open(FH,"<$file");
#open(FH,"<RXSJ_May05.txt");
open(OUT,">$out");
while ($line=<FH>) {
    chop($line);
    if ($line != "#") {
	@data=split(" ",$line);
	$type=$data[5];
	$ra=$data[1];
	$dec=$data[2];
	$red=$data[6];
#	$ra =~ s/ /:/g;
#	$dec =~ s/ /:/g;
	$name="RXSJ".$data[0];
        ($hra,$junk)=split(" ",$ra);
        ($hde,$junk)=split(" ",$dec);
#	if (($hra>11)&&($hra<20)&&($hde>9)) {
	if (($hra>11)&&($hra<20)&&($hde>0)) {
#	    ($z,$alpha,$delta,$ra,$dec)=get_data_ned($name);
#	    print "$name $ra $dec $alpha $delta $z\n";
	    @ra_dat=split(" ",$ra);
	    @dec_dat=split(" ",$dec);
	    $alpha=($ra_dat[0]+$ra_dat[1]/60+$ra_dat[2]/3600)*15;
	    $delta=($dec_dat[0]+$dec_dat[1]/60+$dec_dat[2]/3600);
	    @out_data=get_data_sdss($alpha,$delta);
	    if ($#out_data>0) {
		$u=$out_data[9];
		$g=$out_data[10];
		$r=$out_data[11];
		$i=$out_data[12];
		$z=$out_data[13];
		$Bs=$g+0.217+0.419*($g-$r); 
		$Vs=$g-0.002-0.533*($g-$r); 
		$Rs=$r-0.155-0.089*($g-$r); 
		$ri=$r-$i;
		$rz=$r-$z;
		if ($ri<1.0) {
		    $RIs1=1.01*$ri+0.23;
		} else {
		    $RIs1=0.76*$ri+0.49;
		}
		if ($rzs<2.5) {
		    $RIs2=0.62*$rz+0.25;
		} else {
		    $RIs2=0.37*$rz+0.87;
		}
		$RIs=0.5*($RIs2+$RIs1);
		$BRs=$Bs-$Rs;
		$Is=$Rs-$RIs;
		
		$Bs=substr($Bs,0,5);
		$Rs=substr($Rs,0,5);
		$Is=substr($Is,0,5);
		

		print OUT "$name|$ra|$dec|$red|$type|$u|$g|$r|$i|$z|$Bs|$Rs|$Is\n";
		print "$name|$ra|$dec|$red|$type|$u|$g|$r|$i|$z|$Bs|$Rs|$Is\n";
		print "\n";
#		print "$name|$u|$g|$r|$i|$z|$Bs|$Rs|$Is|$BRs|$RIs\n";
		$n++;
	    } else {
		print OUT "$name|$ra|$dec|$red|$type\n";
	    }
	    print ".";
	    
	}
	

#    exit;
    }
}
close(OUT);
close(FH);


exit;
