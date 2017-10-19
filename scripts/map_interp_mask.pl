#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";




if ($#ARGV<5) {
    print "USE: map_interp_back.pl mask.txt out.model PREFIX_OUT INTERP_MODE GRID_FUNC  GRID_PIX \n";
    print "mask.txt is an 'E3D PT file' MASK=1 excluded\n";
    exit;
}

$slice=$ARGV[0];
$out_model=$ARGV[1];
$prefix=$ARGV[2];
$int_mode=$ARGV[3];
$int_opt=$ARGV[4];
$int_dpix=$ARGV[5];
$over=0;

$ns=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;
$dpix=10e10;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$mask[$ns]=$data[3];
	if ($x_min>$x[$ns]) {
	    $x_min=$x[$ns];
	}
	if ($y_min>$y[$ns]) {
	    $y_min=$y[$ns];
	}
	if ($x_max<$x[$ns]) {
	    $x_max=$x[$ns];
	}
	if ($y_max<$y[$ns]) {
	    $y_max=$y[$ns];
	}
	if ($ns>0) {
	    if (($dpix>abs($x[$ns]-$x[$ns-1]))&&(abs($x[$ns]-$x[$ns-1])!=0)) {
		$dpix=abs($x[$ns]-$x[$ns-1]);
	    }
	}
	$ns++;
    }
}
close(FH);

$nx=abs($x_max-$x_min)/$dpix+1;
$ny=abs($y_max-$y_min)/$dpix+1;
print "GRID=$x_min $x_max $y_min $y_max $dpix [$nx,$ny]\n";
#exit;

$n_row=0;
open(FH,"<$out_model");
while($line=<FH>) {
    chop($line);
    $n_mod=$line;
    for ($i=0;$i<$n_mod;$i++) {
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$model[$n_row][$i]=$data[0];
	$w[$n_row][$i]=$data[1];
	$e_w[$n_row][$i]=$data[2];
	$f[$n_row][$i]=$data[3];
	$e_f[$n_row][$i]=$data[4];
	if ($w[$n_row][$i]!=0) {
	    $dw[$n_row][$i]=2.345*($data[5]/$w[$n_row][$i])*300000;
	    $e_dw[$n_row][$i]=2.345*($data[6]/$w[$n_row][$i])*300000;
	} else {
	    $dw[$n_row][$i]=0.0;
	    $e_dw[$n_row][$i]=0.0;
	    print "WARNNING: Wrong wavelength at $n_row $i (w=$w[$n_now][$i])\n";
	}
	$v[$n_row][$i]=$data[7];
	$e_v[$n_row][$i]=$data[8];

#	print "$n_row $i $dw[$n_row][$i]\n";

	if ($model[$n_row][$i] ne "poly1d") {
	    $q[$n_row][$i]=1;
	    $v_clean[$i][$n_clean[$i]]=$v[$n_row][$i];
	    $n_clean[$i]++;
	}  else {
	    # We assume that there is at least a gaussian
	    # 
	    $cont[$n_row]=int_poly(\@data,$w[$n_row][0]); 
	    $w_clean[$i][$n_clean[$i]]=0;
	    $n_clean[$i]++;
	    $q[$n_row][$i]=1;
	}	
    }
    $n_row++;

}
close(FH);

for ($i=0;$i<$n_mod;$i++) {
    if ($model[0][$i] ne "poly1d") {
	print "$n_clean[$i] Clean final data for system $i\n";    
	my @x_stat;
	for ($j=0;$j<$n_clean[$i];$j++) {
	    $x_stat[$j]=$v_clean[$i][$j];
	}
	$stat = Math::Stat->new(\@x_stat);
	$mean = $stat->average();
	$sig = $stat->stddev();
	print "Sistemic velocity of system $i =$mean+-$sig\n";
    }
}

for ($i=0;$i<$n_mod;$i++) {
    print "$model[0][$i]\n";
    if ($model[0][$i] ne "poly1d") {

	if($i<10) {
	    $num="0".$i;
	} else {
	    $num=$i;
	}
	$outfile=$prefix."_".$num.".fits";
	
	$outfile_f=$prefix."_flux_".$num.".fits";
	$outfile_ef=$prefix."_eflux_".$num.".fits";
	$outfile_w=$prefix."_vel_".$num.".fits";
	$outfile_ew=$prefix."_evel_".$num.".fits";
	$outfile_dw=$prefix."_disp_".$num.".fits";
	$outfile_edw=$prefix."_edisp_".$num.".fits";
	
	
#
# We have to create an interpolated data
#
	
	$ii=0;
	$jj=0;

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $f[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_f");
	$call="interpol -if junk_int.tmp -of ".$outfile_f." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $e_f[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_ef");
	$call="interpol -if junk_int.tmp -of ".$outfile_ef." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $v[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_w");
	$call="interpol -if junk_int.tmp -of ".$outfile_w." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $e_v[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_ew");
	$call="interpol -if junk_int.tmp -of ".$outfile_ew." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $dw[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_dw");
	$call="interpol -if junk_int.tmp -of ".$outfile_dw." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $e_dw[$j][$i]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile_edw");
	$call="interpol -if junk_int.tmp -of ".$outfile_edw." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);



    } else {
	$outfile=$prefix."_cont.fits";
	$ii=0;
	$jj=0;

	open(JUNK,">junk_int.tmp");
	for ($j=0;$j<$n_row;$j++) {
	    if ($mask[$j]==0) {
		printf JUNK "$j $x[$j] $y[$j] $cont[$j]\n";
	    }
	}
	close(JUNK);
	system("rm $outfile");
	$call="interpol -if junk_int.tmp -of ".$outfile." -dp ".$int_dpix." -gf ".$int_mode." -go ".$int_opt." -nc 4";
	system($call);

    }

}

exit;


sub int_poly {
    my $point=$_[0];
    my $wave=$_[1];
    my $suma=0;
    my @data=@$point;
    my $j,$k;
    for ($j=1;$j<$#data;$j=$j+2) {
	$k=($j-1)/2;
	$suma+=$data[$j]*($wave**$k);
#	print "WAVE=$wave $suma $data[$j] $k\n";
    }
    return $suma;
}
