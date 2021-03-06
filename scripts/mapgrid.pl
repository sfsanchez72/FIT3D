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




if ($#ARGV<2) {
    print "USE: mapgrid.pl slice.txt out.model PREFIX_OUT [OVERRIDE]\n";
    exit;
}

$slice=$ARGV[0];
$out_model=$ARGV[1];
$prefix=$ARGV[2];
$over=0;
if ($#ARGV==4) {
    $over=$ARGV[3];
} 


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
	$dw[$n_row][$i]=2.345*$data[5];
	$e_dw[$n_row][$i]=2.345*$data[6];
	$f[$n_row][$i]=2.51*$data[5]*$data[3];
	$e_f[$n_row][$i]=2.51*$data[6]*$data[3]+2.51*$data[5]*$data[4];

	if ($model[$n_row][$i] ne "poly1d") {
	    $q[$n_row][$i]=1;
	    $w_clean[$i][$n_clean[$i]]=$w[$n_row][$i];
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

for ($i=0;$i<($n_mod-1);$i++) {
    print "$n_clean[$i] Clean final data for system $i\n";    
    my @x_stat;
    for ($j=0;$j<$n_clean[$i];$j++) {
	$x_stat[$j]=$w_clean[$i][$j];
    }
    $stat = Math::Stat->new(\@x_stat);
    $mean = $stat->average();
    $sig = $stat->stddev();
    print "Sistemic Wavelength of system $i =$mean+-$sig\n";
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
	$outfile_w=$prefix."_wave_".$num.".fits";
	$outfile_ew=$prefix."_ewave_".$num.".fits";
	$outfile_dw=$prefix."_fwhm_".$num.".fits";
	$outfile_edw=$prefix."_efwhm_".$num.".fits";
	
	
#
# We have to create a grid data
#
	
	$ii=0;
	$jj=0;
	for ($j=0;$j<$n_row;$j++) {
	    $ii=($y[$j]-$y_min)/$dpix;
	    $jj=($x[$j]-$x_min)/$dpix;
	    
	    $af[$ii][$jj]=$f[$j][$i];
	    $ae_f[$ii][$jj]=$e_f[$j][$i];
	    $aw[$ii][$jj]=$w[$j][$i];
	    $ae_w[$ii][$jj]=$e_w[$j][$i];
	    $adw[$ii][$jj]=$dw[$j][$i];
	    $ae_dw[$ii][$jj]=$e_dw[$j][$i];
	    
	}
	system("rm $outfile_f");
	system("rm $outfile_ef");
	system("rm $outfile_w");
	system("rm $outfile_ew");
	system("rm $outfile_dw");
	system("rm $outfile_edw");
	write_img($outfile_f,$nx,$ny,\@af);
	write_wcs($outfile_f,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	write_img($outfile_ef,$nx,$ny,\@ae_f);
	write_wcs($outfile_ef,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	write_img($outfile_w,$nx,$ny,\@aw);
	write_wcs($outfile_w,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	write_img($outfile_ew,$nx,$ny,\@ae_w);
	write_wcs($outfile_ew,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	write_img($outfile_dw,$nx,$ny,\@adw);
	write_wcs($outfile_dw,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	write_img($outfile_edw,$nx,$ny,\@ae_dw);
	write_wcs($outfile_edw,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
    } else {
	$outfile=$prefix."_cont.fits";
	$ii=0;
	$jj=0;
	for ($j=0;$j<$n_row;$j++) {
	    $ii=($y[$j]-$y_min)/$dpix;
	    $jj=($x[$j]-$x_min)/$dpix;
	    
	    $acont[$ii][$jj]=$cont[$j];
	}
	system("rm $outfile");
	write_img($outfile,$nx,$ny,\@acont);
	write_wcs($outfile,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
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
