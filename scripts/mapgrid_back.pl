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
use PDL::Image2D;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";




if ($#ARGV<3) {
    print "USE: mapgrid_back.pl NX NY out.model PREFIX_OUT [OVERRIDE]\n";
    exit;
}

$nx=$ARGV[0];
$ny=$ARGV[1];
$dpix=1;
$out_model=$ARGV[2];
$prefix=$ARGV[3];
$over=0;
if ($#ARGV==4) {
    $over=$ARGV[4];
} 


#$nx=abs($x_max-$x_min)/$dpix+1;
#$ny=abs($y_max-$y_min)/$dpix+1;
print "GRID=$x_min $x_max $y_min $y_max $dpix [$nx,$ny]\n";
#exit;

$n_row=0;
$N_MOD=0;
open(FH,"<$out_model");
while($line=<FH>) {
    chop($line);
    $n_mod=$line;
    if ($n_mod>$N_MOD) {
	$N_MOD=$n_mod;
    }
    for ($i=0;$i<$n_mod;$i++) {
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$model[$n_row][$i]=$data[0];
	$w[$n_row][$i]=$data[1];
	$e_w[$n_row][$i]=$data[2];
	$f[$n_row][$i]=$data[3];
	$e_f[$n_row][$i]=$data[4];
	$dw[$n_row][$i]=2.345*($data[5]);#*300000;
	$e_dw[$n_row][$i]=2.345*($data[6]);#*300000;
	$v[$n_row][$i]=$data[7];
	$e_v[$n_row][$i]=$data[8];

	$AV[$n_row][$i]=$data[9];
	$e_AV[$n_row][$i]=$data[10];

#
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
#	    print "$cont[$n_row]\n";
	}	
    }
    $n_row++;

}
close(FH);

print "N_ROW=$n_row\n";
#for ($i=0;$i<$n_mod;$i++) {
for ($i=0;$i<$N_MOD;$i++) {
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

for ($i=0;$i<$N_MOD;$i++) {
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
	if ($model[0][$i] eq "back") {
	    $outfile_AV=$prefix."_AV_".$num.".fits";
	    $outfile_eAV=$prefix."_eAV_".$num.".fits";
	}
	
	
#
# We have to create a grid data
#
	
#	$ii=0;
#	$jj=0;
	$j=0;
	for ($ii=0;$ii<$nx;$ii++) {
	    for ($jj=0;$jj<$ny;$jj++) {

#	for ($j=0;$j<$n_row;$j++) {
		#$ii=($y[$j]-$y_min)/$dpix;
		#$jj=($x[$j]-$x_min)/$dpix;
#		print "$ii $jj $f[$j][$i]\n";
		
		$af[$jj][$ii]=$f[$j][$i];
		$ae_f[$jj][$ii]=$e_f[$j][$i];
		$aw[$jj][$ii]=$v[$j][$i];
		$ae_w[$jj][$ii]=$e_v[$j][$i];
		$adw[$jj][$ii]=$dw[$j][$i];
		$ae_dw[$jj][$ii]=$e_dw[$j][$i];
		$aAV[$jj][$ii]=$AV[$j][$i];
		$ae_AV[$jj][$ii]=$e_AV[$j][$i];
		$j++;
	    }
	}
	system("rm $outfile_f");
	system("rm $outfile_ef");
	system("rm $outfile_w");
	system("rm $outfile_ew");
	system("rm $outfile_dw");
	system("rm $outfile_edw");
	if ($model[0][$i] eq "back") {
	    system("rm $outfile_AV");
	    system("rm $outfile_AV");
	}
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
	if ($model[0][$i] eq "back") {
	    write_img($outfile_AV,$nx,$ny,\@aAV);
	    write_wcs($outfile_AV,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	    write_img($outfile_eAV,$nx,$ny,\@ae_AV);
	    write_wcs($outfile_eAV,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);
	}
    } else {
	$outfile=$prefix."_cont.fits";
	$ii=0;
	$jj=0;
	$j=0;
	for ($ii=0;$ii<$nx;$ii++) {
	    for ($jj=0;$jj<$ny;$jj++) {

#	for ($j=0;$j<$n_row;$j++) {
	#    $ii=($y[$j]-$y_min)/$dpix;
	#    $jj=($x[$j]-$x_min)/$dpix;
	    
	    $acont[$jj][$ii]=$cont[$j];
	    $j++;
	}
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
