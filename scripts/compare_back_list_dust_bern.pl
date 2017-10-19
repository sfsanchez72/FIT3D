#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Image2D;



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/work1/bernd/r3d_perl/my.pl";


if ($#ARGV<11) {
    print "USE: compare_back_list.pl SPEC1.txt BACK_LIST OUTFILE MASK_LIST REDSHIFT SIGMA_DISP WAVE_SCALE WAVE_WIDTH PLOT Av_min Av_max delta_Av [min max] [wmin wmax]\n";
    exit;
}

$unc_file=$ARGV[0];
$back_list=$ARGV[1];
$outfile=$ARGV[2];
$mask_list=$ARGV[3];
$redshift=$ARGV[4];
$sigma=$ARGV[5];
$wave_scale=$ARGV[6];
$wave_width=$ARGV[7];
$out_file="junk.junk";
$factor=1;
$box=1;
$plot=$ARGV[8];
$smooth=1;
$MIN_CHISQ=1e12;
$Av_min=$ARGV[9];
$Av_max=$ARGV[10];
$delta_Av=$ARGV[11];


$def=0;
if ($#ARGV==13) {
    $min=$ARGV[12];
    $max=$ARGV[13];
    $def=1;
}

if ($#ARGV==15) {
    $min=$ARGV[12];
    $max=$ARGV[13];
    $min_wave=$ARGV[14];
    $max_wave=$ARGV[15];
    $def=2;
}


$nf=0;
open(FH,"<$back_list");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$file[$nf]=$line;
	$nf++;
    }
}
close(FH);


if ($mask_list eq "none") {
    $nmask=0;
} else {
    open(FH,"<$mask_list");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$start_mask[$nmask]=$data[0];
	$end_mask[$nmask]=$data[1];
	$nmask++;
    }
    close(FH);
}




$n_unc=0;
$n_scale=0;
$tmp_flux=0;
$y_min=1e12;
$y_max=-1e12;
open(FH,"<$unc_file");
$i_scale=0;
while($line=<FH>) {
    @data=split(" ",$line);
    $index_unc[$n_unc]=$data[0];
    $wave_unc[$n_unc]=$data[1];
    $flux_unc[$n_unc]=$data[2];    
    if ($flux_unc[$n_unc]<$y_min) {
	$y_min=$flux_unc[$n_unc];
    }
    if ($flux_unc[$n_unc]>$y_max) {
	$y_max=$flux_unc[$n_unc];
    }
#    print "$n_unc $i_scale $wave_unc[$n_unc] $wave_scale\n";
    if ($n_unc>0) {
	if (($wave_unc[$n_unc-1]<=$wave_scale)&&($wave_unc[$n_unc]>$wave_scale)) {
	    $i_scale=$n_unc;	    
	}
	if (($wave_unc[$n_unc-1]<=($wave_scale+($wave_width/2.0)))&&($wave_unc[$n_unc]>($wave_scale-($wave_width/2.0)))) {
	    $n_scale++;
	    $tmp_flux=$tmp_flux+$flux_unc[$n_unc];
	#   print "$wave_unc[$n_unc] $tmp_flux $n_scale\n";
	}
    }
    $masked[$n_unc]=1;
    $w_test=$wave_unc[$n_unc-1];
    for ($j=0;$j<$nmask;$j++) {
	if (($w_test>$start_mask[$j])&&($w_test<$end_mask[$j])) {
	    $masked[$n_unc]=0;	    
	}
    }
    
    if ($def==2) {
	if ($w_test<$min_wave) {
	    $masked[$n_unc]=0;	    
	}
	if ($w_test>$max_wave) {
	    $masked[$n_unc]=0;	    
	}
    }

    $flux_masked[$n_unc]=$flux_unc[$n_unc]*$masked[$n_unc];


    $n_unc++;
}
$scale_flux=$tmp_flux/$n_scale;
#print "$scale_flux\n";
#print "$i_scale $wave_scale $flux_unc[$i_scale]\n";
close(FH);
if ($def==2) {
    $y_min=$min;
    $y_max=$max;
} else {
    $min_wave=$wave_unc[0];
    $max_wave=$wave_unc[$n_unc-1];
}

if ($def==1) {
    $y_min=$min;
    $y_max=$max;
}
#
# We check the median value
#
$median=median(@flux_masked);
open(MOUT,">median.flux");
print MOUT "$median\n";
close(MOUT);


$dpix_unc=$wave_unc[1]-$wave_unc[0];

open(OUT,">$outfile");
$Av=$Av_min;
$min_chi_sq=1e12;
while ($Av<$Av_max) {

    for ($iii=0;$iii<$nf;$iii++) {
	$c_file=$file[$iii];
#    print "$c_file\n";
	
	$n_c=0;
	open(FH,"<$c_file");
	while($line=<FH>) {
	    @data=split(" ",$line);
	    $wave_c[$n_c]=$data[1]*(1+$redshift);
	    $dust_rat[$n_c]=A_l(3.1,$data[1]);	    
	    $flux_c_ini[$n_c]=$data[2];
	    $n_c++;
	}
	$dpix_c=$wave_c[1]-$wave_c[0];
	close(FH);
	
	
#
# We convolve with a gaussian
# 
#$sigma=10;
	$rsigma=$sigma/$dpix_c;
	@flux_c=@flux_c_ini;
	$box=int(3*$rsigma);
	if ($box>3) {
	    for ($i=$box;$i<($n_c-$box);$i++) {
		$norm=0;
		$flux_c[$i]=0;
		for ($j=$i-$box;$j<($i+$box);$j++) {
		    $gaus=exp(-0.5*((($i-$j)/$rsigma)**2));
		    $flux_c[$i]=$flux_c[$i]+$flux_c_ini[$j]*$gaus;
		    $norm=$norm+$gaus;
		}
		$flux_c[$i]=($flux_c[$i]/$norm);#*(10**(-0.4*$Av*$dust_rat[$i]));		
	    }
	}
	$pdl_dust=10**(-0.4*$Av*pdl(@dust_rat));

	$pdl_flux_c=pdl(@flux_c);
	$pdl_flux_c=$pdl_flux_c*$pdl_dust;
#print "$n_unc $n_c\n";
#
# We interpolate
#
#	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), pdl(@flux_c));
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	$val=$out_spec_pdl->at($i_scale);
	$scale=$scale_flux/($out_spec_pdl->at($i_scale));#*(10**(0.4*$Av*A_l(3.1,$wave_scale)));		
#	$scale=1;#$flux_unc[$i_scale]/($out_spec_pdl->at($i_scale));#*(10**(0.4*$Av*A_l(3.1,$wave_scale)));		
#	print "$scale $flux_unc[$i_scale] $val $i_scale\n";
	$out_spec_pdl=$out_spec_pdl*$scale;
	
	$min=1e12;
	$max=-1e12;
	$min_rat=1e30;
	$max_rat=-1e30;
	$chi=0;
	for ($j=0;$j<$n_unc;$j++) {
	    $out_spec[$j]=($out_spec_pdl->at($j));#
	    if ($flux_unc[$j]!=0) {
		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($flux_unc[$j]);
#		print "$j $masked[$j] $wave_unc[$j] $flux_unc[$j] $out_spec[$j] $chi\n";
	    }
    }

    @data=split(/\//,$c_file);
    $name=$data[$#data];

    $chi_sq=(($chi/$n_unc)**0.5);
#	print "$chi $n_unc $chi_sq\n";

    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
	$name_min=$name;
	$Av_MIN=$Av;
#	print "$MIN_CHISQ $name_min $Av_MIN\n";
	for ($j=0;$j<$n_unc;$j++) {
	    $out_spec[$j]=$out_spec_pdl->at($j);		
	    $res_spec[$j]=$flux_unc[$j]-$out_spec[$j];
	    $model_spec[$j]=$out_spec[$j];
	}
	
    }

#    $chi_sq=($chi)**0.5;

    
#    print "$name $chi_sq $scale\n";
    print OUT "$name $chi_sq $scale $Av\n";
    
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Wavelength","Counts","$name $Av $chi_sq");
	pgsci(1);
	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(3);
	pgline($n_unc,\@wave_unc,\@flux_masked);    
	pgsci(1);
	pgclose;
	pgend;
#    print "Press Enter"; <stdin>;
    }
    }    
    $Av=$Av+$delta_Av;
}
close(OUT);
#print "----------------------------\n";
print "$MIN_CHISQ $name_min $Av_MIN\n";

$name_min =~ s/spec_ssp_//;
($AGE,$MET)=split("_",$name_min);
if ($AGE =~ "Myr") {
    $age=$AGE;
    $age =~ s/Myr//;
    $age=$age/1000;
} else {
    $age=$AGE;
    $age =~ s/Gyr//;
}
$met=$MET;
$met =~ s/z/0\./;

open(FH,">compare_back_list_dust.out");
print FH "$MIN_CHISQ $age $met $Av_MIN\n";
close(FH);

#
# We save the residual spectrum
#

open(OUT,">res_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    print OUT "$index_unc[$j] $wave_unc[$j] $res_spec[$j]\n";
}
close(OUT);

open(OUT,">model_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    print OUT "$index_unc[$j] $wave_unc[$j] $model_spec[$j]\n";
}
close(OUT);

exit;

