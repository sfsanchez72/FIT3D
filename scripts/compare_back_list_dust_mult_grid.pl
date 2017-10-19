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
use  PDL::Fit::Linfit;



use PDL::Core;
use PDL::Basic;
use PDL::Exporter;
@ISA    = qw( PDL::Exporter );
use PDL::Options ':Func';
use PDL::Slatec; # For matinv()




sub my_linfit1d {
   my $opthash = ref($_[-1]) eq "HASH" ? pop(@_) : {} ; 
   my %opt = parse( { Weights=>ones(1) }, $opthash ) ;
   barf "Usage: linfit1d incorrect args\n" if $#_<1 or $#_ > 3;
   my ($y, $fitfuncs, $wt) = @_;
   if ($#_ == 1) {
      ($y, $fitfuncs) = @_;
      $x = $y->xvals;
   }
   
#   my $wt = $opt{Weights};
   
   # Internally normalise data
   
   my $ymean = (abs($y)->sum)/($y->nelem);
   $ymean = 1 if $ymean == 0;
   my $y2 = $y / $ymean;
   
   # Do the fit
      
   my $M = $fitfuncs->xchg(0,1);
   my $C = $M->xchg(0,1) x ($M * $wt->dummy(0)) ;
   my $Y = $M->xchg(0,1) x ($y2->dummy(0) * $wt->dummy(0));

   # Fitted coefficients vector

   $a = matinv($C) x $Y;
   
   # Fitted data

   $yfit = ($M x $a)->clump(2); # Remove first dim=1
   
   $yfit *= $ymean; # Un-normalise
   if (wantarray) {
      my $coeff = $a->clump(2);
      $coeff *= $ymean; # Un-normalise
      return ($yfit, $coeff);
   }
   else{
      return $yfit;
   }  
   
}
*my_linfit1d = \&my_linfit1d;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<10) {
    print "USE: compare_back_list_mult.pl SPEC1.txt BACK_LIST.fits OUTFILE MASK_LIST REDSHIFT D_red SIGMA_DISP D_sig AV D_aV PLOT [min max] [wmin wmax] [redshift_elines_to_mask]\n";
    exit;
}

$unc_file=$ARGV[0];
$back_list=$ARGV[1];
$outfile=$ARGV[2];
$mask_list=$ARGV[3];
$redshift=$ARGV[4];
$d_redshift=$ARGV[5];
$sigma=$ARGV[6];
$d_sigma=$ARGV[7];
$Av=$ARGV[8];
$d_Av=$ARGV[9];
$plot=$ARGV[10];

$smooth=1;
$MIN_CHISQ=1e12;


$out_file="junk.junk";
$factor=1;
$box=1;


$def=0;
if ($#ARGV==12) {
    $min=$ARGV[11];
    $max=$ARGV[12];
    $def=1;
}

if ($#ARGV==14) {
    $min=$ARGV[11];
    $max=$ARGV[12];
    $min_wave=$ARGV[13];
    $max_wave=$ARGV[14];
    $def=2;
}

if ($#ARGV==15) {
    $min=$ARGV[11];
    $max=$ARGV[12];
    $min_wave=$ARGV[13];
    $max_wave=$ARGV[14];
    $elines_mask=$ARGV[15];
    $def=2;
}




#

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
    $n_mask_org=$nmask;
    open(FH,"<$elines_mask");
    while($line=<FH>) {
	chop($line);
	if ($line !~ "#") {
	    @data=split(" ",$line);
	    $w_eline[$nline]=$data[0];
	    $nline++;
	}
    }
    close(FH);
}

#
# We read the input spectrum
#
$n_unc=0;
$y_min=1e12;
$y_max=-1e12;
open(FH,"<$unc_file");
$i_scale=0;
while($line=<FH>) {
    @data=split(" ",$line);
    $index_unc[$n_unc]=$data[0];
    $wave_unc[$n_unc]=$data[1];
    $flux_unc[$n_unc]=$data[2];
    if ($flux_unc[$n_unc] eq "nan") {
	$flux_unc[$n_unc]=$flux_unc[$n_unc-1];
    }
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


$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};




#
# We create a kernel
#


print "We start the fit...\n";

open(OUT,">$outfile");

$min_chi_sq=FIT($redshift,$sigma,$Av);
print "($redshift,$sigma,$Av) $min_chi_sq\n";
print "----------------------\n";
$delta_chi=10;
#while (($min_chi_sq>1)&&($delta_chi>0.1)) {
while ($delta_chi>0.01) {

    if ($d_redshift!=0) {
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift+$d_redshift,$sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_redshift=-$d_redshift;
	}
	print "($redshift,$sigma,$Av) $min_chi_sq\n";
	do {
	    $redshift=$redshift+$d_redshift;
	    for ($iter=0;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift+($iter)*$d_redshift,$sigma,$Av);
		$R=$redshift+$iter*$d_redshift;
	    }
	    
	} while ($chi_now[2]<$chi_now[1]);
	
	$redshift=($redshift+2*$d_redshift)-$d_redshift*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);

    }

    if ($d_sigma!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma+$d_sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_sigma=-$d_sigma;
	}
	print "($redshift,$sigma,$Av) $min_chi_sq\n";
	
	do {
	    $sigma=$sigma+$d_sigma;
	    for ($iter=0;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma+($iter)*$d_sigma,$Av);
	    }
	} while ($chi_now[2]<$chi_now[1]);
	
	$sigma=($sigma+2*$d_sigma)-$d_sigma*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);
	
    }
    if ($d_Av!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma,$Av+$d_Av);
	if ($chi_test>$min_chi_sq) {
	    $d_Av=-$d_Av;
	}
	print "($redshift,$sigma,$Av) $min_chi_sq\n";
	
	do {
	    $Av=$Av+$d_Av;
	    for ($iter=0;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma,$Av+($iter)*$d_Av);
	    }
	    
	} while ($chi_now[2]<$chi_now[1]);
	$Av=($Av+2*$d_Av)-$d_Av*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);
	
    }
    $min_chi_sq_end=FIT($redshift,$sigma,$Av);
    $delta_chi=abs(($min_chi_sq_end-$min_chi_sq)/$min_chi_sq);
    print "($redshift,$sigma,$Av) $min_chi_sq\n";

    $d_redshift=$d_redshift/5;
    $d_sigma=$d_sigma/2;
    $d_Av=$d_Av/2;
    #    print "END****\n";
}
close(OUT);
#print "----------------------------\n";
print "CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av REDSHIFT=$redshift SIGMA_DISP=$sigma\n";


open(FH,">compare_back_list_dust_mult.out");
print FH "$MIN_CHISQ $age_min $met_min $Av $redshift $sigma\n";
close(FH);

$call="cp compare_back_list_dust_mult.out ".$outfile;
system($call);

#
# We save the residual spectrum
#

open(OUT,">org_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    print OUT "$index_unc[$j] $wave_unc[$j] $flux_unc[$j]\n";
}
close(OUT);

open(OUT,">res_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    print OUT "$index_unc[$j] $wave_unc[$j] $res_spec[$j]\n";
}
close(OUT);

open(OUT,">model_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    print OUT "$index_unc[$j] $wave_unc[$j] $model_spec_min[$j]\n";
}
close(OUT);

exit;


sub FIT {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av=$_[2];


#####################


#############




for ($j=0;$j<$n_c;$j++) {
    $wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
}
$dpix_c=$wave_c[1]-$wave_c[0];

$rsigma=$sigma/$dpix_c;

for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;
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

    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
#    print "$name[$iii] $age_mod[$iii] $met_mod[$iii]\n";


$box=int(3*$rsigma);
$kernel=zeroes(2*$box);
$norm=0;
$flux_c[$i]=0;
for ($j=0;$j<2*$box;$j++) {
    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
    set($kernel,$j,$gaus);
    $norm=$norm+$gaus;
}
$kernel=$kernel/$norm;

$pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;


    $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
    my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
    ($n_c_out)=$out_spec_pdl->dims;
#
# We store the result!
#
    for ($i=0;$i<$n_unc;$i++) {
#	print "$i/$n_c\n";
	$val=$out_spec_pdl->at($i);
	if ($val eq "nan") {
	    $val=0;
	}
	$model[$i][$iii]=$val*$masked[$i];
	$model_no_mask[$i][$iii]=$val;#*$masked[$i];
	if ($masked[$i]>0) {
	    $error[$i][$iii]=1/abs($val);
	  
#	    $error[$i][$iii]=10;
	} else {
	    $error[$i][$iii]=1e20;
	}

	for ($j=0;$j<$nline;$j++) {
	    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<5) {
		$error[$i][$iii]=1e10;
	    }
	}


    }
 #   print "$iii/$nf files\n";
}

    $pdl_model=zeroes($n_unc,$nf);
    $pdl_error=zeroes($n_unc,$nf);
    
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av*$dust_rat);  
	    $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    $e_val=$error[$i][$j];
	    set($pdl_error,$i,$j,1);
	}
    }
#
# We fit
#
#    print "We fit with Av=$Av ";
    $pdl_flux_masked=pdl(@flux_masked);
    ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;


#
# We remove the models that are negative
#
    $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	$C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

#    print "\n";
    if ($nf_new>0) {
	while ($nf_neg>0) {
#	    print "$nf_new $nf_neg ";
	    my $pdl_model_new=zeroes($n_unc,$nf_new);
	    my $pdl_error_new=zeroes($n_unc,$nf_new);
	    $nf_i=0;
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print " $C ";
		if ($C>0) {
		    for ($i=0;$i<$n_unc;$i++) {
			$val=$pdl_model->at($i,$k);
			set($pdl_model_new,$i,$nf_i,$val);
			$e_val=$pdl_error->at($i,$k);
			set($pdl_error_new,$i,$nf_i,$e_val);		    
		    }
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
#	    print "\n";
	    ($yfit, $coeffs_new) = my_linfit1d $pdl_flux_masked,$pdl_model_new,$pdl_error_new;	
#	    print "FIT $coeffs_new\n";
	    
	    $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print "$C ";
		if ($C>0) {
		    $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
#		    print "$val , ";
		}
	    }
#	    print "\n";
	}
#	    print "$nf_new $nf_neg ";
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print " $C ";
	    }
#	print "\n";

    } else {
	$nf_new=$nf;
	$pdl_model_new=$pdl_model;
	$pdl_error_new=$pdl_error;
    }



#
# We determine the chi**2
#
    @J=$yfit->dims;
    @K=$coeffs->dims;
#    print "DONE ( $coeffs)\n"; <stdin>;
    $chi=0;
    for ($j=0;$j<$n_unc;$j++) {
#	print "$j/$n_unc\n";
	$out_spec[$j]=($yfit->at($j,0));#
	    if ($flux_unc[$j]!=0) {
		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($flux_unc[$j]);
	    }
    }





    $chi_sq=(($chi/($n_unc-$nf))**0.5);
    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
	$Av_MIN=$Av;
	for ($j=0;$j<$n_unc;$j++) {
	    $model_spec[$j]=0;
	    $wave_res=$wave_unc[$j]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av*$dust_rat);  
	    $norm_C=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		$norm_C=$norm_C+$C;
		$model_spec[$j]=$model_spec[$j]+$C*$model_no_mask[$j][$k]*$dust;
		$model_spec_min[$j]=$model_spec[$j];
	    }
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
#	    print "$wave_unc[$j] $flux_unc[$j] $model_spec[$j] $res_spec[$j]\n";
	}
	$age_min=0;
	$met_min=0;
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    $age_min=$age_min+$C*$age_mod[$k]/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
#	    print "$k $C $age_min $age_mod[$k] $met_min $met_mod[$k]\n";
	}
	
    }
    $name="multi";
    $scale="1";
#    $chi_sq=($chi)**0.5;

#    print "CHI_SQ=$chi_sq\n";
#    print "$name $chi_sq $scale\n";
#    print OUT "$name $chi_sq $scale $Av \n";
    
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Wavelength","Counts","$name $Av $chi_sq $age_min $met_min");
	pgsci(1);
	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(8);
	pgline($n_unc,\@wave_unc,\@model_spec_min);    
	pgsci(3);
	pgline($n_unc,\@wave_unc,\@flux_masked);    
	pgsci(1);
	pgclose;
	pgend;
#    print "Press Enter"; <stdin>;
    }






    return $chi_sq;
}
