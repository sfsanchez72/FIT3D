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


if ($#ARGV<11) {
    print "USE: compare_back_list_mult_rss.pl SPEC_rss.fits VAR_rss.fits BACK_LIST.fits OUTFILE AGE_MET_OUTFILE MASK_LIST REDSHIFT SIGMA_DISP PLOT Av_min Av_max delta_Av [min max] [wmin wmax]\n";
    exit;
}

$unc_file=$ARGV[0];
$var_file=$ARGV[1];
$back_list=$ARGV[2];
$outfile=$ARGV[3];
$outfile_AM=$ARGV[4];
$mask_list=$ARGV[5];
$redshift=$ARGV[6];
$sigma=$ARGV[7];
$out_file="junk.junk";
$factor=1;
$box=1;
$plot=$ARGV[8];
$smooth=1;
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
}

#
# We read the input spectra
#

$pdl_flux_unc=rfits($unc_file);
$pdl_var_unc=rfits($var_file);
$pdl_weight=1/$pdl_var_unc;
($n_unc,$NY)=$pdl_flux_unc->dims;
$crpix_unc=$pdl_flux_unc->hdr->{CRPIX1};
$cdelt_unc=$pdl_flux_unc->hdr->{CDELT1};
$crval_unc=$pdl_flux_unc->hdr->{CRVAL1};
#$n_c=int(($wave_unc[$n_unc-1]-$crval)/$cdelt+1-$crpix);
#print "($n_unc,$NY) $crpix_unc $cdelt_unc $crval_unc\n";


#print "PASO!!!!\n";

#$n_unc=0;
$y_min=1e12;
$y_max=-1e12;
for ($j=0;$j<$n_unc;$j++) {
    $wave_unc[$j]=($crval_unc+$cdelt_unc*($j+1-$crpix));
    $masked[$j]=1;
    if ($j>0) {
	$w_test=$wave_unc[$j-1];
	for ($jj=0;$jj<$nmask;$jj++) {
	    if (($w_test>$start_mask[$jj])&&($w_test<$end_mask[$jj])) {
		$masked[$j]=0;	    
	    }
	}
	if ($def==2) {
	    if ($w_test<$min_wave) {
		$masked[$j]=0;	    
	    }
	    if ($w_test>$max_wave) {
		$masked[$j]=0;	    
	    }
	}
    }
#    print "$j/$n_unc\n";
}
$pdl_masked=pdl(@masked);
$pdl_flux_masked=$pdl_flux_unc*$pdl_masked;

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
$dpix_unc=$wave_unc[1]-$wave_unc[0];




$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};
#$n_c=int(($wave_unc[$n_unc-1]-$crval)/$cdelt+1-$crpix);
#print "($n_c,$nf) $crpix $cdelt $crval\n";



for ($j=0;$j<$n_c;$j++) {
    $wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
}
$dpix_c=$wave_c[1]-$wave_c[0];
$rsigma=$sigma/$dpix_c;
#
# We create a kernel
#
$box=int(3*$rsigma);
$kernel=zeroes(2*$box);
$norm=0;
$flux_c[$i]=0;
for ($j=0;$j<2*$box;$j++) {
    $gaus=exp(-0.5*((($i-$box)/$rsigma)**2));    
    set($kernel,$j,$gaus);
    $norm=$norm+$gaus;
}
$kernel=$kernel/$norm;

$pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;



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
#	if ($masked[$i]>0) {
#	    $error[$i][$iii]=1/abs($val);
#	    $error[$i][$iii]=10;
#	} else {
#	    $error[$i][$iii]=1e20;
#	}
    }
    print "$iii/$nf files\n";
}
print "We start the fit...\n";


$pdl_out=$pdl_flux_unc;


open(OUT,">$outfile");
open(AM,">$outfile_AM");
for ($J=0;$J<$NY;$J++) {
    $MIN_CHISQ=1e12;

    print "$J/$NY\n";
    @dims=$pdl_flux_masked->dims;
#    print "DIM=@dims\n";
    $pdl_flux_masked_now=$pdl_flux_masked->slice(":,$J");
    $pdl_flux_unc_now=$pdl_flux_unc->slice(":,$J");
    $pdl_var_now=$pdl_var_unc->slice(":,$J");
    @flux_masked=list($pdl_flux_masked_now);
    @flux_unc=list($pdl_flux_unc_now);
    @var_unc=list($pdl_var_now);

    for ($i=0;$i<$n_unc;$i++) {
	if ($flux_unc[$i] eq "nan") {
	    $flux_unc[$i]=0;
	}
	if ($flux_masked[$i] eq "nan") {
	    $flux_masked[$i]=0;
	}
    }


    $median=median(@flux_unc);
    print OUT "$J $median -1e21\n";

    $Av=$Av_min;
    $min_chi_sq=1e12;
    while ($Av<$Av_max) {
	#print "($n_unc,$nf)\n";
	$pdl_model=zeroes($n_unc,$nf);
	$pdl_error=zeroes($n_unc);
	for ($j=0;$j<$nf;$j++) {
	    for ($i=0;$i<$n_unc;$i++) {
		$wave_res=$wave_unc[$i]/(1+$redshift);
		$dust_rat=A_l(3.1,$wave_res);
		$dust=10**(-0.4*$Av*$dust_rat);  
#	    $val=($model[$i][$j]/$model[$i_scale][$j])*$flux_unc[$i_scale]*$dust;
		$val=$model[$i][$j]*$dust;
		if ($val eq "nan") {
		    $val=0;
		}
		set($pdl_model,$i,$j,$val);
	#	@DDD=$pdl_weight->dims;
	#	print "PASO $i/$n_unc (@DDD)\n";		
		$e_val=$pdl_weight->at($i,$j);
		if ($e_val eq "nan") {
		    $e_val=1e-20;
		}


	#	print "PASO $i/$n_unc\n";
		if ($masked[$i]==1) {
		    set($pdl_error,$i,$e_val);
		} else {
		    set($pdl_error,$i,1e20);
		}
		# print "$j/$nf $wave_res $val ($model[$i][0]) $e_val $flux_unc[$i]\n";
	    }
	}
#
# We fit
#
	print "We fit with Av=$Av ";
	$pdl_flux_masked_fit=pdl(@flux_masked);
	($yfit, $coeffs) = my_linfit1d $pdl_flux_masked_fit,$pdl_model,$pdl_error;
	
	
#
# We remove the models that are negative
#

	$niter=1;
	do {
	    $nf_new=0;
	    $nf_neg=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k);
		if ($C>0) {
		    $nf_new++;
		} else {
		    $nf_neg++;
		}
	    }
	    if ($nf_new>0) {
		$pdl_model_new=zeroes($n_unc,$nf_new);
		$pdl_error_new=zeroes($n_unc);
		$nf_i=0;
		for ($k=0;$k<$nf;$k++) {
		    $C=$coeffs->at($k);
		    if ($C>0) {
			for ($i=0;$i<$n_unc;$i++) {
			    $val=$pdl_model->at($i,$k);
			    set($pdl_model_new,$i,$nf_i,$val);
			    $e_val=$pdl_error->at($i);
			    set($pdl_error_new,$i,$e_val);		    
			}
			$nf_i++;
		    }
		}
		($yfit, $coeffs_new) = my_linfit1d $pdl_flux_masked_fit,$pdl_model_new,$pdl_error_new;	
		
		$nf_neg=0;
		$nf_i=0;
		for ($k=0;$k<$nf;$k++) {
		    $C=$coeffs->at($k);
		    if ($C>0) {
			$val=$coeffs_new->at($nf_i);
			set($coeffs,$k,0,$val);
			$nf_i++;
			if ($val<0) {
			    $nf_neg++;
			}
		    } else {
			set($coeffs,$k,0);
		    }
		}
		
		
	    } else {
		$nf_new=$nf;
		$pdl_model_new=$pdl_model;
		$pdl_error_new=$pdl_error;
	    }
	    $niter++;
	    print "N.C.Neg=$nf_neg NITER=$niter\n";
	} while ($nf_neg>0);
	    
	    
#
# We determine the chi**2
#
	@J=$yfit->dims;
	@K=$coeffs->dims;
#    print "DONE ( $coeffs)\n"; <stdin>;
	$chi=0;
	for ($j=0;$j<$n_unc;$j++) {
#	print "$j/$n_unc\n";
	    $out_spec[$j]=($yfit->at($j));#
		if ($flux_unc[$j]!=0) {
		    $chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/$var_unc[$j];
		}
	}
	
	
	
	
	
#	$chi_sq=(($chi/($n_unc-$nf))**0.5);
	$chi_sq=(($chi/($n_unc-$nf)));
	if ($chi_sq<$MIN_CHISQ) {
#	my @model_spec;
	    $MIN_CHISQ=$chi_sq;
	    $Av_MIN=$Av;
	    for ($j=0;$j<$n_unc;$j++) {
		$model_spec[$j]=0;
		$wave_res=$wave_unc[$j]/(1+$redshift);
		$dust_rat=A_l(3.1,$wave_res);
		$dust=10**(-0.4*$Av*$dust_rat);  
		$norm_C=0;
		for ($k=0;$k<$nf;$k++) {
		    $C=$coeffs->at($k);
		    if ($C>0) {
			$norm_C=$norm_C+$C;
		    }
		    $model_spec[$j]=$model_spec[$j]+$C*$model_no_mask[$j][$k]*$dust;
		}
		$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
#	    print "$wave_unc[$j] $flux_unc[$j] $model_spec[$j] $res_spec[$j]\n";
	    }
	    $age_min=0;
	    $met_min=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k);
		if ($C>0) {
		    $age_min=$age_min+$C*$age_mod[$k]/$norm_C;
		    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		 #   print "$k $C $norm_C $age_min $age_mod[$k] $met_min $met_mod[$k]\n";
		}
	    }
	    
	}
	$name="multi";
	$scale="1";
#    $chi_sq=($chi)**0.5;
	
	print "CHI_SQ=$chi_sq\n";
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
	    pglabel("Wavelength","Counts","$J/$NY $Av $chi_sq");
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
	$Av=$Av+$delta_Av;
    }
    close(OUT);
#print "----------------------------\n";
    print "CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av_MIN\n";
    
    
#    open(FH,">compare_back_list_dust_mult.out");
#    print FH "$MIN_CHISQ $age_min $met_min $Av_MIN\n";
#    close(FH);
    
#
# We save the residual spectrum
#

    print AM "$median AGEMET $age_min 0.0 $met_min 0.0 $Av_MIN 0.0 $redshift 0.0\n";



#    $pdl_res_spec=pdl(@res_spec);
#    $pdl_model_spec=pdl(@model_spec);


#    print "PASO??\n";
#    print "AQUI\n";

    for ($i=0;$i<$n_unc;$i++) {
#	set($pdl_res,$i,$J,$res_spec[$i]);
	set($pdl_out,$i,$J,$model_spec[$i]);
#	print "$i $model_spec[$i] $res_spec[$i]\n";
    }

#    my $t=$pdl_res->slice(":,$J");
#    $t .= $pdl_res_spec;
#    my $t=$pdl_out->slice(":,$J");
#    $t .= $pdl_model_spec;
#    print "AQUI\n";
}

close(AM);
close(OUT);



$pdl_out->wfits("model.rss.fits");

$call="rm -f res.rss.fits";
system($call);
$call="imarith ".$unc_file." - model.rss.fits res.rss.fits";
system($call);
$call="write_img_header.pl res.rss.fits CRPIX1 ".$crpix_unc;
system($call);
$call="write_img_header.pl res.rss.fits CRVAL1 ".$crval_unc;
system($call);
$call="write_img_header.pl res.rss.fits CDELT1 ".$cdelt_unc;
system($call);

    exit;
    
    
