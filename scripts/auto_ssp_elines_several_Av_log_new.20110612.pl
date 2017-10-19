#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#f
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


$vel_light=299792.458;
$red_elines=0.0;

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

   my   $a = matinv($C) x $Y;
   
   # Fitted data

   my $yfit = ($M x $a)->clump(2); # Remove first dim=1
   
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


if ($#ARGV<5) {
    print "USE: auto_ssp.pl SPEC1.txt BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
    print "CONFIG_FILE:\n";
    print "redshift delta_redshift min_redshift max_redshift\n";
    print "sigma delta_sigma min_sigma max_sigma\n";
    print "Av delta_Av min_Av max_Av [Same range for all]\n";
    print "N_SYSTEMS\n";
    print "(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "...\n";
    print "(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX\n";
    print "start_w_peak end_w_peak\n";
    print "wavelength_to_norm width_AA new_back_templates.fits\n";


    exit;
}

$unc_file=$ARGV[0];
$clean_file="clean_".$ARGV[0];


$back_list=$ARGV[1];
$outfile=$ARGV[2];
$out_elines="elines_".$outfile;
$out_single="single_".$outfile;
$out_fit="fit_".$outfile;
$D_SYS_VEL=300;



$mask_list=$ARGV[3];
$config_file=$ARGV[4];
$plot=$ARGV[5];

if ($plot==2) {
    $plot=1;
    $dev_plot=$outfile.".ps/CPS";
    $dev_plot_single="single_".$outfile.".ps/CPS";
} else {
    $dev_plot="/xs";
    $dev_plot_single="/xs";
}

$smooth=1;
$MIN_CHISQ=1e12;


$out_file="junk.junk";
$factor=1;
$box=1;


$def=0;
if ($#ARGV==7) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $def=1;
}

if ($#ARGV==9) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $def=2;
}

if ($#ARGV==10) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $def=2;
}

$input_redshift=0;
if ($#ARGV==11) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $redshift=$input_redshift;
    $def=2;
}



if ($#ARGV==14) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
    $def=2;
}
#print "$#ARGV $d_redshift\n"; 
open(FH,"<$config_file");
$line=<FH>;
chop($line);
($redshift,$d_redshift,$min_redshift,$max_redshift,$DV,$RV,$DS,$RS,$MIN_W,$MAX_W)=split(" ",$line);


$line=<FH>;
chop($line);
($sigma,$d_sigma,$min_sigma,$max_sigma)=split(" ",$line);

if ($#ARGV==18) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sigma=$ARGV[15];
    $d_sigma=$ARGV[16];
    $min_sigma=$ARGV[17];
    $max_sigma=$ARGV[18];
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
    $def=2;
}




$line=<FH>;
chop($line);
($Av_IN,$d_Av_IN,$min_Av,$max_Av)=split(" ",$line);
if ($#ARGV==22) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sigma=$ARGV[15];
    $d_sigma=$ARGV[16];
    $min_sigma=$ARGV[17];
    $max_sigma=$ARGV[18];
    $Av_IN=$ARGV[19];
    $d_Av_IN=$ARGV[20];
    $min_Av=$ARGV[21];
    $max_Av=$ARGV[22];
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
    $def=2;
}



if ($input_redshift>0) {
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
}

if ($d_redshift>0) {
    $D_SYS_VEL=0.5*($input_max_redshift-$input_min_redshift)*300000;
} else {
    $D_SYS_VEL=0;
}


$REDSHIFT=$redshift;
$Av_ini=$Av_IN;

if ($d_redshift!=0) {
    $fit_redshift=1;
} else {
    $fit_redshift=0;
}

print "FIT_RED $fit_redshift $d_redshift $#ARGV\n"; #<stdin>;



$line=<FH>;
chop($line);
$ns=$line;
$start_w_min=1e12;
$end_w_max=-1e12;
for ($i=0;$i<$ns;$i++) {
    $line=<FH>;
    chop($line);
    ($start_w_e,$end_w_e,$mask_e,$config_e,$npoly_e,$mask_poly_e,$nmin_e,$nmax_e)=split(" ",$line);
    $start_w_E[$i]=$start_w_e;
    $end_w_E[$i]=$end_w_e;
    $mask_E[$i]=$mask_e;
    $config_E[$i]=$config_e;
    $npoly_E[$i]=$npoly_e;
    $mask_poly_E[$i]=$mask_poly_e;
    $nmin_E[$i]=$nmin_e;
    $nmax_E[$i]=$nmax_e;
    if ($start_w_e<$start_w_min) {
	$start_w_min=$start_w_e;
    }
    if ($end_w_e>$end_w_max) {
	$end_w_max=$end_w_e;
    }
}
$line=<FH>;
chop($line);
($MIN_DELTA_CHISQ,$MAX_NITER,$CUT_MEDIAN_FLUX)=split(" ",$line);
$ABS_MIN=0.5*$CUT_MEDIAN_FLUX;
#print "$CUT_MEDIAN_FLUX\n"; exit;



$line=<FH>;
chop($line);
($start_w_peak,$end_w_peak)=split(" ",$line);
$line=<FH>;
chop($line);
($wave_norm,$w_wave_norm,$new_back_file)=split(" ",$line);
close(FH);

#print "F=$start_w_peak,$end_w_peak\n"; exit;
#$call="cp ".$config_e." tmp_e.config";
#system($call);


$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};
for ($i=0;$i<$nf;$i++) {
    $Av[$i]=$Av_IN;
    $d_Av[$i]=$d_Av_IN;
}


if ($wave_norm>0) {
    $pdl_flux_c_ini_new=rfits($new_back_file);
    ($n_c_new,$nf_new)=$pdl_flux_c_ini_new->dims;
    $crpix_new=$pdl_flux_c_ini_new->hdr->{CRPIX1};
    $cdelt_new=$pdl_flux_c_ini_new->hdr->{CDELT1};
    $crval_new=$pdl_flux_c_ini_new->hdr->{CRVAL1};
    if ($n_c_new==0) {
	print "STOP. No Templates to fit\n";
	exit;
    } else {
	print "N.Templates (1) = $nf\n";
	print "N.Templates (2) = $nf_new\n";
    }
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
#  We reorder the input spectrum:
#
$call="sort -n  ".$unc_file." > spec_auto.input";
#print "$call\n"; exit;
#$call="cp ".$unc_file." spec_auto.input";
system($call);
#exit;

#
#  We read the input spectrum
#

$n_unc=0;
$y_min=1e12;
$y_max=-1e12;
#open(FH,"<$unc_file");
open(FH,"<spec_auto.input");
$i_scale=0;
$FLUX=0;
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
    $index_unc[$n_unc]=$data[0];
    $wave_unc[$n_unc]=$data[1];
    $flux_unc[$n_unc]=$data[2];
    if ($#data>2) {
	$e_flux_unc[$n_unc]=sqrt(abs($data[3]));
	$color_unc[$n_unc]=$data[4];
    } else {
	$e_flux_unc[$n_unc]=0.0;
	$color_unc[$n_unc]=1;
    }
    if ($flux_unc[$n_unc] eq "nan") {
	$flux_unc[$n_unc]=$flux_unc[$n_unc-1];
    }
    if ($flux_unc[$n_unc]<$y_min) {
	$y_min=$flux_unc[$n_unc];
    }
    if ($flux_unc[$n_unc]>$y_max) {
	$y_max=$flux_unc[$n_unc];
    }
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

#    print "$masked[$n_unc] $w_test $min_wave $max_wave\n";


#    print "$n_unc $flux_unc[$n_unc] $e_flux_unc[$n_unc]\n";
    $flux_masked[$n_unc]=$flux_unc[$n_unc]*$masked[$n_unc];
    if (($masked[$n_unc]==1)&&($flux_unc[$n_unc]!=0)) {
	if ($e_flux_unc[$n_unc]==0) {
	    $e_flux_unc[$n_unc]=sqrt(abs($flux_unc[$n_unc]));
	}
    } else {
	$e_flux_unc[$n_unc]=1e-12;
    }


#    print "$n_unc $flux_unc[$n_unc] $masked[$n_unc] $flux_masked[$n_unc] '$min_wave' '$max_wave' $w_test $e_flux_unc[$n_unc]\n";
    if (($wave_unc[$n_unc]>$min_wave)&&($wave_unc[$n_unc]<$max_wave)) {
	$FLUX=$FLUX+$flux_masked[$n_unc];
    }
    $n_unc++;
    }
}
close(FH);
#exit;

if ($def==2) {
    $y_min=$min;
    $y_max=$max;
} else {
    ($min_wave,$max_wave)=minmax(@wave_unc);
#    $min_wave=$wave_unc[0];
#    $max_wave=$wave_unc[$n_unc-1];
}

#print "$def $min_wave $max_wave\n";
#exit;

if ($def==1) {
    $y_min=$min;
    $y_max=$max;
}
#
# We check the median value
#
$median_flux=median(@flux_masked);
open(MOUT,">median.flux");
print MOUT "$median_flux\n";
close(MOUT);
$dpix_unc=$wave_unc[1]-$wave_unc[0];

$max=3*$median_flux;



#
# We create a kernel
#


#print "We start the fit...\n";

open(OUT,">$outfile");

$med_flux=median(@flux_unc);
$min_chi_sq=FIT($redshift,$sigma,\@Av);

#print "MIN=$min_chi_sq\n"; exit;

print "CUT = $med_flux $CUT_MEDIAN_FLUX\n";
print "$redshift,$sigma\n";
if ($med_flux<$ABS_MIN) {
    $call="rm ".$out_elines." > junk.junk";
    system($call);
    $call="rm ".$outfile." > junk.junk";
    system($call);
    exit;
}
if ($med_flux>$CUT_MEDIAN_FLUX) {

    
    
    
    if ($MIN_W==0) {
	$MIN_W=$min_wave;
    }
    
    if ($MAX_W==0) {
	$MAX_W=$max_wave;
    }
    
#print "$MIN_W,$MAX_W\n"; exit;

################
# REDSHIFT DETERMINATION



    $K=0;
    open(CHI,">chi_r.txt");
    $nr=0;
    my @chi_r;
    my @red_r;
    if ($d_redshift>0) {
	
	
	$RED=$min_redshift;
	do {
	    $RED=$RED+$d_redshift;
	    $chi_now=FIT_Z_SIGMA($RED,$sigma,\@Av);
	    $chi_r[$nr]=$chi_now;
	    $red_r[$nr]=$RED;
	    
	    
	    if (($nr>1)&&($chi_r[$nr-1]<$chi_r[$nr-2])&&($chi_r[$nr-1]<$chi_r[$nr])&&($chi_r[$nr-1]<$min_chi_sq)) {
		$redshift=$red_r[$nr-1];	
		$min_chi_sq=$chi_r[$nr-1];	
		$K=$nr-1;
#		$RED=$max_redshift;
	    }
	    print CHI "$nr $chi_r[$nr] $red_r[$nr] $redshift\n";
	    
	    $nr++;
	}  while ($RED<$max_redshift);
#
# TWO
#
	$nr=0;
	$RED=$redshift-$d_redshift;
	$max_redshift=$redshift+$d_redshift;
	$d_redshift=0.2*$d_redshift;
	do {
	    $RED=$RED+$d_redshift;
	    $chi_now=FIT_Z_SIGMA($RED,$sigma,\@Av);
	    $chi_r[$nr]=$chi_now;
	    $red_r[$nr]=$RED;
	    if (($nr>1)&&($chi_r[$nr-1]<$chi_r[$nr-2])&&($chi_r[$nr-1]<$chi_r[$nr])) {
		$redshift=$red_r[$nr-1];	
		$min_chi_sq=$chi_r[$nr-1];	
		$K=$nr-1;
		$RED=$max_redshift;
	    }
	    print CHI "$nr $chi_r[$nr] $red_r[$nr] $redshift\n";
	    
	    $nr++;
	}  while ($RED<$max_redshift);

	$fit_redshift=0;
	$d_redshift=0;
    } else {
	$fit_redshift=0;
    }
    close(CHI);
    

    print "REDSHIFT = $redshift\n";	

# sigma DETERMINATION
    
    $K=0;
    open(CHI,">chi_s.txt");
    $nr=0;
    my @chi_s;
    my @sigma_s;
    if ($d_sigma>0) {

#    $sigma=$min_sigma;
	$SIGMA=$min_sigma;
	
	do {
	    $SIGMA=$SIGMA+$d_sigma;
	    $chi_now=FIT_Z_SIGMA($redshift,$SIGMA,\@Av);
	    $chi_s[$nr]=$chi_now;
	    $sigma_s[$nr]=$SIGMA;
#	if ($chi_now<$min_chi_sq) {
	    if (($nr>1)&&($chi_s[$nr-1]<$chi_s[$nr-2])&&($chi_s[$nr-1]<$chi_s[$nr])&&($chi_s[$nr-1]<$min_chi_sq)) {
		$sigma=$sigma_s[$nr-1];	
		$min_chi_sq=$chi_s[$nr-1];    
		$K=$nr;
		$SIGMA=$max_sigma;
	    }
	    $sigma=$sigma_s[$K-1];
	    print CHI "$nr $chi_s[$nr] $sigma_s[$nr] $sigma\n";
	    $nr++;
	} while ($SIGMA<$max_sigma);
	
	$SIGMA=$sigma-$d_sigma;
	$max_sigma=$sigma+$d_sigma;
	$d_sigma=0.33*$d_sigma;
	$nr=0;
	do {
	    $SIGMA=$SIGMA+$d_sigma;
	    $chi_now=FIT_Z_SIGMA($redshift,$SIGMA,\@Av);
	    $chi_s[$nr]=$chi_now;
	    $sigma_s[$nr]=$SIGMA;
#	if ($chi_now<$min_chi_sq) {
	    if (($nr>1)&&($chi_s[$nr-1]<$chi_s[$nr-2])&&($chi_s[$nr-1]<$chi_s[$nr])) {
		$sigma=$sigma_s[$nr-1];	
		$min_chi_sq=$chi_s[$nr-1];    
		$K=$nr;
		$SIGMA=$max_sigma;
	    }
	    print CHI "$nr $chi_s[$nr] $sigma_s[$nr] $sigma\n";
	    $nr++;
	} while ($SIGMA<$max_sigma);
	
	$fit_sigma=0;
	$d_sigma=0;
    } else {
	$fit_sigma=0;
    }
    close(CHI);
    print "SIGMA = $sigma\n";
} else {
#
# Below the cut!
#
    for ($i=0;$i<$nf;$i++) {
	$Av[$i]=0;
	$d_Av[$i]=0;
    }

#    $chi_now=FIT_Z_SIGMA($RED,$sigma,\@Av);
}

#$call="cat chi_s.txt";
#system($call);

#exit;



#$call="table_plot.pl chi_r.txt 0 1 'N.Iter' 'Chi' 101/xs";
#system($call);
#$call="table_plot.pl chi_r.txt 2 1 'Redshift' 'Chi' 102/xs";
#system($call);
#exit;

#
# 
#

@e_flux_unc=median_filter(2.345,\@res_spec);
for ($i=0;$i<$n_unc;$i++) {
    $e_flux_unc[$i]=2.345*abs($e_flux_unc[$i]);
}


#$redshift=$redshift-$d_redshift;
#print "REDSHIFT = $redshift\n";



my @new_flux_unc;
for ($i=0;$i<$n_unc;$i++) {
    $res_flux_unc[$i]=$flux_unc[$i]-$model_spec_min[$i];
    $new_flux_unc[$i]=$flux_unc[$i];
}
@median=median_filter(51,\@res_flux_unc);

open(CLEAN,">$clean_file");
for ($i=0;$i<$n_unc;$i++) {
    $flux_unc_clean[$i]=$flux_unc[$i]-$median[$i];
    print CLEAN "$i $wave_unc[$i] $flux_unc_clean[$i] $new_flux_unc[$i] $median[$i]\n";
}
close(CLEAN);



$fit_redshift=0;

#print "REDSHIFT = $redshift\n";
$redshift_abs=$redshift;
#exit;

#exit;

#print "($redshift,$sigma,\@Av) $min_chi_sq\n";
#print "----------------------\n";
$delta_chi=10;
#while (($min_chi_sq>1)&&($delta_chi>0.1)) {
$NITER=0;
$niter_tmp_max=10;


while (($MIN_CHISQ>$MIN_DELTA_CHISQ)&&($NITER<$MAX_NITER)) {
    if ($NITER==1) {
	$MIN_CHISQ=1e12;
    }


    if ($med_flux>$CUT_MEDIAN_FLUX) {

	$min_chi_sq=FIT($redshift,$sigma,\@Av);	
	for ($k=0;$k<$nf;$k++) {
	    my $C=$coeffs->at($k,0);
	    if ($C==0) {
		$Av[$k]=0;
		$d_Av[$k]=0;
	    }
	    if (($d_Av[$k]!=0)&&($C!=0)) {
		my @Av_test=@Av;
		$Av_test[$k]=$Av[$k]+$d_Av[$k];
#FIT_Z_SIGMA
		$chi_test=FIT($redshift,$sigma,\@Av_test);
		if ($chi_test>$min_chi_sq) {
		    $d_Av[$k]=-$d_Av[$k];
		}
		$niter_tmp=0;
		do {
		    $Av[$k]=$Av[$k]+$d_Av[$k];
		    $chi_now[0]=$min_chi_sq;
		    #print "$k 0 $chi_now[0] @Av\n";
		    for ($iter=1;$iter<3;$iter++) {
			$Av_test[$k]=$Av[$k]+$iter*$d_Av[$k];
			if ($Av_test[$k]<0) {
			    $Av_test[$k]=0.0;
			}
			$chi_now[$iter]=FIT($redshift,$sigma,\@Av_test);
		    }
		    $min_chi_sq=$chi_now[1];
		    $niter_tmp++;
		    
		    if ($Av[$k]>$max_Av) {
			$Av[$k]=$max_Av;
			$d_Av[$k]=-0.5*abs($d_Av[$k]);
		    }
		    if ($Av[$k]<$min_Av) {
			$Av[$k]=$min_Av;
			$d_Av[$k]=-0.5*abs($d_Av[$k]);
		    }
		    
		    
		} while (($chi_now[2]<($chi_now[1]-$MIN_DELTA_CHISQ))&&(($chi_now[1]-$chi_now[2])>0.01)&&($niter_tmp<$niter_tmp_max));
		
		
		if ($chi_now[2]<$MIN_CHISQ) {
		    $Av[$k]=($Av[$k]+2*$d_Av[$k])-$d_Av[$k]*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);
		if ($Av[$k]<0) {
		    $Av[$k]=0;
		}
		    
		}
		if ($Av[$k]>$max_Av) {
		    $Av[$k]=$max_Av;
		    $d_Av[$k]=-abs($d_Av[$k]);
		}
		if ($Av[$k]<$min_Av) {
		    $Av[$k]=$min_Av;
		    $d_Av[$k]=+abs($d_Av[$k]);
		}
		
		
	    }
	}
	
	
	$min_chi_sq_end=FIT($redshift,$sigma,\@Av);

	if ($min_chi_sq>0) {
	    $delta_chi=abs(($min_chi_sq_end-$min_chi_sq)/$min_chi_sq);
	}
	
    }
    $NITER++;


    $wpeak=6562;
    $Fpeak=-1e12;


    if ($ns>0) {

    $ks=0;
    $SYS_VEL=$vel_light*$REDSHIFT;
    $call="rm ".$out_elines." > junk.junk";
    my @REN;
    my @e_REN;
    system($call);
    for ($is=0;$is<$ns;$is++) {
	if ($red_elines>0) {
	    $SYS_VEL=$vel_light*$red_elines;
	    if ($is==0) {
		$SYS_VEL_MAX=$vel_light*$red_elines+$D_SYS_VEL;
		$SYS_VEL_MIN=$vel_light*$red_elines-$D_SYS_VEL;
	    } else {
		$SYS_VEL_MAX=$vel_light*$red_elines+$D_SYS_VEL;
		$SYS_VEL_MIN=$vel_light*$red_elines-$D_SYS_VEL;
	    }
	} else {
	    $SYS_VEL=$vel_light*$REDSHIFT;
	    if ($is==0) {
		$SYS_VEL_MAX=$vel_light*$REDSHIFT+$D_SYS_VEL;
		$SYS_VEL_MIN=$vel_light*$REDSHIFT-$D_SYS_VEL;
	    } else {
		$SYS_VEL_MAX=$vel_light*$REDSHIFT+$D_SYS_VEL;
		$SYS_VEL_MIN=$vel_light*$REDSHIFT-$D_SYS_VEL;
	    }
	}
	$start_w_e=$start_w_E[$is];
	$end_w_e=$end_w_E[$is];
	$mask_e=$mask_E[$is];
	$config_e=$config_E[$is];
	$npoly_e=$npoly_E[$is];
	$mask_poly_e=$mask_poly_E[$is];
	$nmin_e=$nmin_E[$is];
	$nmax_e=$nmax_E[$is];

#	$SYS_VEL_MAX=$vel_light*$redshift-10;
#	$SYS_VEL_MIN=$vel_light*$redshift+10;
	
	print "CONF=$config_e\n";
	open(IN_CONF,"<$config_e");
	open(OUT_CONF,">tmp_e.config");
	$n_conf=0;
	while($line=<IN_CONF>) {
	    chop($line);
	    if ($n_conf==2) {
		@data=split(" ",$line);
		$wpeak_res=$data[0];
		#
		# We look for the redshift
		#
		$w_min_peak=$wpeak_res*(1+$redshift)*0.95;
		$w_max_peak=$wpeak_res*(1+$redshift)*1.01;
		my $med=median(@res_spec);
		my $sig=sigma(@res_spec);
		$SYS_VEL_INI=$vel_light*$redshift;
		$dSYS_VEL=700;
		for ($i=1;$i<$n_unc-1;$i++) {
		    $W=$wave_unc[$i];
		    $F_OUT=$res_spec[$i];
		    if (($W>$w_min_peak)&&($W<$w_max_peak)&&($F_OUT>($med+2*$sig))) {
			if (($F_OUT>$res_spec[$i-1])&&($F_OUT<$res_spec[$i+1])) {
			    $SYS_VEL_NOW=($W/$wpeak_res-1)*$vel_light;
			    if (abs($SYS_VEL_NOW-$SYS_VEL_INI)<$dSYS_VEL) {
				$dSYS_VEL=abs($SYS_VEL_NOW-$SYS_VEL_INI);
				#print "NEW SYS_VEL ($wpeak_res) $SYS_VEL ->";
				#$SYS_VEL=$SYS_VEL_NOW;
				#$SYS_VEL_MIN=$SYS_VEL_NOW*0.95;
				#$SYS_VEL_MAX=$SYS_VEL_NOW*1.05;
				#my $red=$W/$wpeak_res-1;
				#print "$SYS_VEL (z=$red)\n";
			    }
			}
		    }
		    
		}

	    }
	    




	    if (($n_conf==5)) {
		print OUT_CONF "$SYS_VEL   1   $SYS_VEL_MIN  $SYS_VEL_MAX  -1\n";	    
		print  "LINE: $SYS_VEL   1   $SYS_VEL_MIN  $SYS_VEL_MAX  -1\n";	    
	    } else {
		print OUT_CONF "$line\n";
	    }
	    $n_conf++;
	}
	close(OUT_CONF);
	close(IN_CONF);
	
	open(RES,">sub_spec.txt");
	for ($i=0;$i<$n_unc;$i++) {
	    $ID=$i+1;
	    $W=$wave_unc[$i];
	    $F_OUT=$res_spec[$i];
	    print RES "$ID  $W  $F_OUT\n";
	    $NF++;
	}
	close(RES);
	

	$nmin_e=int(0.1*$n_unc);
	$nmax_e=int(0.9*$n_unc);
###############################
# Low order polynomical!
	$call="spec2img.pl sub_spec.txt res_spec.fits 1";
	system($call);
	$call="polyfit_mask.pl res_spec.fits 6 res_poly.fits ".$mask_poly_e." ".$nmin_e." ".$nmax_e." /null";
	system($call);
	$call="imarith.pl res_spec.fits - res_poly.fits new_res.fits";
	system($call);
	$call="img2spec.pl new_res.fits 0 new_res.txt";
	system($call);
	$out_fit_now=$out_fit.".".$start_w_e."_".$end_w_e.".ps/CPS";
	$box=int($sigma*6);
	
#	$call="fit_spec_back.pl sub_spec.txt none ".$start_w_e." ".$end_w_e." none 0 ".$mask_e." tmp_e.config ".$out_fit_now." > junk.junk";

	$call="fit_spec_back.pl sub_spec.txt none ".$start_w_e." ".$end_w_e." none 0 ".$mask_e." tmp_e.config ".$out_fit_now." > junk.junk";

	system($call);
	print "$call\n";
	print "DONE $is\n";
#	print "Press Enter"; <stdin>;
	$call="cat out.fit_spectra >> ".$out_elines;
	system($call);


#	<stdin>;
	#
	#
	# We read the emission lines information
	#
	open(KS,"<out.fit_spectra");
	$lineKS=<KS>;
	chop($lineKS);
	while($lineKS=<KS>) {
	    chop($lineKS);
	    if ($lineKS =~ "eline") {
		$infoKS[$ks]=$lineKS;
		@dat_now=split(" ",$lineKS);
		if ($ks==0) {
		    if (($dat_now[3]>10)&&($fit_redshift!=0)) {
#			$redshift=$dat_now[7]/$vel_light;
		    }
#		    $red_elines=$dat_now[7]/$vel_light;
		}
		$REN[$ks]=$dat_now[7]/$vel_light;
		$e_REN[$ks]=$dat_now[8]/$dat_now[7];
#		print "---- $ks $dat_now[7] $dat_now[8] $REN[$ks] $e_REN[$ks]\n";
		$ks++;
	    }
	}
	close(KS);
	print "REDSHIFT_ELINES = $red_elines\n";
	#<stdin>;
    }



#    $d_redshift=0.002;
#    $fit_redshift=0;
#    exit;


    #
    # We do the FIT with the model FIXED!
    #
    open(OUT_CONF,">tmp_e.config");
    print OUT_CONF "0 $ks 1 0.5\n";
    for ($is=0;$is<$ks;$is++) {
	print OUT_CONF "eline\n";
	@info=split(" ",$infoKS[$is]);
	print OUT_CONF "$info[1]  0    0   0   -1\n";
	print OUT_CONF "$info[3]  0    0   0   -1\n";
	print OUT_CONF "$info[5]  0    0   0   -1\n";
	print OUT_CONF "$info[7]  0    0   0   -1\n";
	print OUT_CONF "0  0    0   0   -1\n";
	print OUT_CONF "0  0    0   0   -1\n";
	print OUT_CONF "0  0    0   0   -1\n";
	print OUT_CONF "0  0    0   0   -1\n";
	print OUT_CONF "0  0    0   0   -1\n";
    }
    close(OUT_CONF);    

    $red_elines=mean(@REN);
    $s_red_elines=sigma(@REN);

    $min_e_REN=0.003;
    for ($is=0;$is<$ks;$is++) {
	if (($e_REN[$is]>0)&&($e_REN[$is]<$min_e_REN)) {
	    $red_elines=$REN[$is];
	    $min_e_REN=$e_REN[$is];
	}
#	    print "***** $is/$ks $red_elines $REN[$is] $e_REN[$is]\n";

    }


    print "mean_REN = $red_elines +- $s_red_elines\n";
    $out_fit_now=$out_fit."_all.ps/CPS";
    $call="fit_spec_back.pl sub_spec.txt none ".$start_w_min." ".$end_w_max." none 0 none tmp_e.config ".$out_fit_now." > junk.junk";
    system($call);
    print "$call\n";
    #
    # Now we subtract the emission lines
    # 

    open(MOD,"<out_mod_res.fit_spectra");
    $NF=0;
    @model_joint=@model_spec_min;
    @res_joint=@res_spec;
#    @res_clean=@res_spec;
#    @flux_clean=@flux_unc;
    $nc=0;
    $chi_joint=0;
    $ncc=0;
    while($line=<MOD>) {
	chop($line);
	@data=split(" ",$line);
	$WMOD=$data[0];
	$FMOD=$data[2];
	while ((abs($data[0]-$wave_unc[$NF])>0.1)&&($NF<$n_unc-1)) {
	    $NF++;
	} 
	$flux_now=$flux_unc[$NF]-$FMOD;	    
	$flux_masked[$NF]=$flux_now*$masked[$NF];
	$model_joint[$NF]=$model_spec_min[$NF]+$FMOD;
	$res_joint[$NF]=$res_spec[$NF]-$FMOD;
	if ($masked[$NF]==1) {
	    $flux_clean[$nc]=$flux_unc[$NF];
	    $res_clean[$nc]=$res_joint[$NF];
	    $wave_clean[$nc]=$wave_unc[$NF];
	    if ($model_joint[$NF]!=0) {
		$chi_joint=$chi_joint+(($flux_unc[$NF]-$model_joint[$NF])**2)/abs($model_joint[$NF]);	    
		$ncc++;
	    }
#	    print "$nc $wave_clean[$nc] $res_clean[$nc]\n";
	    $nc++;
	}

    }
    close(MOD);

#    print "PASO!!!!\n";



    if ($ncc!=0) {
	$chi_joint=($chi_joint/($ncc))**0.5;
    }

#    $chi2=($chi2/($NFREE))**0.5;


    }
    if ($ns>0) {
	$rms=sigma(@res_clean);
	$med_flux=median(@flux_clean);
    } else {
	$rms=sigma(@res_spec);
	$med_flux=median(@flux_unc);
	$chi_joint=$MIN_CHISQ;
    }
    open(OUT_MOD,">mod_joint_spec.txt");
    open(OUT_RES,">res_joint_spec.txt");
    for ($j=0;$j<$n_unc;$j++) {
	if ($color_unc[$j]==1) {
	    print OUT_MOD "$index_unc[$j] $wave_unc[$j] $model_joint[$j]\n";
	    print OUT_RES "$index_unc[$j] $wave_unc[$j] $res_joint[$j]\n";
	}
    }
    close(OUT_RES);
    close(OUT_MOD);

    @e_flux_unc=median_filter(2.345,\@res_joint_spec);

    for ($j=0;$j<$n_unc;$j++) {
	if ($e_flux_unc[$j]==0) {
	    $e_flux_unc[$j]=1;
	} else {
	    $e_flux_unc[$j]=10*2.345*abs($e_flux_unc[$j]);
	}
    }

#    print "Press Enter"; <sdtin>;
    #
    #


 
   print "--------------------------------------------------------------\n";
#    print "@Av\n";
    $TEMP=FIT($redshift,$sigma,\@Av);
#    print "@Av\n";
#print "MSP CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av_min REDSHIFT=$redshift SIGMA_DISP=$sigma RMS=$rms MED_FLUX=$med_flux AGE_mass=$age_min_mass MET_mass=$met_min_mass\n";

    print "MSP CHISQ=$chi_joint AGE=$age_min MET=$met_min AV=$Av_min REDSHIFT=$redshift SIGMA_DISP=$sigma RMS=$rms MED_FLUX=$med_flux AGE_mass=$age_min_mass MET_mass=$met_min_mass\n";


#    @med_flux=median_filter(11,\@res_joint);
#    for ($j=0;$j<$n_unc;$j++) {
#	$e_flux_unc[$j]=1/($median_flux+($res_joint[$j]-$med_flux[$i])**2);
#    }





    #    print "END****\n";
}
close(OUT);


#print "----------------------------\n";



open(FH,">compare_back_list_dust_mult.out");
print FH "#MIN_CHISQ age_min met_min Av redshift sigma median_FLUX redshift_ssp med_flux rms_residual age_min_mass met_min_mass SYS_VEL $back_list \n";
#print FH "$MIN_CHISQ $age_min $met_min $Av_min $redshift $sigma $FLUX $redshift_abs $med_flux $rms $age_min_mass $met_min_mass $SYS_VEL\n";
if ($chi_joint==0) {
    $chi_joint=1;
}
print FH "$chi_joint $age_min $met_min $Av_min $redshift $sigma $FLUX $redshift_abs $med_flux $rms $age_min_mass $met_min_mass $SYS_VEL\n";
close(FH);

$call="cp compare_back_list_dust_mult.out ".$outfile;
system($call);







#
# We save the residual spectrum
#

open(OUT,">org_spec.txt");
open(OUT2,">org_phot.txt");
for ($j=0;$j<$n_unc;$j++) {
    if ($color_unc[$j]==1) {
	print OUT "$index_unc[$j] $wave_unc[$j] $flux_unc[$j]\n";
    } else {
	print OUT2 "$index_unc[$j] $wave_unc[$j] $flux_unc[$j] $e_flux_unc[$j]\n";
    }
}
close(OUT2);
close(OUT);

open(OUT,">res_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    if ($color_unc[$j]==1) {
	print OUT "$index_unc[$j] $wave_unc[$j] $res_spec[$j]\n";
    }
}
close(OUT);

$NBOX=int(9*$sigma);
$N1=int(0.25*$n_unc);
$N2=int(0.75*$n_unc);
$call="smooth_spec_clip.pl res_spec.txt tmp_spec.txt ".$NBOX." /null 2 ".$N1." ".$N2;
#print "$call\n";
system($call);
$call="spec_arith.pl res_spec.txt - tmp_spec.txt gas_spec.txt";
system($call);


open(OUT,">model_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    if ($color_unc[$j]==1) {
	print OUT "$index_unc[$j] $wave_unc[$j] $model_spec_min[$j]\n";
    }
}
close(OUT);


open(OUT,">no_gas_spec.txt");
for ($j=0;$j<$n_unc;$j++) {
    my $dat=$model_spec_min[$j]+$res_joint[$j];
    if ($color_unc[$j]==1) {
	print OUT "$index_unc[$j] $wave_unc[$j] $dat\n";
    }
}
close(OUT);

    if ($wave_norm>0) {
#	print "*************** FIT_SINGLE SSP ******************\n";
	$MIN_CHISQ=1e12;
#	$plot=$plot_old;
#	print "Press Enter"; <stdin>;
	
	$min_chi_sq=FIT_SINGLE($redshift,$sigma,\@Av);
	$rms_single=sigma(@res_clean);
	$med_flux_single=median(@flux_clean);
	print "SSP CHISQ=$min_chi_sq AGE=$age_min MET=$met_min AV=@Av REDSHIFT=$redshift SIGMA_DISP=$sigma RMS=$rms_single MED_FLUX=$med_flux_single\n";
#	print "Press Enter"; <stdin>;
	open(FH,">$out_single");
	print FH "#MIN_CHISQ age_min met_min Av redshift sigma median_FLUX redshift_ssp  med_flux rms_residual age_min_mass met_min_mass $new_back_file \n";
	print FH "$min_chi_sq $age_min $met_min $Av_min $redshift $sigma $FLUX $red_abs $red_em $med_flux_single  $rms_single  $age_min $met_min\n";
	close(FH);


	system($call);
}
#    <stdin>;

    
    


$out_coeffs="coeffs.out";;
$out_all="all_".$outfile;
$call="joint_loop_table.pl ".$outfile." ".$out_coeffs." ".$out_all;
system($call);




exit;


sub FIT {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
#    print "AV FIT = '@AV' '$Av_NOW' '$_'\n";
    my $k,$j,$i,$iter;

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

#####################


#############




    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
#	    print "$j $dpix_c_val[$j] $crval $cdelt\n";
	}
	
    }
    $dpix_c_val[0]=$dpix_c_val[1];
$dpix_c=$wave_c[1]-$wave_c[0];

#$dpix_c=median($wave_c[1]-$wave_c[0];

#($dpix_c,$dpix_max)=minmax(@dpix_c_val);
#$dpix_c=median(@dpix_c_val);
#print "DPIX= $dpix_c\n"; exit;

$rsigma=$sigma/$dpix_c;

#    print "$sigma $dpix_c\n";

for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;    
    $name_min =~ s/.dat//;
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
    $header="NORM".$iii;
#    $ml[$iii]=$pdl_flux_c_ini->hdr->{$header};

 $val_ml=$pdl_flux_c_ini->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }




#    print "$name[$iii] $age_mod[$iii] $met_mod[$iii]\n";
    
    
    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    $flux_c[$i]=0;
    for ($j=0;$j<2*$box+1;$j++) {
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
	    $error[$i]=0.01*abs($e_flux_unc[$i]);
	} else {
	    $error[$i]=0;
	}


	
#	if ($NITER==0) {
	    for ($j=0;$j<$nline;$j++) {
		if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
#		    $error[$i]=0.1*$e_flux_unc[$i];
#		    $error[$i]=30;#*$e_flux_unc[$i];
		}
	    }
#	}

#		    print "$i $error[$i]\n";

#	print "$i $wave_unc[$i] $error[$i] $e_flux_unc[$i]\n";

    }
 #   print "$iii/$nf files\n";
}

    $pdl_model=zeroes($n_unc,$nf);
    $pdl_error=zeroes($n_unc,$nf);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    $e_val=$error[$i];
	    $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,$j,1/(abs($val_now)**2));
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
    my $nf_new=0;
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
#	    print "PASO $nf_new $nf_neg \n";
#	    if ($nf_new>0) {
		my $pdl_model_new=zeroes($n_unc,$nf_new);
		my $pdl_error_new=zeroes($n_unc,$nf_new);
#	    }
	    $nf_i=0;
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print " $C ";
#		print "$C $k PASO\n";
		if ($C>0) {
		    for ($i=0;$i<$n_unc;$i++) {
			$val=$pdl_model->at($i,$k);
			set($pdl_model_new,$i,$nf_i,$val);
#			print "$i PASO\n";
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
#	    print "PASO\n";
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
		if ($nf_new==0) {
		    $nf_neg=0;
		}
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



#    print "PASO\n";

#
# We determine the chi**2
#
    @J=$yfit->dims;
    @K=$coeffs->dims;
#    print "DONE ( $coeffs)\n"; <stdin>;
    $chi=0;
    $chi2=0;
    $NFREE=0;
    for ($j=0;$j<$n_unc;$j++) {
	$out_spec[$j]=($yfit->at($j,0));#
	$chi_sec[$j]=0;
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)) {
		#print "$j $flux_unc[$j]-$out_spec[$j]\n";
#0000000000000000000000000000000
#		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
		$chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
#		$chi2=$chi2+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($e_flux_unc[$j]);
		$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($out_spec[$j]);
#		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)*$error[$j];
		$NFREE++;
	    }
    }

    $chi_sq=$chi;
    if ($NFREE>0) {
	$chi_sq=($chi_sq/($NFREE))**0.5;
    }
#    $chi2=($chi2/($NFREE))**0.5;

#    $chi_sq=(($chi/($n_unc-$nf))**0.5);
#    $chi_sq=(($chi/($n_unc-1))**0.5);
    if ($chi_sq<$MIN_CHISQ) {
#    print "PASO $chi_sq $MIN_CHISQ\n";
	$MIN_CHISQ=$chi_sq;
	@Av_MIN=@Av;
	for ($j=0;$j<$n_unc;$j++) {
	    $model_spec[$j]=0;
	    $wave_res=$wave_unc[$j]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $norm_C=0;
	    $norm_C_mass=0;
	    for ($k=0;$k<$nf;$k++) {
		$dust=10**(-0.4*$Av[$k]*$dust_rat);  
		$C=$coeffs->at($k,0);
		$norm_C=$norm_C+$C;
		$norm_C_mass=$norm_C_mass+$C*$ml[$k];
		$model_spec[$j]=$model_spec[$j]+$C*$model_no_mask[$j][$k]*$dust;
		$model_spec_min[$j]=$model_spec[$j];
	    }
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
#	    print "$wave_unc[$j] $flux_unc[$j] $model_spec[$j] $res_spec[$j]\n";
	}
	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;
	open(C,">coeffs.out");
	print C "ID   AGE     MET    COEFF   Norm.Coeff  M/L   AV\n";
	print "---------------------------------------------\n";
	print "ID AGE     MET    COEFF   Norm.Coeff   M/L    AV\n";
	print "---------------------------------------------\n";
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
#	    if ($C==0) {
#		$d_Av[$k]=0;
#	    }
#	    print C "$k $age_mod[$k] $met_mod[$k] $C $CN $ml[$k] $Av[$k]\n";
	    printf(C "%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C,$CN,$ml[$k],$Av[$k]);
	    printf("%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C,$CN,$ml[$k],$Av[$k]);
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }
#		print "COEFF=$k $C $norm_C_mass $age_min $age_mod[$k] $met_min $met_mod[$k] $ml[$k] $age_min_mass $met_min_mass\n";

	}
	print "---------------------------------------------\n";
	close(C);
	$age_min=10**($age_min);
#	$met_min=10**($met_min);
	$age_min_mass=10**($age_min_mass);
#	$met_min_mass=10**($met_min_mass);
#	print "$age_min $met_min\n";
#	print "**********************************\n";

	#exit;
 #	$age_min=apr_n($age_min,4);
 #	$met_min=apr_n($met_min,6);
#	$age_min_mass=apr_n($age_min_mass,4);
#	$met_min_mass=apr_n($met_min_mass,6);
	
    }
    $name=$unc_file.", ";
    $scale="1";
#    $chi_sq=($chi)**0.5;

#    print "CHI_SQ=$chi_sq\n";
#    print "$name $chi_sq $scale\n";
#    print OUT "$name $chi_sq $scale $Av \n";
    
    if ($plot==1) {
	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
	#pglabel("Wavelength","Counts","$name $Av $chi_sq $age_min $met_min $redshift");
#	print "PP = @Av\n";
	pglabel("Wavelength","Counts","X=$chi_sq ($chi_joint) T=$age_min ($age_min_mass) Z=$met_min ($met_min_mass) Av=@Av z=$redshift sigma=$sigma");
	pgsch(0.5);
	pgsci(1);
	for ($in=0;$in<$n_unc;$in++) {
	    my $X=$wave_unc[$in];
	    my $Y=$flux_unc[$in];
	    my $eY=$e_flux_unc[$in];
	    pgsci($color_unc[$in]);    
	    my $s=16*$color_unc[$in]-15;
#    pgpoint($n,\@wave,\@flux,1);
	    #pgpoint(1,[$X],[$Y],$s);
	    #pgerrb(2,1,[$X],[$Y],[$eY],0.4);
	    #pgerrb(4,1,[$X],[$Y],[$eY],0.4);
	    pgsci(3);
	    my $Y=$flux_masked[$in];
	    if ($s<3) {
		pgpoint(1,[$X],[$Y],$s);
	    }
	}


	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(8);
	pgline($n_unc,\@wave_unc,\@model_spec_min);    
#	pgsci(3);
#	pgline($n_unc,\@wave_unc,\@flux_masked);    
	pgsci(5);
	pgline($n_unc,\@wave_unc,\@res_spec);    
	if ($nc>0) {
	    pgsci(6);
	    pgline($nc,\@wave_clean,\@res_clean);    
	}
	pgsci(1);

#	    pgsci(7);
#	    pgline($n_unc,\@wave_unc,\@error);    

	pgsci(13);
	pgsls(2);
	for ($j=0;$j<$nf;$j++) {
	    $C=$coeffs->at($j,0);
	    if ($C>0) {
		pgsci(13);
		my $model_now=$pdl_model->slice(",($j)");
		$model_now=$C*$model_now;
		@a_model_now=list($model_now);
		pgline($n_unc,\@wave_unc,\@a_model_now);    

	    }
	    if ($norm_C != 0) {
		$CN=$C/$norm_C;
	    } else {
		$CN=$C;
	    }
	    $pos_x=$min_wave*1.02;
	    $pos_y=$y_max-0.05*($j+1)*($y_max-$y_min);
	    pgsch(0.8);           # Set character height
	    if ($C>0) {
		pgsci(1);
	    } else {
		pgsci(15);
	    }
	    $CN=apr($CN);
	    pgptxt($pos_x,$pos_y,0,0,"$CN $age_mod[$j] $met_mod[$j]");
	    
#	    print "$pos_x $pos_y\n"; <stdin>;
	    pgsch(1.2);           # Set character height
	}
	pgsls(1);


	pgsci(1);


	pgclose;
	pgend;
#    print "Press Enter"; <stdin>;
    }


#
#
#




    return $chi_sq;
}


sub FIT_SINGLE {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;

    my $i,$j,$k;

#
# We renormalize the models
#


    for ($j=0;$j<$n_c_new;$j++) {
	$wave_c_new[$j]=($crval_new+$cdelt_new*($j+1-$crpix_new))*(1+$redshift);
    }
    $dpix_c=$wave_c_new[1]-$wave_c_new[0];
    $rsigma=$sigma/$dpix_c;


    my $CHI_MIN=1e12;

    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    $flux_c[$i]=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;
    
    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini_new,$kernel;
    

    for ($iii=0;$iii<$nf_new;$iii++) {
	$header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini_new->hdr->{$header};
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
#****** CONVOLUTION
	

	$pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c_new), $pdl_flux_c);
	($n_c_out)=$out_spec_pdl->dims;
	

	
	
	for ($i=0;$i<$n_unc;$i++) {
	    if ($masked[$i]>0) {
		$error[$i]=0.01/abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }
	    $val=$out_spec_pdl->at($i);
	    $model[$i][$iii]=$val;#*$masked[$i];
#	    if ($NITER==0) {
		for ($j=0;$j<$nline;$j++) {
		    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
#		    $error[$i]=(30/(1+$INTER))*$e_flux_unc[$i];
		    $error[$i]=0.1/abs($e_flux_unc[$i]);
#			$error[$i]=10*$e_flux_unc[$i];
		    }
		}
	    }
#	}
	
    }



    $pdl_model=zeroes($n_unc,$nf_new);
    $pdl_error=zeroes($n_unc,$nf_new);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf_new;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
#	    if ($j==0) {
		$wave_res=$wave_unc[$i]/(1+$redshift);
		$dust_rat=A_l(3.1,$wave_res);
		$dust=10**(-0.4*$Av[$j]*$dust_rat);  
		$a_dust[$i]=$dust;
#	    }
#	    print "DUST = $val $model[$i][$j] $a_dust[$i]\n";
	    $val=$model[$i][$j]*$a_dust[$i];
	    set($pdl_model,$i,$j,$val);
	    $e_val=$error[$i];
	    set($pdl_error,$i,$j,1);
	}
    }
    $pdl_flux_masked=pdl(@flux_masked);


#    ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
    #
    # We normalize the models
    #
 

#    my $nx1=int(($wave_norm-$w_wave_norm-$crval3)/$cdelt3);
#    my $nx2=int(($wave_norm+$w_wave_norm-$crval3)/$cdelt3);
#    my $pdl_sec_cen=$pdl_flux_masked->slice("$nx1:$nx2");
#    my $SUM=sum($pdl_sec_cen);
    
    my $pdl_chi=zeroes($nf_new);
    my $pdl_e=pdl(@error);
#    print "$pdl_e\n";
    $J_MIN=-1;
   for ($j=0;$j<$nf_new;$j++) {
       my $pdl_model_tmp=$pdl_model->slice(":,$j");
#       my $pdl_sec_cen=$pdl_model_tmp->slice("$nx1:$nx2");
#       my $SUM_NOW=sum($pdl_sec_cen);
#       my $NORM=$SUM/$SUM_NOW;
#       $pdl_model_tmp=$pdl_model_tmp*$NORM;
       ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model_tmp,$pdl_e;


       my $pdl_chi_now=$pdl_masked*(($pdl_flux_masked-$yfit)**2)*$pdl_e;
       my $chi_sq=sum($pdl_chi_now);
       $chi_sq=($chi_sq/($n_unc-1))**0.5;
       set($pdl_chi,$j,$chi_sq);
#       print "$j $CHI_MIN $chi_sq\n";

       if ($CHI_MIN>$chi_sq) {
	   $CHI_MIN=$chi_sq;
	   $J_MIN=$j;
	   $NORM_J=$NORM;


#	   @out_spec=list($pdl_model_tmp);
#	   @model_spec_min=list($pdl_model_tmp);
	   @out_spec=list($yfit);
	   @model_spec_min=list($yfit);
	   my $pdl_res=pdl(@flux_unc)-$pdl_model_tmp;
	   @res_spec=list($pdl_res);
	   $JJ=0;
	   for ($jj=0;$jj<$n_unc;$jj++) {
	       if ($masked[$jj]>0) {
		   $res_clean[$JJ]=$res_spec[$jj];
		   $flux_clean[$JJ]=$flux_unc[$jj];
		   $JJ++;
	       }
	   }




	if ($plot==1) {

	    pgbegin(0,$dev_plot_single,1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	    pgsch(1.2);           # Set character height
	    pglabel("Wavelength","Counts","SSP X=$CHI_MIN T=$age_mod[$J_MIN] Z=$met_mod[$J_MIN] z=$redshift sigma=$sigma Av=@Av");
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

   }
#    <stdin>;
#    $pdl_model_tmp=$pdl_model->slice(":,$J_MIN");
#    $pdl_model_tmp=$pdl_model_tmp*$NORM_J;



#    print "SSP CHI=$CHI_MIN $age_mod[$J_MIN] $met_mod[$J_MIN]\n";

    if ($CHI_MIN<$MIN_CHISQ) {    
	$MIN_CHISQ=$CHI_MIN;
	@Av_MIN=@Av;
	$age_min=$age_mod[$J_MIN];
	$met_min=$met_mod[$J_MIN];

    }
    return($CHI_MIN);
}



################ FIT_NEW #########################


sub FIT_Z_SIGMA {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $k,$j,$i,$iter;

    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }
    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];
    $rsigma=$sigma/$dpix_c;

# CONVOLUTION
    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;
    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
# END CONVOLUTION

    
    $pdl_dust=zeroes($i,$j);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    if ($Av[$j]!=0) {
		$dust=10**(-0.4*$Av[$j]*$dust_rat);  
		set($pdl_dust,$i,$j,$dust);
	    } else {
		set($pdl_dust,$i,$j,1);
	    }
	}
    }


    my $model_now=zeroes($n_unc);    
    for ($iii=0;$iii<$nf;$iii++) {    
	$pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	$C=$coeffs->at($iii,0);
	my $dust_spec=$out_spec_pdl*$pdl_dust->slice(",$iii");
	$model_now=$model_now+$C*$dust_spec;
    }

    $chi=0;
    $chi2=0;
    $NFREE=0;
    for ($j=0;$j<$n_unc;$j++) {
#	print "$j $wave_unc[$j] '$MIN_W,$MAX_W'\n";
	if (($wave_unc[$j]>$MIN_W)&&($wave_unc[$j]<$MAX_W)) {
	    $out_spec[$j]=$model_now->at($j,0);
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)) {
		$chi_j[$j]=$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
		$chi=$chi+$chi_j[$j];
		$NFREE++;
	    }
	}
    }

   $chi_sq=$chi;
    if ($NFREE>0) {
	$chi_sq=($chi_sq/($NFREE))**0.5;
    }

    return $chi_sq;
}
