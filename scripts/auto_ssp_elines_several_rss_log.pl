#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
# MOD: 1/3/2010

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


if ($#ARGV<5) {
    print "USE: auto_ssp.pl RSS.FITS BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift min_flux_fit_all min_flux_fit_continuum]\n";
    print "CONFIG_FILE:\n";
    print "redshift delta_redshift min_redshift max_redshift\n";
    print "sigma delta_sigma min_sigma max_sigma\n";
    print "Av delta_Av min_Av max_Av\n";
    print "N_SYSTEMS\n";
    print "(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "...\n";
    print "(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX\n";
    print "start_w_peak end_w_peak\n";
    print "wavelength_to_norm width_AA new_back_templates.fits\n";
close(FH);
    exit;
}

#$dev="/xterm";
$dev="/xs";
$sum_min=0.01;
$cube_file=$ARGV[0];
if ($ARGV[0] !~ "rss") {
    $call="cp ".$cube_file." rss_".$cube_file;
    system($call);
    $cube_file="rss_".$cube_file;
}



$clean_file=$cube_file;
$clean_file="clean_".$cube_file;
$org_file="org_".$cube_file;
$res_file="res_".$cube_file;
$res_ssp_file="res_ssp_".$cube_file;
$mod_file="mod_".$cube_file;




$back_list=$ARGV[1];
$outfile=$ARGV[2];
$out_elines="elines_".$outfile;
$out_single="single_".$outfile;
$call="rm ".$outfile;
system($call);
$call="rm ".$out_elines;
system($call);
$call="rm ".$out_single;
system($call);


$mask_list=$ARGV[3];
$config_file=$ARGV[4];
$plot=$ARGV[5];

	if ($plot==1) {
	    $dev="/xs";
	} else {
	    $dev="/null";
	}

$smooth=1;



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
    $def=2;
}

$sum_min=0.05;
$sum_min_fit=0.05;
if ($#ARGV==15) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sum_min=$ARGV[15];
    $def=2;
}

$sum_min_fit=$sum_min;
if ($#ARGV==16) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sum_min=$ARGV[15];
    $sum_min_fit=$ARGV[16];
    $def=2;
}

open(FH,"<$config_file");
$line=<FH>;
chop($line);
($redshift,$d_redshift,$min_redshift,$max_redshift)=split(" ",$line);
$line=<FH>;
chop($line);
if ($input_redshift>0) {
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
}
($sigma,$d_sigma,$min_sigma,$max_sigma)=split(" ",$line);
$line=<FH>;
chop($line);
($Av,$d_Av,$min_Av,$max_Av)=split(" ",$line);
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
$line=<FH>;
chop($line);
($start_w_peak,$end_w_peak)=split(" ",$line);
$line=<FH>;
chop($line);
($wave_norm,$w_wave_norm,$new_back_file)=split(" ",$line);
close(FH);

$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};



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
#print "$fit_new\n";
#<stdin>;
#print "F=$start_w_peak,$end_w_peak\n"; exit;
#$call="cp ".$config_e." tmp_e.config";
#system($call);



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
$pdl_cube=rfits($cube_file);
$pdl_clean=$pdl_cube;
($nz,$nx)=$pdl_cube->dims;
$ny=1;
$crpix3=$pdl_cube->hdr->{CRPIX1};
$cdelt3=$pdl_cube->hdr->{CDELT1};
if ($cdelt3==0) {
    $cdelt3=$pdl_cube->hdr->{CD1_1};
}
$crval3=$pdl_cube->hdr->{CRVAL1};
$h=$pdl_cube->gethdr;

$mod_cube=zeroes($nz,$nx);
$res_cube=zeroes($nz,$nx);
$res_ssp_cube=zeroes($nz,$nx);
$org_cube=zeroes($nz,$nx);

$chi_map=zeroes($nx,$ny);
$flux_map=zeroes($nx,$ny);
$age_map=zeroes($nx,$ny);
$met_map=zeroes($nx,$ny);
$dust_map=zeroes($nx,$ny);
$vel_abs_map=zeroes($nx,$ny);
$vel_em_map=zeroes($nx,$ny);
$disp_abs_map=zeroes($nx,$ny);
$disp_em_map=zeroes($nx,$ny);



$age_file=$cube_file;
$age_file =~ s/rss/age/;
$dust_file=$cube_file;
$dust_file =~ s/rss/dust/;
$met_file=$cube_file;
$met_file =~ s/rss/met/;
$chi_file=$cube_file;
$chi_file =~ s/rss/chi/;
$flux_file=$cube_file;
$flux_file =~ s/rss/flux/;
$vel_abs_file=$cube_file;
$vel_abs_file =~ s/rss/vel_abs/;
$disp_abs_file=$cube_file;
$disp_abs_file =~ s/rss/disp_abs/;
$vel_em_file=$cube_file;
$vel_em_file =~ s/rss/vel_em/;
$disp_em_file=$cube_file;
$disp_em_file =~ s/rss/disp_em/;


($redshift_ini,$d_redshift_ini,$min_redshift_ini,$max_redshift_ini)=($redshift,$d_redshift,$min_redshift,$max_redshift);
$sigma_ini=$sigma;
$d_sigma_ini=$d_sigma;
$Av_ini=$Av;
$d_Av_ini=$d_Av;

for ($ix=0;$ix<$nx;$ix++) {
    for ($iy=0;$iy<$ny;$iy++) {

	$MIN_CHISQ=1e12;
	($redshift,$d_redshift,$min_redshift,$max_redshift)=($redshift_ini,$d_redshift_ini,$min_redshift_ini,$max_redshift_ini);
	$sigma=$sigma_ini;
	$d_sigma=$d_sigma_ini;
	$Av=$Av_ini;
	$d_Av=$d_Av_ini;
	
	$pdl_slice=$pdl_cube->slice(":,($ix)");
#	print "$pdl_slice\n";
	$iz0=int(0.25*$nz);
	$iz1=int(0.75*$nz);
	$pdl_slice2=$pdl_slice->slice("$iz0:$iz1");
	$sum=sum($pdl_slice2)/($iz1-$iz0);
	print "$ix/$nx, $iy/$ny $sum/$sum_min\n"; #<stdin>;
	if ($sum>$sum_min) {



	$n_unc=0;
	$y_min=1e12;
	$y_max=-1e12;
	$i_scale=0;
	$FLUX=0;
	for ($n_unc=0;$n_unc<$nz;$n_unc++) {
	    $index_unc[$n_unc]=$n_unc+1;
	    $wave_unc[$n_unc]=$crval3+$cdelt3*(($crpix3-1)+$n_unc);
	    $flux_unc[$n_unc]=$pdl_slice->at($n_unc);

	    if ($flux_unc[$n_unc] eq "nan") {
		$flux_unc[$n_unc]=$flux_unc[$n_unc-1];
		$masked[$n_unc]=0;
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

#	    if ($flux_unc[$n_unc]==0) {
#		$masked[$n_unc]=0;
#	    }


	    $flux_masked[$n_unc]=$flux_unc[$n_unc]*$masked[$n_unc];
	    if (($masked[$n_unc]==1)&&($flux_unc[$n_unc]!=0)) {
		$e_flux_unc[$n_unc]=1/abs($flux_unc[$n_unc]);
	    } else {
		$e_flux_unc[$n_unc]=0;
	    }
	    if (($wave_unc[$n_unc]>$min_wave)&&($wave_unc[$n_unc]<$max_wave)) {
		$FLUX=$FLUX+$flux_masked[$n_unc];
	    }

	}


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






#
# We create a kernel
#


#print "We start the fit...\n";

#open(OUTFILE,">$outfile");

	$plot_old=$plot;
$min_chi_sq=FIT($redshift,$sigma,$Av);



################
# We look for the redshift;

	    if ($sum>$sum_min_fit) {
		if ($d_redshift!=0) {
		    $RED=$min_redshift;
		    while ($RED<$max_redshift) {
			$RED=$RED+$d_redshift;
			$chi_now=FIT($RED,$sigma,$Av);
			if ($chi_now<$min_chi_sq) {
			    $redshift=$RED;	
			    $min_chi_sq=$chi_now;
			}
#			print "PASO = $RED $redshift $chi_now $min_chi_sq\n";

		    }
		    
		    $fit_redshift=0;
		    $d_redshift=0.0;
		}
		    $VEL_red=$redshift*$vel_light;
	    } else {
		$fit_redshift=0;
		$d_redshift=0.00;
		$d_sigma=0;
		$d_Av=0;
	    }

 
		my @new_flux_unc;
		for ($i=0;$i<$n_unc;$i++) {
		    $res_flux_unc[$i]=$flux_unc[$i]-$model_spec_min[$i];
		    $new_flux_unc[$i]=$flux_unc[$i];
		    $res_masked[$i]=$res_flux_unc[$i]*$masked[$i];
		}

	$median_FLUX=median(@flux_unc);
	$sigma_FLUX=sigma(@res_masked);

	$max=4*$median_FLUX;

#@median=median_filter(51,\@res_flux_unc);
		
		open(CLEAN,">$clean_file");
		for ($i=0;$i<$n_unc;$i++) {
		    $flux_unc[$i]=$flux_unc[$i]-$median[$i];
		    set($pdl_clean,$i,$ix,$flux_unc[$i]);
		    print CLEAN "$i $wave_unc[$i] $flux_unc[$i] $new_flux_unc[$i] $median[$i]\n";
		}
		close(CLEAN);





#exit;

#exit;

#print "($redshift,$sigma,$Av) $min_chi_sq\n";
#print "----------------------\n";
$delta_chi=10;
#while (($min_chi_sq>1)&&($delta_chi>0.1)) {
$NITER=0;
	$niter_tmp_max=50;

	$MIN_CHISQ=1e12;
while (($MIN_CHISQ>$MIN_DELTA_CHISQ)&&($NITER<$MAX_NITER)) {
    $plot=1;
    if ($NITER==1) {
	$MIN_CHISQ=1e12;
    }
    if (($d_redshift!=0)&&($fit_redshift==1)) {
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift+$d_redshift,$sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_redshift=-$d_redshift;
	}
	$niter_tmp=0;
	do {
	    $chi_now[0]=$min_chi_sq;
	    $redshift=$redshift+$d_redshift;
	    for ($iter=1;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift+($iter)*$d_redshift,$sigma,$Av);
		$RED=$redshift+$iter*$d_redshift;
	    }
	    $min_chi_sq=$chi_now[1];
	    $niter_tmp++;
	} while (($chi_now[2]<$chi_now[1])&&($niter_tmp<$niter_tmp_max));
	$redshift=($redshift+2*$d_redshift)-$d_redshift*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);

	if ($redshift>$max_redshift) {
	    $redshift=$max_redshift;
	    $d_redshift=-0.5*abs($d_redshift);
	}
	if ($redshift<$min_redshift) {
	    $redshift=$min_redshift;
	    $d_redshift=-0.5*abs($d_redshift);
	}
    }
    $fit_redshift=0;
	$red_abs=$redshift;
print "REDSHIFT = $redshift, VEL = $VEL_red\n";




    if ($d_sigma!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma+$d_sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_sigma=-$d_sigma;
	}
#	print "SIGMA ($redshift,$sigma,$Av) $min_chi_sq\n";
	$niter_tmp=0;
	do {
	    $chi_now[0]=$min_chi_sq;
	    $sigma=$sigma+$d_sigma;

	    for ($iter=0;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma+($iter)*$d_sigma,$Av);
	    }
	    $min_chi_sq=$chi_now[1];
	    $niter_tmp++;
	} while (($chi_now[2]<$chi_now[1])&&($niter_tmp<$niter_tmp_max));
#	} while ($chi_now[2]<$chi_now[1]);
	
	$sigma=($sigma+2*$d_sigma)-$d_sigma*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);

    if ($sigma>$max_sigma) {
	$sigma=$max_sigma;
#0.5*($max_sigma+$min_sigma);
	$d_sigma=-0.5*abs($d_sigma);
    }
    if ($sigma<$min_sigma) {
	$sigma=$min_sigma;
#0.5*($max_sigma+$min_sigma);
	$d_sigma=0.5*abs($d_sigma);
    }	
    }

    if ($d_Av!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma,$Av+$d_Av);
 	if ($chi_test>$min_chi_sq) {
	    $d_Av=-$d_Av;
	}
#	print "AV ($redshift,$sigma,$Av) $min_chi_sq\n";
	$niter_tmp=0;
	do {
	    $Av=$Av+$d_Av;
	    $chi_now[0]=$min_chi_sq;

	    for ($iter=1;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma,$Av+($iter)*$d_Av);
	    }
	    $min_chi_sq=$chi_now[1];
	    $niter_tmp++;
	} while (($chi_now[2]<$chi_now[1])&&($niter_tmp<$niter_tmp_max));	    
#	} while ($chi_now[2]<$chi_now[1]);
	$Av=($Av+2*$d_Av)-$d_Av*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);
	
 	if ($Av>$max_Av) {
	    $Av=$max_Av;
	    $d_Av=0;#-0.5*abs($d_Av);
	}
	if ($Av<$min_Av) {
	    $Av=$min_Av;
		$d_Av=0;#0.5*abs($d_Av);
	}
    }
    $min_chi_sq_end=FIT($redshift,$sigma,$Av);
    $delta_chi=abs(($min_chi_sq_end-$min_chi_sq)/$min_chi_sq);
#    print "END ($redshift,$sigma,$Av) $min_chi_sq $delta_chi\n";

    $d_redshift=0;#0.05*$d_redshift;
    $d_sigma=$d_sigma/2;
    $d_Av=$d_Av/2;
    $NITER++;
    if ($NITER>1) {
	$d_Av=0;
    }
    if ($NITER>2) {
	$d_sigma=0;
    }



#($start_w_e,$end_w_e,$mask_e,$config_e,$npoly_e,$mask_poly_e)=split(" ",$line);
    # We fit now the emission lines!
    #

    $wpeak=6562;
    $Fpeak=-1e12;


    $ks=0;
    $SYS_VEL=$vel_light*$redshift;

    for ($is=0;$is<$ns;$is++) {
	if ($is==0) {
	    $SYS_VEL_MAX=$vel_light*$redshift*1.2;
	    $SYS_VEL_MIN=$vel_light*$redshift*0.8;
	} else {
	    $SYS_VEL_MAX=$vel_light*$redshift*1.05;
	    $SYS_VEL_MIN=$vel_light*$redshift*0.95;
	}
	$start_w_e=$start_w_E[$is];
	$end_w_e=$end_w_E[$is];
	$mask_e=$mask_E[$is];
	$config_e=$config_E[$is];
	$npoly_e=$npoly_E[$is];
	$mask_poly_e=$mask_poly_E[$is];
	$nmin_e=$nmin_E[$is];
	$nmax_e=$nmax_E[$is];

	
	print "CONF=$config_e\n";
	open(IN_CONF,"<$config_e");
	open(OUT_CONF,">tmp_e.config");
	$n_conf=0;
	while($line=<IN_CONF>) {
	    chop($line);
	    if ($n_conf==2) {
		@data=split(" ",$line);
		$wpeak_res=$data[0];
	    }
	    if (($n_conf==5)) {
		print OUT_CONF "$SYS_VEL   1   $SYS_VEL_MIN  $SYS_VEL_MAX  -1\n";	    
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
#	$call="spec2img.pl sub_spec.txt res_spec.fits 1";
#	system($call);
#	$call="polyfit_mask.pl res_spec.fits 6 res_poly.fits ".$mask_poly_e." ".$nmin_e." ".$nmax_e;
 #	print "$call\n";
#	exit;
#	system($call);
#	$call="imarith.pl res_spec.fits - res_poly.fits new_res.fits";
#	system($call);
#	$call="img2spec.pl new_res.fits 0 new_res.txt";
#	system($call);
#	$call="fit_spec_back.pl new_res.txt none ".$start_w_e." ".$end_w_e." none 0 ".$mask_e." tmp_e.config";


       $call="fit_spec_back.pl sub_spec.txt none ".$start_w_e." ".$end_w_e." none 0 ".$mask_e." tmp_e.config ".$dev." > junk.junk";
	system($call);
	print "$call\n";
	print "DONE $is\n";
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
		if ($is==0) {
		    @dat_now=split(" ",$lineKS);
		    if ($dat_now[3]>$sum_min) {
			if ($fit_redshift==1) {
			    $redshift=$dat_now[7]/$vel_light;
			}
			$SYS_VEL=$dat_now[7];
		    }
		}
		$ks++;
	    }
	}
	close(KS);
	print "REDSHIFT_ELINES = $redshift , VEL = $SYS_VEL\n";
	$red_em=$redshift;
#	<stdin>;
    }

    if ($d_redshift_ini!=0) {
	$d_redshift=0.00001;
	$fit_redshift=1;
    }
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
#	if ($is==0) {
#	    set($vel_em_map,$ix,$iy,$info[7]);
#	    set($disp_em_map,$ix,$iy,($info[5]/$info[1])*$info[7]);
#	}
    }
    close(OUT_CONF);    


    $call="fit_spec_back.pl sub_spec.txt none ".$start_w_min." ".$end_w_max." none 0 none tmp_e.config ".$dev." > junk.junk";
    system($call);
    print "$call\n";
    #
    # Now we subtract the emission lines
    # 

    open(MOD,"<out_mod_res.fit_spectra");
    $NF=0;
    @model_joint=@model_spec_min;
    @res_joint=@res_spec;
    while($line=<MOD>) {
	chop($line);
	@data=split(" ",$line);
	$WMOD=$data[0];
	$FMOD=$data[2];
	while (abs($data[0]-$wave_unc[$NF])>0.1) {
	    $NF++;
	} 
	$flux_now=$flux_unc[$NF]-$FMOD;	    
	$flux_masked[$NF]=$flux_now*$masked[$NF];
	$model_joint[$NF]=$model_spec_min[$NF]+$FMOD;
	$res_joint[$NF]=$res_spec[$NF]-$FMOD;
    }
    close(MOD);


    open(OUT_MOD,">mod_joint_spec.txt");
    open(OUT_RES,">res_joint_spec.txt");
    for ($j=0;$j<$n_unc;$j++) {
	print OUT_MOD "$index_unc[$j] $wave_unc[$j] $model_joint[$j]\n";
	print OUT_RES "$index_unc[$j] $wave_unc[$j] $res_joint[$j]\n";
    }
    close(OUT_RES);
    close(OUT_MOD);

#print "PASO\n";	<stdin>;


    for ($j=0;$j<$n_unc;$j++) {
	set($mod_cube,$j,$ix,$model_joint[$j]);
	set($res_cube,$j,$ix,$res_joint[$j]);
	set($res_ssp_cube,$j,$ix,$res_spec[$j]);
	set($org_cube,$j,$ix,$flux_unc[$j]);
    }	


#    print "Press Enter"; <sdtin>;
    #
    #



    #    print "END****\n";
}
close(OUT);


#print "----------------------------\n";
print "CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av REDSHIFT=$redshift SIGMA_DISP=$sigma  AGE_mass=$age_min_mass MET_mass=$met_min_mass\n";

#	<stdin>;
open(FH,">compare_back_list_dust_mult.out");
print FH "$nx $ny $ix $iy $MIN_CHISQ $age_min $met_min $Av $redshift $sigma $FLUX $red_abs $red_em $median_FLUX $sigma_FLUX  $age_min_mass $met_min_mass\n";
close(FH);

       $call="cat compare_back_list_dust_mult.out >> ".$outfile;
system($call);

	if ($wave_norm>0) {
	print "*************** FIT_SINGLE SSP ******************\n";
	$MIN_CHISQ=1e12;
	$plot=$plot_old;
#	print "Press Enter"; <stdin>;



	$min_chi_sq=FIT_SINGLE($redshift,$sigma,$Av);
	$rms_single=sigma(@res_clean);
	$med_flux_single=median(@flux_clean);


	print "CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av REDSHIFT=$redshift SIGMA_DISP=$sigma\n";
#	print "Press Enter"; <stdin>;
	open(FH,">compare_back_list_dust_mult.out");
print FH "$nx $ny $ix $iy $MIN_CHISQ $age_min $met_min $Av $redshift $sigma $FLUX $red_abs $red_em $med_flux_single $rms_single $age_min $met_min\n";
close(FH);

       $call="cat compare_back_list_dust_mult.out >> ".$out_single;
	system($call);
	}


	open(OELINES,">>$out_elines");
	print OELINES "$nx $ny $ix $iy\n";
	close(OELINES);
	$call="cat out.fit_spectra >> ".$out_elines;
	system($call);

#	print "Press Enter "; <stdin>;

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


$res_cube->wfits($res_file);
$org_cube->wfits($org_file);
$mod_cube->wfits($mod_file);
#	print "PRESS "; <stdin>;
	}

	
    }
}

$pdl_clean->wfits($clean_file);


exit;


sub FIT_SINGLE_OLD {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av=$_[2];

#
# We renormalize the models
#


    for ($j=0;$j<$n_c_new;$j++) {
	$wave_c_new[$j]=($crval_new+$cdelt_new*($j+1-$crpix_new))*(1+$redshift);
    }
    $dpix_c=$wave_c_new[1]-$wave_c_new[0];
    $rsigma=$sigma/$dpix_c;


    my $CHI_MIN=1e12;

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
#    $header="NORM".$iii;
#    $ml[$iii]=$pdl_flux_c_ini->hdr->{$header};
    
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

    $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
    my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c_new), $pdl_flux_c);
    ($n_c_out)=$out_spec_pdl->dims;

    

    for ($i=0;$i<$n_unc;$i++) {
	if ($masked[$i]>0) {
	    $error[$i]=$e_flux_unc[$i];
	} else {
	    $error[$i]=0;
	}
	$val=$out_spec_pdl->at($i);
	$model[$i][$iii]=$val;#*$masked[$i];
	if ($NITER==0) {
	    for ($j=0;$j<$nline;$j++) {
		if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
		    $error[$i]=0.1*$e_flux_unc[$i];
		}
	    }
	}
    }
    
}

    $pdl_model=zeroes($n_unc,$nf_new);
    $pdl_error=zeroes($n_unc,$nf_new);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf_new;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av*$dust_rat);  
	    $val=$model[$i][$j]*$dust;
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
 

    my $nx1=int(($wave_norm-$w_wave_norm-$crval3)/$cdelt3);
    my $nx2=int(($wave_norm+$w_wave_norm-$crval3)/$cdelt3);
#    print "$crval_new $cdelt_new $wave_norm $w_wave_norm $nx1 $nx2\n";
    my $pdl_sec_cen=$pdl_flux_masked->slice("$nx1:$nx2");
    my $SUM=sum($pdl_sec_cen);
    
    my $pdl_chi=zeroes($nf_new);
    my $pdl_e=pdl(@error);
#    print "$pdl_e\n";
    $J_MIN=-1;
   for ($j=0;$j<$nf_new;$j++) {
       my $pdl_model_tmp=$pdl_model->slice(":,$j");
       my $pdl_sec_cen=$pdl_model_tmp->slice("$nx1:$nx2");
       my $SUM_NOW=sum($pdl_sec_cen);
       my $NORM=$SUM/$SUM_NOW;
       $pdl_model_tmp.=$pdl_model_tmp*$NORM;
       my $pdl_chi_now=$pdl_masked*(($pdl_flux_masked-$pdl_model_tmp)**2)*$pdl_e;
;#/sqrt($pdl_flux_masked*$pdl_flux_masked);
       my $chi_sq=sum($pdl_chi_now);
       $chi_sq=($chi_sq/($n_unc-1))**0.5;
       set($pdl_chi,$j,$chi_sq);
#       print "$chi_sq $CHI_MIN $age_mod[$j] $met_mod[$j] $j $J_MIN\n";
#    print "SSP $chi_sq CHI=$CHI_MIN $age_mod[$J_MIN] $met_mod[$J_MIN]\n";
       if ($CHI_MIN>$chi_sq) {
	   $CHI_MIN=$chi_sq;
	   $J_MIN=$j;
       }

   }
#    <stdin>;
    $pdl_model_tmp=$pdl_model->slice(":,$J_MIN");
    @out_spec=list($pdl_model_tmp);


#    print "SSP CHI=$CHI_MIN $age_mod[$J_MIN] $met_mod[$J_MIN]\n";

    if ($CHI_MIN<$MIN_CHISQ) {    
	$MIN_CHISQ=$CHI_MIN;
	$Av_MIN=$Av;
	for ($j=0;$j<$n_unc;$j++) {
	    $model_spec[$j]=0;
	    $wave_res=$wave_unc[$j]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av*$dust_rat);  
	    $model_spec[$j]=$out_spec[$j];
	    $model_spec_min[$j]=$model_spec[$j];
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
#	    print "$j $flux_unc[$j] $model_spec[$j] $res_spec[$j]\n";
	}
	$age_min=$age_mod[$J_MIN];
	$met_min=$met_mod[$J_MIN];
	
    }

#    <stdin>;
    if ($plot==1) {
	pgbegin(0,$dev,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Wavelength","Counts","SSP X=$CHI_MIN T=$age_min Z=$met_min Av=$Av z=$redshift sigma=$sigma");
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

#    print "Press Enter"; <stdin>;
    return($CHI_MIN);
}

sub FIT {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av=$_[2];

#    print "Av= $Av $min_Av-$max_Av\n";
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

    $header="NORM".$iii;
#    $ml[$iii]=$pdl_flux_c_ini->hdr->{$header};

 $val_ml=$pdl_flux_c_ini->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }



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
	    $error[$i][$iii]=1/abs($val+0.001)+$chi_sec[$i];
	  
#	    $error[$i][$iii]=10;
	} else {
	    $error[$i][$iii]=1e20;
	}

	if ($NITER==0) {
	    for ($j=0;$j<$nline;$j++) {
		if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
		    $error[$i][$iii]=1e10;
#		print "$wave_unc[$i] $error[$i][$iii]\n";
		}
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
	$chi_sec[$j]=0;

	    if ($flux_unc[$j]!=0) {
		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($flux_unc[$j]);
		$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($flux_unc[$j]);

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
	    $norm_C_mass=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		$norm_C=$norm_C+$C;
		$norm_C_mass=$norm_C_mass+$C*$ml[$k];
		$model_spec[$j]=$model_spec[$j]+$C*$model_no_mask[$j][$k]*$dust;
		$model_spec_min[$j]=$model_spec[$j];
	    }
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    

	}
	$age_min=0;
	$met_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);

	    $age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
#	    $age_min=$age_min+$C*$age_mod[$k]/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
#	    print "COEFF = $k $ml[$k] $norm_C_mass $age_min_mass $met_min_mass\n";
	    if ($norm_C_mass>0) {
	#	$age_min_mass=$age_min_mass+$C*$ml[$k]*$age_mod[$k]/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
	    }

	}
	$age_min=10**($age_min);
#	$met_min=10**($met_min);
	$age_min_mass=10**($age_min_mass);

#	$age_min=apr_n($age_min,4);
#	$met_min=apr_n($met_min,6);
#	$age_min_mass=apr_n($age_min_mass,4);
#	$met_min_mass=apr_n($met_min_mass,6);

	
    }
    $name=$cube_file.", (".$ix.",".$iy.") (".$nx.",".$ny.") ";
    $scale="1";
#    $chi_sq=($chi)**0.5;

#    print "CHI_SQ=$chi_sq\n";
#    print "$name $chi_sq $scale\n";
#    print OUT "$name $chi_sq $scale $Av \n";


    
    if ($plot==1) {
	pgbegin(0,$dev,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Wavelength","Counts","MSP X=$chi_sq T=$age_min($age_min_mass) Z=$met_min($met_min_mass) Av=$Av z=$redshift sigma=$sigma");
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





#    print "Press Enter"; <stdin>;
    return $chi_sq;
}


sub FIT_SINGLE {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av=$_[2];

#
# We renormalize the models
#


    for ($j=0;$j<$n_c_new;$j++) {
	$wave_c_new[$j]=($crval_new+$cdelt_new*($j+1-$crpix_new))*(1+$redshift);
    }
    $dpix_c=$wave_c_new[1]-$wave_c_new[0];
    $rsigma=$sigma/$dpix_c;


    my $CHI_MIN=1e12;

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
#	print "$iii $age_mod[$iii] $met_mod[$iii]\n";
	
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

	$pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c_new), $pdl_flux_c);
	($n_c_out)=$out_spec_pdl->dims;
	
	
	
	for ($i=0;$i<$n_unc;$i++) {
	    if ($masked[$i]>0) {
		$error[$i]=$e_flux_unc[$i];
	    } else {
		$error[$i]=0;
	    }
	    $val=$out_spec_pdl->at($i);
	    $model[$i][$iii]=$val;#*$masked[$i];
	    if ($NITER==0) {
		for ($j=0;$j<$nline;$j++) {
		    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
			$error[$i]=0.1*$e_flux_unc[$i];
		    }
		}
	    }
	}
	
    }

    $pdl_model=zeroes($n_unc,$nf_new);
    $pdl_error=zeroes($n_unc,$nf_new);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf_new;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    if ($j==0) {
		$wave_res=$wave_unc[$i]/(1+$redshift);
		$dust_rat=A_l(3.1,$wave_res);
		$dust=10**(-0.4*$Av*$dust_rat);  
		$a_dust[$i]=$dust;
	    }
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
#       my $pdl_chi_now=$pdl_masked*(($pdl_flux_masked-$pdl_model_tmp)**2)*$pdl_e;
       my $pdl_chi_now=$pdl_masked*(($pdl_flux_masked-$yfit)**2)*$pdl_e;
       my $chi_sq=sum($pdl_chi_now);
       $chi_sq=($chi_sq/($n_unc-1))**0.5;
       set($pdl_chi,$j,$chi_sq);


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
	    pgbegin(0,$dev,1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	    pgsch(1.2);           # Set character height
	    pglabel("Wavelength","Counts","SSP X=$CHI_MIN T=$age_mod[$J_MIN] Z=$met_mod[$J_MIN] Av=$Av z=$redshift sigma=$sigma");
	    pgsci(1);
	    pgline($n_unc,\@wave_unc,\@flux_unc);    
	    pgsci(2);
	    pgline($n_unc,\@wave_unc,\@out_spec);    
#	    pgsci(8);
#	    pgline($n_unc,\@wave_unc,\@model_spec_min);    
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
	$Av_MIN=$Av;
	$age_min=$age_mod[$J_MIN];
	$met_min=$met_mod[$J_MIN];

    }
    return($CHI_MIN);
}



