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
    print "USE: auto_ssp.pl CUBE.FITS BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift min_flux_fit_all min_flux_fit_continuum]\n";
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
close(FH);
    exit;
}

#$dev="/xterm";
$RAND=int(rand(30));
$dev=$RAND."/xs";
$sum_min=0.05;
$cube_file=$ARGV[0];


$back_list=$ARGV[1];
$outfile=$ARGV[2];
$out_elines="elines_".$outfile;
$call="rm ".$outfile;
system($call);
$call="rm ".$out_elines;
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
close(FH);

$line=<FH>;
chop($line);
($FIX_vel,$in_vel_map,$FIX_disp,$in_disp_map,$FIX_dust,$in_dust_map)=split(" ",$line);
close(FH);

if ($FIX_vel==1) {
    $pdl_in_redshift=rfits($in_vel_map);
    $pdl_in_redshift=$pdl_in_redshift/$vel_light;
    $input_d_redshift=0;
    $d_redshift=0;
}
if ($FIX_disp==1) {
    $pdl_in_sigma=rfits($in_disp_map);
    $d_sigma=0;
}
if ($FIX_dust==1) {
    $pdl_in_Av=rfits($in_dust_map);
    $d_Av=0;
}





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


$pdl_flux_c_ini=rfits($back_list);
($n_c,$nf)=$pdl_flux_c_ini->dims;
$crpix=$pdl_flux_c_ini->hdr->{CRPIX1};
$cdelt=$pdl_flux_c_ini->hdr->{CDELT1};
$crval=$pdl_flux_c_ini->hdr->{CRVAL1};

#
# We read the input spectrum
#
$pdl_cube=rfits($cube_file);
($nx,$ny,$nz)=$pdl_cube->dims;
$crpix3=$pdl_cube->hdr->{CRPIX3};
$cdelt3=$pdl_cube->hdr->{CDELT3};
if ($cdelt3==0) {
    $cdelt3=$pdl_cube->hdr->{CD3_3};
}
$crval3=$pdl_cube->hdr->{CRVAL3};
$h=$pdl_cube->gethdr;

$mod_cube_file="mod_".$outfile.".fits";
$res_cube_file="res_".$outfile.".fits";

$mod_cube=zeroes($nx,$ny,$nz);
$res_cube=zeroes($nx,$ny,$nz);

$mod_cube->sethdr( $h );
$res_cube->sethdr( $h );


$chi_map=zeroes($nx,$ny);
$flux_map=zeroes($nx,$ny);
$age_map=zeroes($nx,$ny);
$met_map=zeroes($nx,$ny);
$age_mass_map=zeroes($nx,$ny);
$met_mass_map=zeroes($nx,$ny);
$dust_map=zeroes($nx,$ny);
$vel_abs_map=zeroes($nx,$ny);
$vel_em_map=zeroes($nx,$ny);
$disp_abs_map=zeroes($nx,$ny);
$disp_em_map=zeroes($nx,$ny);

#$NAMEFILE=$cube_file;
@NAMES=split(/\//,$cube_file);
$NAMEFILE=$NAMES[$#NAMES];
#print "NAME= $NAMEFILE\n";



$age_file=$NAMEFILE;
$age_file =~ s/cube/age/;
$age_mass_file=$NAMEFILE;
$age_mass_file =~ s/cube/age_mass/;
$dust_file=$NAMEFILE;
$dust_file =~ s/cube/dust/;
$met_file=$NAMEFILE;
$met_file =~ s/cube/met/;
$met_mass_file=$NAMEFILE;
$met_mass_file =~ s/cube/met_mass/;
$chi_file=$NAMEFILE;
$chi_file =~ s/cube/chi/;
$flux_file=$NAMEFILE;
$flux_file =~ s/cube/flux/;
$vel_abs_file=$NAMEFILE;
$vel_abs_file =~ s/cube/vel_abs/;
$disp_abs_file=$NAMEFILE;
$disp_abs_file =~ s/cube/disp_abs/;
$vel_em_file=$NAMEFILE;
$vel_em_file =~ s/cube/vel_em/;
$disp_em_file=$NAMEFILE;
$disp_em_file =~ s/cube/disp_em/;


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
	
	$pdl_slice=$pdl_cube->slice("($ix),($iy),:");
	$iz0=int(0.25*$nz);
	$iz1=int(0.75*$nz);
	$pdl_slice2=$pdl_slice->slice("$iz0:$iz1");
	$sum=sum($pdl_slice2)/($iz1-$iz0);
#	print "$ix/$nx, $iy/$ny $sum/$sum_min\n";
	if ($sum>$sum_min) {
	print "$ix/$nx, $iy/$ny $sum/$sum_min\n";


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
	if ($FIX_vel==1) {
	    $redshift=$pdl_in_redshift->at($ix,$iy);
	}
	if ($FIX_disp==1) {
	    $sigma=$pdl_in_sigma->at($ix,$iy);
	}
	if ($FIX_dust==1) {
	    $Av=$pdl_in_Av->at($ix,$iy);
	}
	

	$min_chi_sq=FIT($redshift,$sigma,$Av);

################
# We look for the redshift;

	if ($FIX_vel==0) {
	    if ($sum>$sum_min_fit) {
		
		$RED=$min_redshift;
		while ($RED<$max_redshift) {
		    $RED=$RED+$d_redshift;
		    $chi_now=FIT($RED,$sigma,$Av);
		    if ($chi_now<$min_chi_sq) {
			$redshift=$RED;	
			$min_chi_sq=$chi_now;
		    }
		}

		$fit_redshift=1;
		$d_redshift=0.00005;
		$VEL_red=$redshift*$vel_light;
	    } else {
		$fit_redshift=0;
		$d_redshift=0.01;
		$d_sigma=0;
		$d_Av=0;
	    }
	}


		my @new_flux_unc;
		for ($i=0;$i<$n_unc;$i++) {
		    $res_flux_unc[$i]=$flux_unc[$i]-$model_spec_min[$i];
		    $new_flux_unc[$i]=$flux_unc[$i];
		}
#@median=median_filter(51,\@res_flux_unc);
		





#exit;

#exit;

#print "($redshift,$sigma,$Av) $min_chi_sq\n";
#print "----------------------\n";
$delta_chi=10;
#while (($min_chi_sq>1)&&($delta_chi>0.1)) {
$NITER=0;
#	print "***** $MIN_CHISQ>$MIN_DELTA_CHISQ\n";
while (($MIN_CHISQ>$MIN_DELTA_CHISQ)&&($NITER<$MAX_NITER)) {
	if ($FIX_vel==0) {

    if (($d_redshift!=0)&&($fit_redshift==1)) {
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift+$d_redshift,$sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_redshift=-$d_redshift;
	}
	do {
	    $chi_now[0]=$min_chi_sq;
	    $redshift=$redshift+$d_redshift;
	    for ($iter=1;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift+($iter)*$d_redshift,$sigma,$Av);
		$RED=$redshift+$iter*$d_redshift;
	    }
	    $min_chi_sq=$chi_now[1];
	} while ($chi_now[2]<($chi_now[1]-$MIN_DELTA_CHISQ));
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
#print "REDSHIFT = $redshift, VEL = $VEL_red\n";




    if ($d_Av!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma,$Av+$d_Av);
	if ($chi_test>$min_chi_sq) {
	    $d_Av=-$d_Av;
	}
#	print "AV ($redshift,$sigma,$Av) $min_chi_sq\n";
	
	do {
	    $Av=$Av+$d_Av;
	    $chi_now[0]=$min_chi_sq;

	    for ($iter=1;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma,$Av+($iter)*$d_Av);
	    }
	    $min_chi_sq=$chi_now[1];
	    
	} while ($chi_now[2]<($chi_now[1]-$MIN_DELTA_CHISQ));
	$Av=($Av+2*$d_Av)-$d_Av*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);
	
	if ($Av>$max_Av) {
	    $Av=$max_Av;
		$d_Av=-0.5*abs($d_Av);
	}
	if ($Av<$min_Av) {
	    $Av=$min_Av;
		$d_Av=-0.5*abs($d_Av);
	}
    }


    if ($d_sigma!=0) {
	
	$min_chi_sq=FIT($redshift,$sigma,$Av);
	$chi_test=FIT($redshift,$sigma+$d_sigma,$Av);
	if ($chi_test>$min_chi_sq) {
	    $d_sigma=-$d_sigma;
	}
#	print "SIGMA ($redshift,$sigma,$Av) $min_chi_sq\n";
	
	do {
	    $chi_now[0]=$min_chi_sq;
	    $sigma=$sigma+$d_sigma;

	    for ($iter=0;$iter<3;$iter++) {
		$chi_now[$iter]=FIT($redshift,$sigma+($iter)*$d_sigma,$Av);
	    }
	    $min_chi_sq=$chi_now[1];
	} while ($chi_now[2]<($chi_now[1]-$MIN_DELTA_CHISQ));
	
	$sigma=($sigma+2*$d_sigma)-$d_sigma*(($chi_now[2]-$chi_now[1])/($chi_now[2]-2*$chi_now[1]+$chi_now[0])+0.5);

    if ($sigma>$max_sigma) {
	$sigma=$max_sigma;
	$d_sigma=-0.5*abs($d_sigma);
    }
    if ($sigma<$min_sigma) {
	$sigma=$min_sigma;
	$d_sigma=-0.5*abs($d_sigma);
    }	
    }




    $min_chi_sq_end=FIT($redshift,$sigma,$Av);
    $delta_chi=abs(($min_chi_sq_end-$min_chi_sq)/$min_chi_sq);
#    print "END ($redshift,$sigma,$Av) $min_chi_sq $delta_chi\n";

    $d_redshift=0;#0.05*$d_redshift;
    $d_sigma=$d_sigma/2;
    $d_Av=$d_Av/2;
    $NITER++;
    if ($NITER>2) {
	$d_Av=0;
    }
    if ($NITER>3) {
	$d_sigma=0;
    }
	} else {
	    $min_chi_sq=FIT($redshift,$sigma,$Av);
	    $NITER++;
	}


#($start_w_e,$end_w_e,$mask_e,$config_e,$npoly_e,$mask_poly_e)=split(" ",$line);
    # We fit now the emission lines!
    #

    $wpeak=6562;
    $Fpeak=-1e12;
#    open(SPEC,"<res_spec.txt");
#    while($line=<SPEC>) {
#	chop($line);
#	@data=split(" ",$line);
#	if (($data[1]>$start_w_peak)&&($data[1]<$end_w_peak)) {
#	    if ($data[2]>$Fpeak) {
#		$wpeak=$data[1];
#		$Fpeak=$data[2];
#		print "FPEAK = $wpeak $Fpeak\n";
#	    }
#	}
#    }
#    close(SPEC);


    if ($ns>0) {

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

	
#	print "CONF=$config_e\n";
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
	

       $call="fit_spec_back.pl sub_spec.txt none ".$start_w_e." ".$end_w_e." none 0 ".$mask_e." tmp_e.config ".$dev." > junk.junk";
	system($call);
#	print "$call\n";
#	print "DONE $is\n";
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
#			$redshift=$dat_now[7]/$vel_light;
			$red_em=$dat_now[7]/$vel_light;
			$SYS_VEL=$dat_now[7];
		    }
		}
		$ks++;
	    }
	}
	close(KS);
#	print "REDSHIFT_ELINES = $red_em , VEL = $SYS_VEL\n";
#	$red_em=$redshift;
#	<stdin>;
    }

    $d_redshift=0.0001;
    $fit_redshift=1;
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
	if ($is==0) {
	    set($vel_em_map,$ix,$iy,$info[7]);
	    set($disp_em_map,$ix,$iy,($info[5]/$info[1])*$info[7]);
	}
    }
    close(OUT_CONF);    


    $call="fit_spec_back.pl sub_spec.txt none ".$start_w_min." ".$end_w_max." none 0 none tmp_e.config ".$dev." > junk.junk";
    system($call);
#    print "$call\n";
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
    } else {
	for ($j=0;$j<$n_unc;$j++) {
	    $flux_masked[$j]=$flux_now*$masked[$j];
	    $model_joint[$j]=$model_spec_min[$j];
	    $res_joint[$j]=$res_spec[$j];
	}

    }

    open(OUT_MOD,">mod_joint_spec.txt");
    open(OUT_RES,">res_joint_spec.txt");
    for ($j=0;$j<$n_unc;$j++) {
	print OUT_MOD "$index_unc[$j] $wave_unc[$j] $model_joint[$j]\n";
	print OUT_RES "$index_unc[$j] $wave_unc[$j] $res_joint[$j]\n";
    }
    close(OUT_RES);
    close(OUT_MOD);

    for ($j=0;$j<$n_unc;$j++) {
	set($mod_cube,$ix,$iy,$j,$model_joint[$j]);
	set($res_cube,$ix,$iy,$j,$res_joint[$j]);
    }	


#    print "Press Enter"; <sdtin>;
    #
    #

#    $redshift=$red_abs;

    #    print "END****\n";
}
close(OUT);
#print "----------------------------\n";
#print "CHISQ=$MIN_CHISQ AGE=$age_min MET=$met_min AV=$Av REDSHIFT=$redshift SIGMA_DISP=$sigma\n";


open(FH,">compare_back_list_dust_mult.out");
print FH "$nx $ny $ix $iy $MIN_CHISQ $age_min $met_min $Av $redshift $sigma $FLUX $red_abs $red_em $median_FLUX $sigma_FLUX  $age_min_mass $met_min_mass\n";
close(FH);
	set($chi_map,$ix,$iy,$MIN_CHISQ);
	set($flux_map,$ix,$iy,$FLUX);
	set($age_map,$ix,$iy,$age_min);
	set($met_map,$ix,$iy,$met_min);
	set($age_mass_map,$ix,$iy,$age_min_mass);
	set($met_mass_map,$ix,$iy,$met_min_mass);
	set($dust_map,$ix,$iy,$Av);
	set($vel_abs_map,$ix,$iy,$redshift*$vel_light);
	set($disp_abs_map,$ix,$iy,($sigma*$redshift*$vel_light/4863));




	
$call="cat compare_back_list_dust_mult.out >> ".$outfile;
system($call);

	open(OELINES,">>$out_elines");
	print OELINES "$nx $ny $ix $iy\n";
	close(OELINES);
	$call="cat out.fit_spectra >> ".$out_elines;
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


#	print "PRESS "; <stdin>;
	}

	
    }
	$age_map->wfits($age_file);
	$met_map->wfits($met_file);
	$age_mass_map->wfits($age_mass_file);
	$met_mass_map->wfits($met_mass_file);
	$dust_map->wfits($dust_file);
	$chi_map->wfits($chi_file);
	$flux_map->wfits($flux_file);
	$vel_abs_map->wfits($vel_abs_file);
	$disp_abs_map->wfits($disp_abs_file);
	$vel_em_map->wfits($vel_em_file);
	$disp_em_map->wfits($disp_em_file);



}



$mod_cube->wfits($mod_cube_file);
$res_cube->wfits($res_cube_file);

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


    $header="NORM".$iii;
    $ml[$iii]=$pdl_flux_c_ini->hdr->{$header};
    


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
	    $error[$i][$iii]=1/(abs($val+0.001)+$chi_sec[$i]);
	  
#	    $error[$i][$iii]=10;
	} else {
	    $error[$i][$iii]=1e20;
	}

	for ($j=0;$j<$nline;$j++) {
	    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
		$error[$i][$iii]=1e10;
#		print "$wave_unc[$i] $error[$i][$iii]\n";
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
    if (($nf_new>0)&&($n_unc>0)) {
	while (($nf_neg>0)&&($nf_new>0)) {
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
	    if ($norm_C>0) {
#		$age_min=$age_min+$C*$age_mod[$k]/$norm_C;
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;


		if ($norm_C_mass>0) {
		    $age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		    #$age_min_mass=$age_min_mass+$C*$ml[$k]*$age_mod[$k]/$norm_C_mass;
		    $met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		}
	    }
	}
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
    }
    $name=$NAMEFILE.", (".$ix.",".$iy.") (".$nx.",".$ny.") ";
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
	pglabel("Wavelength","Counts","$name $Av $chi_sq $age_min ($age_min_mass) $met_min ($met_min_mass) $redshift $sigma");
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
