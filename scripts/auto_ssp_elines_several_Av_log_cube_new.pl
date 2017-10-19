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

   

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<5) {
    print "USE: auto_ssp_elines_several_Av_log_cube.pl SPEC1.RSS.fits BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
    print "CONFIG_FILE:\n";
    print "redshift delta_redshift min_redshift max_redshift DELTA_VEL RANGE_VEL DELTA_SIGMA(%) RANGE_SIGMA(%)\n";
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

open(TMP,"<$ARGV[4]");
$line=<TMP>;
chop($line);
@DATA=split(" ",$line);
close(TMP);

if ($#DATA==7) {
    $DV=$DATA[4];
    $RV=$DATA[5];
    $DS=$DATA[6];
    $RS=$DATA[7];
} else {
    $DV=10;
    $RV=30;
    $DS=0.05;
    $RS=0.2;

}

if ($#DATA==9) {
    $DV=$DATA[4];
    $RV=$DATA[5];
    $DS=$DATA[6];
    $RS=$DATA[7];
    $MIN_W=$DATA[8];
    $MAX_W=$DATA[8];
} else {
    $DV=10;
    $RV=30;
    $DS=0.05;
    $RS=0.2;
    $MIN_W=3800;
    $MAX_W=4300;
}

#print "$#DATA\n @DATA\n"; exit;


$infile=$ARGV[0];
print "Reading the Cube $infile\n";
$pdl=rfits($infile);
($nx,$ny,$nz)=$pdl->dims();
$crval=$pdl->hdr->{"CRVAL3"};
$cdelt=$pdl->hdr->{"CDELT3"};
$crpix=$pdl->hdr->{"CRPIX3"};




print "($nx,$ny,$nz)\n";
$outfile=$ARGV[2];
$h=$pdl->gethdr;
$mod_cube_file="mod_".$outfile.".fits";
$res_cube_file="res_".$outfile.".fits";
$gas_cube_file="gas_".$outfile.".fits";
$no_gas_cube_file="no_gas_".$outfile.".fits";

$mod_cube=zeroes($nx,$ny,$nz);
$res_cube=zeroes($nx,$ny,$nz);
$gas_cube=zeroes($nx,$ny,$nz);
$no_gas_cube=zeroes($nx,$ny,$nz);

$mod_cube->sethdr( $h );
$res_cube->sethdr( $h );
$gas_cube->sethdr( $h );
$no_gas_cube->sethdr( $h );






for ($j=0;$j<$#ARGV+1;$j++) {
    $PARAM[$j]=$ARGV[$j];
}
$NP=$#ARGV+1;
$PARAM[2]=" tmp.out";
$elines="elines_".$ARGV[2];
$call="rm -f ".$elines;
$coeffs="coeffs_".$ARGV[2];
$call="rm -f ".$elines;
system($call);
$call="rm -f auto.log";
system($call);
open(OUT,">$ARGV[2]");
open(FL,"<F.limit");
$FL=<FL>;
chop($FL);
close(FL);

$i_med=int($nx/2);
$j_med=int($ny/2);

$k1=int(0.4*$nz);
$k2=int(0.6*$nz);


my $a=$pdl->slice(",,$k1:$k2");
my $c=average($a->xchg(0,2));
my $map=$c->xchg(0,1);
$i1=$i_med-int(0.3*$i_med);
$i2=$i_med+int(0.3*$i_med);
$j1=$j_med-int(0.3*$j_med);
$j2=$j_med+int(0.3*$j_med);

$val_med=$map->at($i_med,$j_med);
for ($i=$i1;$i<$i2;$i++) {
    for ($j=$j1;$j<$j2;$j++) {
	$val=$map->at($i_med,$j_med);
	if (($val>$val_med)&&($val<1e30)) {
	    $val_med=$val;
	    $i_med=$i;
	    $j_med=$j;
	}
    }
}
print "PEAK at $i_med,$j_med\n";


#####################
# LOOP
#####################
print "DONE\n";
print "Processing the data\n";
LOOP();
print "DONE\n";

close(OUT);

print "Writting results\n";
$mod_cube->wfits($mod_cube_file);
$res_cube->wfits($res_cube_file);
$gas_cube->wfits($gas_cube_file);
$no_gas_cube->wfits($no_gas_cube_file);
print "DONE\n";


exit;



sub add_file_mod {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($mod_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}




sub add_file_res {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($res_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}

sub add_file_gas {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($gas_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}

sub add_file_no_gas {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($no_gas_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}


sub get_nearest_old {
    my $I=$_[0];
    my $J=$_[1];
    my $dist=1e12;
    my $file=$ARGV[2];
    my $line;
    my @near;
    open(IN,"<$file");
    while($line=<IN>) {
        chop($line);
	my @data=split(" ",$line);
	my $ix=$data[2];
	my $iy=$data[3];
	my $dist_now=sqrt(($I-$ix)**2+($J-$iy)**2);
	if (($dist_now<$dist)&&($dist_now>0)) {
	    $dist=$dist_now;
	    $near[0]=$data[8];
	    $near[1]=$data[9];
	}
    }
    close(IN);
    return @near;    
}

sub get_nearest {
    my $I=$_[0];
    my $J=$_[1];
    my $dist=1e12;
    my $file=$ARGV[2];
    my $line;
    my @near;
    open(IN,"<$file");
    my $sum=0;
    while($line=<IN>) {
        chop($line);
	my @data=split(" ",$line);
	my $ix=$data[2];
	my $iy=$data[3];
	my $chi=$data[4];
	my $dist_now=sqrt(($I-$ix)**2+($J-$iy)**2);
	if (($dist_now<2)&&($dist_now>0)&&($chi>0)) {
	    $dist=$dist_now;
	    $near[0]=$near[0]+$data[8]/(1+$dist);#/$chi;
	    $near[1]=$near[1]+$data[9]/(1+$dist);#/$chi;
	    $near[2]=$near[2]+$data[7]/(1+$dist);#/$chi;
	    $sum=$sum+1/(1+$dist);#/$chi;
	}
    }
    if ($sum>0) {
	$near[0]=$near[0]/$sum;
	$near[1]=$near[1]/$sum;
	$near[2]=$near[2]/$sum;
    }
    close(IN);
    return @near;    
}


sub get_nearest_last {
    my $I=$_[0];
    my $J=$_[1];
    my $dist=1e12;
    my $file=$ARGV[2];
    my $line;
    my @near;
    open(IN,"<$file");
    my $sum=0;
    while($line=<IN>) {
        chop($line);
	my @data=split(" ",$line);
	my $ix=$data[2];
	my $iy=$data[3];
	my $chi=$data[4];
	$near[0]=$data[8];#/$chi;
	$near[1]=$data[9];#/$chi;
    }

    close(IN);
    return @near;    
}



sub DO {
    my $i=$_[0];
    my $j=$_[1];
    if (($i>0)&&($i<$nx)&&($j>0)&&($j<$ny)) {

    if ($k>0) {
	@near=get_nearest($i,$j);
	$PARAM[11]=$near[0];
	if ($PARAM[12]!=0) {
	    $PARAM[12]=$DV/300000;
	}
	$PARAM[13]=$near[0]-$RV/300000;
	$PARAM[14]=$near[0]+$RV/300000;
	if ($ARGV[13]>$PARAM[13]) {
	    $PARAM[13]=$ARGV[13];
	}
	if ($ARGV[14]<$PARAM[14]) {
	    $PARAM[14]=$ARGV[14];
	}

	$PARAM[15]=$near[1];
	if ($PARAM[16]!=0) {
	    $PARAM[16]=$DS*$PARAM[15];
	}
	$PARAM[17]=$near[1]*(1-$RS);
	$PARAM[18]=$near[1]*(1+$RS);
	if ($ARGV[17]>$PARAM[17]) {
	    $PARAM[17]=$ARGV[17];
	}
	if ($ARGV[18]<$PARAM[18]) {
	    $PARAM[18]=$ARGV[18];
	}

	if ($noo>0) {
	    if ($PARAM[18]>1.3*$mean_sigma) {
		$PARAM[18]=1.3*$mean_sigma;
	    }
	}

	$PARAM[19]=$near[2];
	if ($PARAM[20]!=0) {
	    $PARAM[20]=$DS*$PARAM[19];
	}
	$PARAM[21]=$near[2]*(1-$RS);
	$PARAM[22]=$near[2]*(1+$RS);
	if ($ARGV[21]>$PARAM[21]) {
	    $PARAM[21]=$ARGV[21];
	}
	if ($ARGV[22]<$PARAM[22]) {
	    $PARAM[22]=$ARGV[22];
	}


    }
    $k++;
    $nz_med=int($nz/2);
    $F=$pdl->at($i,$j,$nz_med);
    if (($F>$FL)&&($F<1e30)) {
	$call="rm -f spec.auto.txt";
	system($call);
	$call="rm -f tmp.out";
	system($call);
	open(SPEC,">spec.auto.txt");
	for ($I=0;$I<$nz;$I++) {
	    $wave=$crval+$cdelt*($I+1-$crpix);
	    $val=$pdl->at($i,$j,$I);
	    print SPEC "$I $wave $val\n";
	}
	close(SPEC);
#	$call="get_spec_cube.pl ".$infile." ".$i." ".$j." spec.auto.txt";
#	system($call);
	$call="auto_ssp_elines_several_Av_log_new.pl spec.auto.txt ";
	for ($jj=1;$jj<$NP;$jj++) {
	    $call=$call." ".$PARAM[$jj];
	}
	$call=$call." >> auto.log";
	system($call);
#	print "$call\n"; <stdin>; exit;
	open(FH,"<tmp.out");
	$line=<FH>;
	if ($i==0) {
	    print OUT "$line";
	}
	    $line=<FH>;
	chop($line);
	print OUT "$nx $ny $i $j $line\n";
	close(FH);
	open(ELINES,">>$elines");
	$ne=0;
	open(IN,"<elines_tmp.out");
	while($in=<IN>) {
	    if ($in =~ "eline") {
		$line_in[$ne]=$in;
		$ne++;
	    }
	}
	close(IN);
	print ELINES "$nx $ny $i $j\n";
	print ELINES "$ne\n";
	for ($ii=0;$ii<$ne;$ii++) {
	    print ELINES "$line_in[$ii]";
	}
	close(ELINES);
	$call="cat coeffs.out >> ".$coeffs;
	system($call);
	
	
	
	add_file_mod($i,$j,"org_spec.txt");
	add_file_res($i,$j,"res_joint_spec.txt");
	add_file_gas($i,$j,"gas_spec.txt");
	add_file_no_gas($i,$j,"no_gas_spec.txt");
	

	
	
	
	
	print "$i/$nx $j/$ny DONE\n";
    } else {
	print "$i/$nx $j/$ny PASSED BY\n";
    }
    
    return 1;
    } else {
	return 0;
    }
	
}


sub LOOP_SQ {
# 1st Quadrant
#
    $k=0;
    for ($j=$j_med;$j<$ny;$j++) {
	for ($i=$i_med;$i<$nx;$i++) {
	    DO($i,$j);
	}
    }
    print "1st Quadrant\n";
    
# 2nd Quadrant
#
#$k=0;
    for ($j=$j_med-1;$j>-1;$j--) {
	for ($i=$i_med;$i<$nx;$i++) {
	    DO($i,$j);
	}
    }
    print "2nd Quadrant\n";
    
# 3nd Quadrant
#
#$k=0;
    for ($j=$j_med-1;$j>-1;$j--) {
	for ($i=$i_med-1;$i>-1;$i--) {
	    DO($i,$j);
	}
    }
    print "3rd Quadrant\n";
# 4th Quadrant
#
#$k=0;
    for ($j=$j_med;$j<$ny;$j++) {
	for ($i=$i_med-1;$i>-1;$i--) {
	    DO($i,$j);
	}
    }
    print "4th Quadrant\n";
    
    $done=1;
    return $done;
}

sub LOOP {

    my $nmax=($nx-1)*($ny-1);
    my $n_now=0;
    my $delta=2;
    $k=0;
    my $js=$j_med;
    my $is=$i_med;
    my $sig=0;
    $i=$is;
    $j=$js;
    $s=DO($i,$j);
    $n_now=$now+$s;
    while ($n_now<$nmax) {
	if ($sig==0) {
	    for ($j=$js+1;$j<$js+$delta;$j++) {
		$s=DO($i,$j);
		$n_now=$n_now+$s;
	    }
	    $j--;
	    for ($i=$is+1;$i<$is+$delta;$i++) {
		$s=DO($i,$j);
		$n_now=$n_now+$s;
	    }
	    $sig=1;
	    $i--;
	} else {
	    for ($j=$js-1;$j>$js-$delta;$j--) {
		$s=DO($i,$j);
		$n_now=$n_now+$s;
	    }
	    $j++;
	    for ($i=$is-1;$i>$is-$delta;$i--) {
		$s=DO($i,$j);
		$n_now=$n_now+$s;
	    }
	    $sig=0;
	    $i++;
	}
	$delta++;
	$js=$j;
	$is=$i;
	$gas_cube->wfits($gas_cube_file);
	$noo=0;
	open(OO,"<$outfile");
	while($oo=<OO>) {
	    chop($oo);
	    if ($oo !~ "#") {
		@doo=split(" ",$oo);
		$poo[$noo]=$doo[9];
		if ($poo[$noo]>0) {
		    $noo++;
		}
	    }
	}
	close(OO);
	$mean_sigma=mean(@poo);
	print "$n_now/$nmax\n";
    }
    $done=1;
    return $done;
}
