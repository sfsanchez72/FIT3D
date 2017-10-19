#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<7) {
    print "USE: extract_gauss_external.pl RAW.fits Spectral_axis[0/1] CEN.fits FWHM.fits OUTPUT.fits ALLOWED_SHIFT_CEN plot nplot\n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
#$trace_file=$ARGV[2];
$cen_file=$ARGV[2];
$fwhm_file=$ARGV[3];
$out_file=$ARGV[4];
#$nx_min=$ARGV[6];
#$nx_max=$ARGV[7];
$shift=$ARGV[5];
$plot=$ARGV[6];
$nplot=$ARGV[7];
$command="P";
print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
@tmp_array=read_img($infile);
print "Done\n";
if ($spec_axis==0) {
    @in_array=@tmp_array;
    $nx=$naxis[0];
    $ny=$naxis[1];
} else {
    for ($j=0;$j<$naxis[1];$j++) {
	for ($i=0;$i<$naxis[0];$i++) {
	    $in_array[$i][$j]=$tmp_array[$j][$i];
	}
    }
    $nx=$naxis[1];
    $ny=$naxis[0];
}
print "DONE\n";
print "Reading the Centroid an FWHM files\n";

@a_fwhm2=read_img($fwhm_file);
@a_cen=read_img($cen_file);
$nax_tr=read_naxes($fwhm_file);   
@naxis_tr=@$nax_tr;
#@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx_tr=$naxis_tr[0];
$ny_tr=$naxis_tr[1];

if ($nx_tr!=$nx) {
    print "The X-dimension of the trace file ($trace_file), $nx_tr\n";
    print "does not match that of the original image, $nx\n";
    exit;
}

print "Smoothing the tracing data...\n";
$sum=0;
for ($j=0;$j<11;$j++) {
    for ($i=0;$i<11;$i++) {
	$r=sqrt(($i-4)**2+($j-4)**2);
	$kernel[$i][$j]=exp(-0.5*($r/3)**2);
	$sum=$sum+$kernel[$i][$j];
#	print "$i $j $kernel[$i][$j]\n";
    }
}
$pdl_kernel=pdl(@kernel);
$pdl_kernel=$pdl_kernel/$sum;



for ($i=0;$i<$nx;$i++) {
    my @ar;
    for ($j=0;$j<$ny_tr;$j++) {	
	$ar[$j]=$a_fwhm2[$j][$i];
    }
    $mean=median(@ar);
    $sigma=sigma(@ar);
    for ($j=0;$j<$ny_tr;$j++) {	
	if (abs($ar[$j]-$mean)>$sigma) {
	    $a_fwhm2[$j][$i]=$mean;	    
#	    print "$a_fwhm2[$j][$i]=$mean\n";
	}
    }
}

$pdl_fwhm2=pdl(@a_fwhm2);
$smoothed=med2df($pdl_fwhm2,21,21,{Boundary => Reflect});
$smoothed2=conv2d($smoothed,$pdl_kernel,{Boundary => Reflect});
$smoothed3=med2df($smoothed2,5,5,{Boundary => Reflect});
for ($j=0;$j<$ny_tr;$j++) {	
    for ($i=0;$i<$nx;$i++) {
	$a_fwhm2[$j][$i]=$smoothed2->at($i,$j);
#	print "$i,$j $a_fwhm2[$j][$i]\n";
	}
}
system("rm mfwhm.fits");
$smoothed->wfits("mfwhm.fits");
$pdl_cen=pdl(@a_cen);
$smoothed=med2df($pdl_cen,3,3,{Boundary => Reflect});
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny_tr;$j++) {	
	$val=$smoothed->at($i,$j);
	$delta=abs($val-$a_cen[$j][$i]);
	if ($delta>$shift) {
	    $a_cen[$j][$i]=$smoothed->at($i,$j);
	}
    }
}

$pdl_cen2=pdl(@a_cen);
system("rm mcen.fits");
$pdl_cen2->wfits("mcen.fits");

print "DONE\n";
#
#
#
print "New Gaussian Determination\n";
#$command="P";
#$plot=1;
#for ($i=12;$i<$nx;$i++) {
for ($i=0;$i<$nx;$i++) {
#for ($i=$nx_min;$i<$nx_max;$i++) {
    $min=1e12;
    $max=-1e12;
    my @a;
    my @sig;
    my @xa;
    my @cen;
    my @fwhm;
    my @back;



#    $pdl_a=pdl(@a);
#    $pdl_xa=pdl(@xa);
#    $pdl_cen=pdl(@cen);
#    $pdl_fwhm=pdl(@fwhm);

#
# We find Alpha and Beta
#
    $pdl_val=zeroes($ny,$ny_tr);
#    my $pdl_beta=zeroes($ny_tr,$ny_tr);
    for ($jj=0;$jj<$ny_tr;$jj++) {
	$beta[$jj]=0;
	for ($j=0;$j<$ny;$j++) {	
	    $s=1;
#	    if ($in_array[$j][$i]>0) {
#		$s=sqrt(abs($in_array[$j][$i])); 
#	    } else {
#		$s=0.1;
#	    }
#	    $val[$jj][$j]=exp(-0.5*(($j*1.0-$a_cen[$jj][$i])/($a_fwhm2[$jj][$i]/2.345))**2);
	    $val=exp(-0.5*(($j*1.0-$a_cen[$jj][$i])/($a_fwhm2[$jj][$i]/2.345))**2);
#+$back[$k]
#	    if (abs($j-$a_cen[$jj][$i])<5) {
#		print "TEST=$j $a_cen[$jj][$i] $a_fwhm2[$jj][$i] $val[$jj][$j]\n";
#	    }
#my_gauss1d($j*1.0,$a_cen[$jj][$i],$a_fwhm2[$jj][$i])/$s;
#	    print $jj $j, $";
	    set($pdl_val,$j,$jj,$val);
	    $beta[$jj]=$beta[$jj]+$in_array[$j][$i]*$val/$s;
	    if ($max<$in_array[$j][$i]) {
		$max=$in_array[$j][$i];
	    }
	    if ($min>$in_array[$j][$i]) {
		$min=$in_array[$j][$i];
	    }
	}
#    print "$jj/$ny_tr\n";
#    <stdin>;
#	print "$jj $beta[$jj]\n";
    }

#    $pdl_beta=pdl(@beta);
#    $pdl_val=pdl(@val);
#    @dims=$pdl_val->dims;
#    print "*=@dims ($ny_tr,$ny)\n";
#    @dims=$pdl_val_test->dims;
#    print "*=@dims ($ny_tr,$ny)\n";
#    @dims=$pdl_beta->dims;
#    print "*=@dims\n";

#    $pdl_alpha=matmult($pdl_val,transpose($pdl_val));
    $pdl_trans=transpose($pdl_val);
    $pdl_val->wfits("pdl_val.fits");
    $pdl_trans->wfits("pdl_trans.fits");
#    $pdl_alpha=matmult($pdl_trans,$pdl_val);
    $pdl_alpha=matmult($pdl_val,$pdl_trans);
    $pdl_alpha->wfits("pdl_alpha.fits");
    eval {
	$pdl_in_alpha=matinv($pdl_alpha);
    };
#    for ($kk=0;$kk<$ny_tr;$kk++) {
#	for ($jj=0;$jj<$ny_tr;$jj++) {
#	    $v=$pdl_alpha->at($jj,$kk);
#	    print "$jj, $kk, $v\n";
#	}
#    }

#
# We find the solution
#
    my @new_pk;
    for ($j=0;$j<$ny_tr;$j++) {
	$new_pk[$j]=0;
	for ($k=0;$k<$ny_tr;$k++) {
#	    $v=$pdl_alpha->at($k,$j);
#	    $v=$pdl_in_alpha->at($k,$j);
	    $v=$pdl_in_alpha->at($j,$k);
#	    print "$j, $k, $v\n";
	    if ($v>1e-300) {
		$in_v=$v;
		$new_pk[$j]=$new_pk[$j]+$beta[$k]*$in_v;#[$k][$j];
#		print "$j $k $v $in_v $beta[$k] $new_pk[$j]\n";
	    }
	}
	$a_new_flux[$j][$i]=1.069*$new_pk[$j]*$a_fwhm2[$j][$i];
#	print "$j $a_cen[$j][$i] $a_fwhm2[$j][$i] $new_pk[$j]\n";
    }
#    print "$i/$nx\n";
    if (($plot==1)&&($command ne "A")) {

#	for ($jj=0;$jj<$ny_tr;$jj++) {
#	    for ($j=0;$j<$ny;$j++) {
#		if (abs($j-$a_cen[$jj][$i])<5) {
#		    $v=$pdl_val->at($j,$jj);
#		    print "$jj $j $a_cen[$jj][$i] $a_fwhm2[$jj][$i] $v\n";
#		}
#	    }
#	<stdin>;
#	}


	for ($j=0;$j<$ny;$j++) {
#	    $a_mod0[$j]=0;
	    $a_mod[$j]=0;
	    $xa[$j]=$j;
	    for ($k=0;$k<$ny_tr;$k++) {	
		$v=$pdl_val->at($j,$k);
#		$a_mod0[$j]=$a_mod0[$j]+$a_pk[$k][$i]*$v;#my_gauss1d($j,$a_cen[$k][$i],$a_fwhm2[$k][$i]);	    
		$a_mod[$j]=$a_mod[$j]+$new_pk[$k]*$v;#my_gauss1d($j,$a_cen[$k][$i],$a_fwhm2[$k][$i]);	    
	    }
	    $y[$j]=$in_array[$j][$i];
	    $x[$j]=$j;
	}
#	pgask(0);
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgsubp(1,$nplot);
	for ($ii=0;$ii<$nplot;$ii++) {
	    pgsci(1);
	    pgenv(($ny/$nplot)*$ii,($ny/$nplot)*($ii+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","$k/$ny_tr");
	    pgpoint($ny,\@x,\@y,3);
#	    pgsci(5);
#	    pgline($ny,\@x,\@a_mod0);
	    pgsci(8);
	    pgline($ny,\@x,\@a_mod);

	}
	pgclos();
	pgend();
	if ($command ne "a") {
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}
    }
    print "$i/$nx\n";
}



print "Writting the Aperture extracted flux $out_file\n";
system("rm $out_file");
write_fits($out_file,[$nx,$ny_tr],2,\@a_new_flux);




print "DONE\n";
exit;


sub my_gauss1d {
    my $x=@_[0];
    my $cen=@_[1];
    my $fwhm=@_[2];
    my $valgauss=exp(-0.5*(($x-$cen)/($fwhm/2.345))**2);
    return $valgauss;
}

