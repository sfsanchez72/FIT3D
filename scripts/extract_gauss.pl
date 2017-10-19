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
#use PDL::Fit::Linfit;
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<8) {
    print "USE: extract_gauss.pl RAW.fits Spectral_axis[0/1] TRACE.fits WIDTH OUTPUT.fits NX_MIN NX_MAX plot nplot\n";
    exit;
}

$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$trace_file=$ARGV[2];
$aperture=$ARGV[3];
#$fix=$ARGV[4];
$out_file=$ARGV[4];
$nx_min=$ARGV[5];
$nx_max=$ARGV[6];
$plot=$ARGV[7];
$nplot=$ARGV[8];
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
print "Reading file $trace_file ";
$nax_tr=read_naxes($trace_file);   
@naxis_tr=@$nax_tr;
@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx_tr=$naxis_tr[0];
$ny_tr=$naxis_tr[1];

if ($nx_tr!=$nx) {
    print "The X-dimension of the trace file ($trace_file), $nx_tr\n";
    print "does not match that of the original image, $nx\n";
    exit;
}

print "Determining the parameters of the gaussians...\n";
for ($i=$nx_min;$i<$nx_max;$i++) {
#
# We create a cut for each peak
#
#$i=int($nx/2);
for ($k=0;$k<$ny_tr;$k++) {	
    $xk[$k]=$k;
    $range=int(2*$aperture);
    $j_min=$peak_y_max_new[$k][$i]-$range;
    $j_max=$peak_y_max_new[$k][$i]+$range;
	if ($j_min<0) {
	    $j_min=0;
	}
    if ($j_max>=$ny) {
	$j_max=$ny-1;
    }
    $min=1e12;
    $max=-1e12;
    my @a;
    my @xa;
    $na=0;
    for ($j=$j_min;$j<$j_max;$j++) {	
	$a[$na]=$in_array[$j][$i];
	if ($min>$a[$na]) {
	    $min=$a[$na];
	}
	if ($max<$a[$na]) {
	    $max=$a[$na];
	}
	$xa[$na]=$j;
	$na++;
    }
    $pdl_xa=pdl(@xa);
    $pdl_a=pdl(@a);
    ($pcen, $ppk, $pfwhm2, $pback, $perr, $pfit)=fitgauss1d($pdl_xa,$pdl_a);
    $cen[$k]=$pcen->slice(0)->sclr;
    $pk[$k]=$ppk->slice(0)->sclr;
    $fwhm2[$k]=$pfwhm2->slice(0)->sclr;
    $back[$k]=$pback->slice(0)->sclr;
    $err[$k]=$perr->slice(0)->sclr;

    if ($fwhm2[$k]>50) {
	$fwhm2[$k]=50;
    }

    $a_fwhm2[$k][$i]=$fwhm2[$k];
    $a_cen[$k][$i]=$cen[$k];
    $a_flux[$k][$i]=1.069*$pk[$k]*$fwhm2[$k];
    $a_pk[$k][$i]=$pk[$k];
    $a_back[$k][$i]=$back[$k];
    


#    print "($cen, $pk, $fwhm2, $back, $err, $fit)";
#    ($cen, $pk, $fwhm2, $back, $err, $fit) 
    if (($plot==1)&&($command ne "A")) {
	my @a_mod;
	for ($j=0;$j<$na;$j++) {
	    $a_mod[$j]=$pk[$k]*exp(-0.5*(($xa[$j]-$cen[$k])/($fwhm2[$k]/2.345))**2)+$back[$k]
	}


	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($j_min,$j_max,$min,$max,0,0);
	pglabel("Y-axis","Counts","$k/$ny_tr");
	pgpoint($na,\@xa,\@a,3);
	pgsci(8);
	pgline($na,\@xa,\@a_mod);
	pgclos();
	if ($command ne "a") {
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}
    }
}
#$command="P";
$med_fwhm2=median(@fwhm2);
$sig_fwhm2=sigma(@fwhm2);
if ($i==50*int($i/50)) {
    print "($i/$nx) FWHM=$med_fwhm2+-$sig_fwhm2\n";
}
#$npoly=3;
#$xk_pdl = pdl(@xk);
#$fwhm2_pdl = pdl(@fwhm2);
#($s_y,$coeff) = fitpoly1d $xk_pdl,$fwhm2_pdl,$npoly;
#for ($j=0;$j<$npoly;$j++) {
#    $c[$j]=$coeff->slice($j)->sclr;
#}
#for ($j=0;$j<$ny_tr;$j++) {
#    $fwhm2_out[$j]=0;
#    for ($kk=0;$kk<$npoly;$kk++) {
#	$fwhm2_out[$j]=$fwhm2_out[$j]+$c[$kk]*($j**$kk);
#    }
#}
@fwhm2_out=median_filter(6,\@fwhm2);
if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$ny_tr,$med_fwhm2-3*$sig_fwhm2,$med_fwhm2+3*$sig_fwhm2,0,0);
	pglabel("Spectral Peaks","FWHM","");
	pgpoint($ny_tr,\@xk,\@fwhm2,3);
	pgsci(2);
	pgline($ny_tr,\@xk,\@fwhm2_out);
	pgsci(8);
	pgline(2,[0,$ny_tr],[$med_fwhm2,$med_fwhm2]);
	pgclos();
	if ($command ne "a") {
	    print "Press Enter"; 
	    $command=<stdin>;
	    chop($command);
	}
    }

#print "Press Enter"; <stdin>;
}

for ($i=0;$i<$nx_min;$i++) {
    for ($k=0;$k<$ny_tr;$k++) {	
	$a_fwhm2[$k][$i]=$a_fwhm2[$k][$nx_min];
	$a_cen[$k][$i]=$a_cen[$k][$nx_min];
	$a_flux[$k][$i]=$a_flux[$k][$nx_min];
	$a_back[$k][$i]=$a_back[$k][$nx_min];    
    }
}

for ($i=$nx_max;$i<$nx;$i++) {
    for ($k=0;$k<$ny_tr;$k++) {	
	$a_fwhm2[$k][$i]=$a_fwhm2[$k][$nx_max-1];
	$a_cen[$k][$i]=$a_cen[$k][$nx_max-1];
	$a_flux[$k][$i]=$a_flux[$k][$nx_max-1];
	$a_back[$k][$i]=$a_back[$k][$nx_max-1];    
    }
}

#$a_mfwhm2=median_smooth(@a_fwhm2,3,3);
print "DONE\n";
print "Writting fwhm.fits, cen.fits & bank.fits\n";
system("rm fwhm.fits");
#write_fits("fwhm.fits",[$nx,$ny_tr],2,\@a_mfwhm2);
write_fits("fwhm.fits",[$nx,$ny_tr],2,\@a_fwhm2);
system("rm cen.fits");
write_fits("cen.fits",[$nx,$ny_tr],2,\@a_cen);
system("rm back.fits");
write_fits("back.fits",[$nx,$ny_tr],2,\@a_back);
system("rm flux.fits");
write_fits("flux.fits",[$nx,$ny_tr],2,\@a_flux);
$pdl_fwhm2=pdl(@a_fwhm2);
$smoothed=med2df($pdl_fwhm2,5,5,{Boundary => Reflect});
for ($j=0;$j<$ny_tr;$j++) {	
    for ($i=0;$i<$nx;$i++) {
	$a_fwhm2[$j][$i]=$smoothed->at($i,$j);
	}
}
print "DONE\n";
#
#
#
print "New Gaussian Determination\n";
#$command="P";
#$plot=1;
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
#	print "$j $new_pk[$j]\n";
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
	    $a_mod0[$j]=0;
	    $a_mod[$j]=0;
	    $xa[$j]=$j;
	    for ($k=0;$k<$ny_tr;$k++) {	
		$v=$pdl_val->at($j,$k);
		$a_mod0[$j]=$a_mod0[$j]+$a_pk[$k][$i]*$v;#my_gauss1d($j,$a_cen[$k][$i],$a_fwhm2[$k][$i]);	    
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
	    pgsci(5);
	    pgline($ny,\@x,\@a_mod0);
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
    if ($i==50*int($i/50)) {
	print "($i/$nx)\n";
    }
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

