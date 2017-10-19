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
use PDL::Func;
use PDL::Math;

use PDL::Fit::Linfit;




$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<9) {
    print "USE: disp_cor_error.pl EXTRACTED.fits CRVAL CDELT APERTURE NY_SPECTRA WITDH NPOLY OUTPUT.fits DISPERSION.txt plot [NMAX]\n";
    exit;
}

$infile=$ARGV[0];
$crval=$ARGV[1];
$cdelt=$ARGV[2];
$aperture=$ARGV[3];
$ny_cen=$ARGV[4];
$ny_width=$ARGV[5];
$npoly=$ARGV[6];
$out_file=$ARGV[7];
$disp_file=$ARGV[8];
$plot=$ARGV[9];
$nmax=0;
if ($#ARGV==10) {
    $nmax=$ARGV[10];
}
$command="";
$ask_file=0;
if ($#ARGV==11) {
    $nmax=$ARGV[10];
    $command="F";    
    $ask_file=1;
    $file_load=$ARGV[11];
    $plot=0;
}

print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";


#
# We construct the checking spectrum
#
$j_min=$ny_cen-$ny_width;
$j_max=$ny_cen+$ny_width;
if ($j_min<0) {
    $j_min=0;
}
if ($j_max>$ny) {
    $j_max=$ny;
}
$min=1e12;
$max=-1e12;
my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    $k=0;
    for ($j=$j_min;$j<$j_max;$j++) {
	$cut[$k]=$in_array[$j][$i];
	$k++;
    }
    $w[$i]=$i*1.0;
    $spec[$i]=median(@cut);
    if ($min>$spec[$i]) {
	$min=$spec[$i];
    }
    if ($max<$spec[$i]) {
	$max=$spec[$i];
    }
}

@borders_ini=(0,$nx,$min,$max);
@borders=(0,$nx,$min,$max);

#
# Main Loop to select the features.
#
$QUIT="N";
my @xp;
my @yp;
my @wp;
my @mp;
$np=0;
if ($plot==1) {
    $LOOP="Q";
}

do {
    print "a => Unzoom\n";
    print "x => X-Zoom\n";
    print "y => Y-Zoom\n";
    print "m => Mark a peak\n";
    print "d => delete a peak\n";
    print "r => Move right\n";
    print "l => Move left\n";
    print "R => Load a file\n";
    print "F => Load and recenter a file\n";
    print "W => Save a file\n";    
    print "Q => Quit\n";

    while($QUIT eq "N") {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($borders[0],$borders[1],$borders[2],$borders[3],0,1);
	pgsci(1);
	pgsci(7);
	pgline($nx,\@w,\@spec);
	pgsch(0.7);
	pgsci(5);
	for ($i=0;$i<$np;$i++) {
	    if ($mp[$i]==1) {
		pgpoint(1,[$xp[$i]],[$yp[$i]],2);
		pgptxt($xp[$i],1.5*$yp[$i],0,0.5,$wp[$i]);
	    }
	}
	pgsch(1.2);
	pgsci(1);
	
	while ($command ne "X") {
	    if ($ask_file==0) {
		pgband(7,0,50,50,$x,$y,$command);
	    }
	    if ($command eq "x") {
		$x_cen=0.5*($borders[1]-$borders[0])+$borders[0];
		$x_delta=0.5*($borders[1]-$borders[0]);
		$x_cen=$x;
		$x_delta=0.75*$x_delta;
		$borders[0]=$x_cen-$x_delta;
		$borders[1]=$x_cen+$x_delta;
		$command="X";
	    }
	    if ($command eq "y") {
	    $y_cen=0.5*($borders[3]-$borders[2])+$borders[2];
	    $y_delta=0.5*($borders[3]-$borders[2]);
	    $y_cen=$y;
	    $y_delta=0.75*$y_delta;
	    $borders[2]=$y_cen-$y_delta;
	    $borders[3]=$y_cen+$y_delta;
	    $command="X";
	}
	if ($command eq "r") {
	    $x_cen=0.5*($borders[1]-$borders[0])+$borders[0];
	    $x_delta=0.5*($borders[1]-$borders[0]);
	    $x_cen=$x_cen+0.25*$x_delta;
	    $borders[0]=$x_cen-$x_delta;
	    $borders[1]=$x_cen+$x_delta;
	    $command="X";
	}
	if ($command eq "l") {
	    $x_cen=0.5*($borders[1]-$borders[0])+$borders[0];
	    $x_delta=0.5*($borders[1]-$borders[0]);
	    $x_cen=$x_cen-0.25*$x_delta;
	    $borders[0]=$x_cen-$x_delta;
	    $borders[1]=$x_cen+$x_delta;
	    $command="X";
	}
	if ($command eq "a") {
	    @borders=@borders_ini;
	    $command="X";
	}
	if ($command eq "q") {
	    $command="X";
	    $QUIT="Y";
	}
	if ($command eq "R") {
	    print "Identification file?"; $file=<stdin>;
	    chop($file);
	    $np=0;
	    open(FH,"<$file");
	    while($line=<FH>) {
		@data=split(" ",$line);
		$xp[$np]=$data[0];
		$wp[$np]=$data[1];
		$mp[$np]=1;
		$yp[$np]=$spec[int($data[0])];
		$np++;
	    }
	    close(FH);
	    $command="X";
	}
	if ($command eq "F") {
	    if ($ask_file==0) {
		print "Identification file (force)?"; $file=<stdin>;		
		chop($file);
	    } else {
		$file=$file_load;
#		$ask_file=1;
	    }

	    $np=0;
	    $nt=0;
	    open(FH,"<$file");
	    while($line=<FH>) {
		@data=split(" ",$line);
		$xt[$nt]=$data[0];
		$wt[$nt]=$data[1];
		$mt[$nt]=1;
		$yt[$nt]=$spec[int($data[0])];
		$nt++;
	    }
#	    $np=$nt;
	    close(FH);
	    $command="X";

#
#   We look for the new maximum
# 
	    $np=0;
	    for ($j=0;$j<$nt;$j++) {
		my $x=$xt[$j];
		my $NNN=3;
		my $i_min=int($x)-$aperture;
		my $i_max=int($x)+$aperture;
		my $maximum=-1e12;		
		for ($i=$i_min;$i<$i_max;$i++) {
		    if (($spec[$i]>$maximum)&&($spec[$i]>$spec[$i-1])&&($spec[$i]>$spec[$i+1])) {
			$i_peak=$i;
			$maximum=$spec[$i];
		    }
		}
		if ($maximum>0) {
		    $a=$i_peak-1;
		    $b=$i_peak;
		    $c=$i_peak+1;
		    $fa=-$spec[$a];
		    $fb=-$spec[$b];
		    $fc=-$spec[$c];
		    $den=($fc-2*$fb+$fa);
		    if ($den!=0) {
			$x_peak=$c-($b-$a)*(($fc-$fb)/$den+0.5);
		    } else {
			$x_peak=$i_peak;
		    }

		    $x_peak_mean=0;
		    $peak_mean=0;
		    $n_width=3;
		    $i=int($x);
		    for ($jj=$i-$n_width;$jj<=$i+$n_width;$jj++) {
			$x_peak_mean=$x_peak_mean+$spec[$jj]*($jj);
			$peak_mean=$peak_mean+$spec[$jj];
		    }
		    $x_peak_mean=$x_peak_mean/$peak_mean;
#		    print "X_PEAK = $x_peak $x_peak_mean\n";


#		    $x_peak=0.5*($x_peak+$x_peak_mean);

#		    if ($np==0) {
#			$x_peak=$x_peak+0.25;
#		    }

		    $xp[$np]=$x_peak;
		    $yp[$np]=$spec[$i_peak];
		    $mp[$np]=1;
		    $wp[$np]=$wt[$j];
		    #print "$np $xp[$np] $yp[$np] $wp[$np]\n";
		    

		    $np++;
		}
		pgsch(0.7);
		pgsci(5);
		for ($i=0;$i<$np;$i++) {
		    if ($mp[$i]==1) {
			pgpoint(1,[$xp[$i]],[$yp[$i]],2);
			pgptxt($xp[$i],1.5*$yp[$i],0,0.5,$wp[$i]);
		    }
		}
		pgsch(1.2);		    
	    }
	    $command="X";

	    if ($ask_file==1) {
		$QUIT="Y";
		$LOOP="N";
	    }

	}
    
   

	if ($command eq "W") {
	    print "Saving Identification file?"; $file=<stdin>;
	    chop($file);
	    open(FH,">$file");
	    for ($kk=0;$kk<$np;$kk++) {
		if ($mp[$kk]==1) {
		    print FH "$xp[$kk] $wp[$kk]\n";
		}
	    }
	    close(FH);
	    $command="X";
	}

	if ($command eq "V") {
	    for ($kk=0;$kk<$np;$kk++) {
		if ($mp[$kk]==1) {
		    print"$xp[$kk] $wp[$kk]\n";
		}
	    }
	    $command="X";
	}

	if ($command eq "m") {
	    $xs=$x;
	    $ys=$spec[int($xs)];
	    pgsci(8);
	    pgpoint(1,[$xs],[$ys],2);
	    pgsci(1);
	    $i_min=int($xs)-$aperture;
	    $i_max=int($xs)+$aperture;
	    $i_peak=int($xs);
	    for ($i=$i_min;$i<$i_max;$i++) {
		if (($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])) {
		    $i_peak=$i;
		}
	    }
	    $a=$i_peak-1;
	    $b=$i_peak;
	    $c=$i_peak+1;
	    $fa=-$spec[$a];
	    $fb=-$spec[$b];
	    $fc=-$spec[$c];
	    $den=($fc-2*$fb+$fa);
	    if ($den!=0) {
		$x_peak=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	    } else {
		$x_peak=$i_peak;
	    }
	    $xp[$np]=$x_peak;
	    $yp[$np]=$spec[$i_peak];
	    $mp[$np]=1;
	    pgsch(0.7);
	    pgsci(5);
	    for ($i=0;$i<$np;$i++) {
		if ($mp[$i]==1) {
		    pgpoint(1,[$xp[$i]],[$yp[$i]],2);
		    pgptxt($xp[$i],1.5*$yp[$i],0,0.5,$wp[$i]);
		}
	    }
	    pgsch(1.2);

	    print "Wavelength at $x_peak?"; $data=<stdin>;
	    chop($data);
	    $wp[$np]=$data;
	    $np++;
	}

	if ($command eq "d") {
	    $xs=$x;
	    $ys=$spec[int($xs)];
	    $i_min=int($xs)-$aperture;
	    $i_max=int($xs)+$aperture;
	    $i_peak=int($xs);
	    for ($i=$i_min;$i<$i_max;$i++) {
		if (($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])) {
		    $i_peak=$i;
		}
	    }
	    for ($i=0;$i<$np;$i++) {
		if (abs($i_peak-int($xp[$i]))<2) {
		    $mp[$i]=0;
		}
#		print "$i_peak $xp[$i] $mp[$i]\n";
	    }
	    pgsch(0.7);
	    for ($i=0;$i<$np;$i++) {
		if ($mp[$i]==1) {
		    pgsci(5);
		    pgpoint(1,[$xp[$i]],[$yp[$i]],2);
		    pgptxt($xp[$i],1.5*$yp[$i],0,0.5,$wp[$i]);
		} else {
		    pgsci(2);
		    pgpoint(1,[$xp[$i]],[$yp[$i]],22);
		}
	    }
	    pgsci(1);
	    pgsch(1.2);
	}

	

	}
#	print "$command $QUIT\n";
	    $command="";
    pgclos();
    pgend();

    }

#print "DONE\n";
#
# We determine now the pixels in the new coordinate system:
#
$n=0;
$x_max=0;
for ($i=0;$i<$np;$i++) {
    if ($mp[$i]==1) {
	$x_old[$n]=$xp[$i];	
	$x_new[$n]=($wp[$i]-$crval)/$cdelt;
	$y_old[$n]=$yp[$i];
#	print "WWW $wp[$i] '$crval' '$cdelt' $x_new[$n] $i $n\n";
	if ($x_max<$x_new[$n]) {
	    $x_max=$x_new[$n]+10;
	}
	$n++;
    }
}

    if ($npoly==0) {
	my $spline=new Math::Spline(\@x_new,\@x_old);
	$i_old_max=$spline->evaluate($nx);
	$n_cut=0;
	$i=0;
	while ($n_cut==0) {
	    $a_i_new[$i]=$i;
	    $a_i_old[$i]=$spline->evaluate($i*1.000);
	    if ($a_i_old[$i]>$nx) {
		$n_cut=$i;
	    }	    
	    $i++;	    
	}
	
	if ($nmax>0) {
	    $n_cut=$nmax;
	}
	print "N_MAX=$n_cut\n";
#
# RMS
#
	$rms=0;
	for ($i=0;$i<$n;$i++) {
	    $x_new_old[$i]=$spline->evaluate($x_new[$i]);
	    $rms_now[$i]=substr(abs($x_new_old[$i]-$x_old[$i]),0,5);
	    $rms=$rms+($x_new_old[$i]-$x_old[$i])**2;
	}
	$rms=sqrt($rms/$n);
	print "RMS=$rms\n";
    } else {
#	print "NPOLY=$npoly\n";
	if ($npoly>0) {
	    $pdl_x_old=pdl(@x_old);
	    $pdl_x_new=pdl(@x_new);

	    $pdl_y_old=pdl(@y_old);
	    $pdl_func=zeroes($n,$npoly);
	    $pdl_error=zeroes($n,$npoly);
	    for ($i=0;$i<$npoly;$i++) {
		my $pdl_now=$pdl_func->slice(":,$i");
		$pdl_now .= ($pdl_x_new)**$i;
		my $pdl_now=$pdl_error->slice(":,$i");
		$pdl_now .= 1/sqrt(abs($pdl_y_old));
	    }
	    ($s_y,$COEFF) = linfit1d ($pdl_x_new,$pdl_x_old,$pdl_func, {Weights => $pdl_error});
	    $coeff=$COEFF->slice(":,(0)");

	    for ($j=0;$j<$npoly;$j++) {
		$c[$j]=$coeff->at($j);
	    }
	    $i_old_min=int($c[0]);
	    $i_old_max=0;
	for ($kk=0;$kk<$npoly;$kk++) {
	    $i_old_max=$i_old_max+$c[$kk]*($nx**$kk);
	}    
	    $n_cut=0;
	    $i=0;
	    while ($n_cut==0) {
		$a_i_new[$i]=$i;
		$a_i_old[$i]=0;
		for ($kk=0;$kk<$npoly;$kk++) {
		    $a_i_old[$i]=$a_i_old[$i]+$c[$kk]*(($i*1.000)**$kk);
		}
		if ($a_i_old[$i]>$nx) {
		    $n_cut=$i;
		}
		$i++;
	    }
	    if ($nmax>0) {
		$n_cut=$nmax;
	    }
	    print "N_MAX=$n_cut\n";
#
# RMS
#
	    $rms=0;
	    for ($i=0;$i<$n;$i++) {
		$x_new_old[$i]=0;
		for ($kk=0;$kk<$npoly;$kk++) {
		    $x_new_old[$i]=$x_new_old[$i]+$c[$kk]*(($x_new[$i])**$kk);
		}
		$rms_now[$i]=substr(abs($x_new_old[$i]-$x_old[$i]),0,5);
		$rms=$rms+($x_new_old[$i]-$x_old[$i])**2;
	    }
	    $rms=sqrt($rms/$n);
	    print "RMS=$rms\n";
	} else {
# LEGENDRE
	    $npoly2=abs($npoly);
	    $pdl_x_old=pdl(@x_old);
	    $pdl_x_new=pdl(@x_new);
	    $pdl_y_old=pdl(@y_old);
	    $pdl_func=zeroes($n,$npoly2+1);
	    $pdl_error=zeroes($n,$npoly2+1);
	    for ($i=0;$i<$npoly2+1;$i++) {
		my $pdl_now=$pdl_func->slice(":,$i");
		$pdl_now .= legendre($i,$pdl_x_new);
		my $pdl_now=$pdl_error->slice(":,$i");
		$pdl_now .= 1/sqrt(abs($pdl_y_old));
	    }

	    ($s_y,$COEFF) = linfit1d ($pdl_x_new,$pdl_x_old,$pdl_func, {Weights => $pdl_error});
	    $coeff=$COEFF->slice(":,(0)");
#	    ($s_y,$coeff) = fit_legendre($pdl_x_new,$pdl_x_old,$npoly2);
	    for ($j=0;$j<$npoly2+1;$j++) {
		$c[$j]=$coeff->at($j);		
	    }
	    $i_old_min=int($c[0]);
	    $i_old_max=0;
	    for ($kk=0;$kk<$npoly2+1;$kk++) {
#	    for ($kk=0;$kk<$npoly2;$kk++) {
		$i_old_max=$i_old_max+$c[$kk]*legendre($kk,$nx);
	    }    
	    $n_cut=0;
	    $i=0;
	    while ($n_cut==0) {
		$a_i_new[$i]=$i;
		$a_i_old[$i]=0;
		my $val_now=0;
		for ($kk=0;$kk<$npoly2+1;$kk++) {
#		    $a_i_old[$i]=$a_i_old[$i]+legendre($kk,($i*1.000))*($coeff->at($kk));
		    $val_now=$val_now+legendre($kk,($i))*($coeff->at($kk));
		    $CC=$coeff->at($kk);
#		    print "$kk $CC\n";
#		    $a_i_old[$i]=apr_n($a_i_old[$i],3);
		}
		$a_i_old[$i]=$val_now;
		if ($a_i_old[$i]>$nx) {
		    $n_cut=$i;
		}

#		print "$i $a_i_old[$i] $S_Y\n";
		$i++;
	    }
	    if ($nmax>0) {
		$n_cut=$nmax;
	    }
	    print "N_MAX=$n_cut\n";
#
# RMS
#
	    $rms=0;
	    for ($i=0;$i<$n;$i++) {
		$x_new_old[$i]=0;
		for ($kk=0;$kk<$npoly2+1;$kk++) {
		    $x_new_old[$i]=$x_new_old[$i]+$c[$kk]*legendre($kk,$x_new[$i]);
		}
		$rms_now[$i]=substr(abs($x_new_old[$i]-$x_old[$i]),0,5);
		$rms=$rms+($x_new_old[$i]-$x_old[$i])**2;
#		print  "$i $x_new_old[$i] $x_old[$i] $rms_now[$i]\n";

#		$S_Y=$s_y->at($i);		
	    }
	    $rms=sqrt($rms/$n);
	    print "RMS=$rms\n";

	}
    }


if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv(0,$nx,0,$n_cut,0,0);
    pglabel("X old","X new","Dispersion solution");
    pgsci(1);
    pgpoint($n,\@x_old,\@x_new,2);
    pgsci(2);
    pgline($n_cut,\@a_i_old,\@a_i_new);
    pgsch(0.7);
    for ($i=0;$i<$n;$i++) {
	pgsci(5);
	pgptxt($x_old[$i],$x_new[$i],0,0.5,$rms_now[$i]);
    }
    pgsci(1);
    pgclos();
    pgend();
    print "Press 'Q' to Quit or any other key to repeat";
    $command=<stdin>;
    chop($command);
    if ($command eq "N") {
	print "New Order of the Polynomical function?";
	$npoly=<stdin>;
	chop($npoly);
	$QUIT="N";
	$command="";
    }
    if ($command eq "Q") {
	$LOOP="E";
    } else {
	$QUIT="N";
	$command="";
    }
}

} while($LOOP eq "Q");
#
# We save the dispersion correction
#
print "Writting the dispersion solution at $disp_file\n";
open(FH,">$disp_file");
print FH "#New_index Wavelength OLD_index\n";
for ($i=0;$i<$n_cut;$i++) {
    $wave_n=$crval+$cdelt*$i;
#    $a_i_old[$i]=apr_n($a_i_old[$i],6);
    print FH "$i $wave_n $a_i_old[$i]\n";
}
close(FH);
print "DONE\n";

#exit;
#
# We correct for the new dispersion
#
my @out_array;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$spec[$i]=$in_array[$j][$i];
	$wave[$i]=$i;
    }  
    my $pdl_wave=pdl(@wave);
    my $pdl_spec=pdl(@spec);
    my $pdl_a_i_old=pdl(@a_i_old);
    my $INTER = init PDL::Func( Interpolate => "Hermite");
#    my $INTER = init PDL::Func( Interpolate => "Linear");
#
#    $INTER->set( x => $pdl_wave, y => $pdl_spec, bc => { start => [1,0], end => [1,0] } );              
    $INTER->set( x => $pdl_wave, y => $pdl_spec);              
#    $INTER->set( bc => { start => [5], end => [5] } );              
    $INTER->set( bc => { start => [3], end => [3] } );              
#    $INTER->set( bc => { start => [2,-1], end => [2,-1] } );              



#    my $out_spec_pdl = $INTER->interpolate($pdl_a_i_old);             

#    my  $out_spec_pdl=my_interpol(pdl(@a_i_old), pdl(@wave), pdl(@spec));
    my $out_spec_pdl = interpol(pdl(@a_i_old), pdl(@wave), pdl(@spec));
    @n_out=$out_spec_pdl->dims;
    $min=1e12;
    $max=-1e12;
    for ($i=0;$i<$n_cut;$i++) {
	$wave_new[$i]=$crval+$cdelt*$i;
	if ($i<$n_out[0]) {
	    $out_spec[$i]=$out_spec_pdl->at($i);
	} else {
	    $out_spec[$i]=$out_spec_pdl->at($n_out[0]-1);
	}
	$out_array[$j][$i]=$out_spec[$i];
	if ($min>$spec[$i]) {
	    $min=$spec[$i];
	}
	if ($max<$spec[$i]) {
	    $max=$spec[$i];
	}
    }

    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($crval,$crval+$cdelt*$n_cut,$min,$max,0,0);
	pglabel("Wavelength","Counts","$j/$ny");
	pgsci(1);
	pgline($n_cut,\@wave_new,\@out_spec);
	pgclos();
	pgend();
	print "Press Enter (A=Automatic)";
	$command=<stdin>;
	chop($command);
	if ($command eq "A") {
	    $plot=0;
	}
    }    
}


print "Writting the $out_file\n";
system("rm $out_file");
write_rss($out_file,$n_cut,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;


sub my_interpol() {
    my $pdl_a_i_old=$_[0];
    my $pdl_wave=$_[1];
    my $pdl_spec=$_[2];
    my $n,$i,$j;
    ($n)=$pdl_a_i_old->dims();;
    my $pdl_out=zeroes($n);
    for ($i=0;$i<$n;$i++) {
	print "PASO $i/$n\n";
	my $jmin=0;
	if ($pdl_a_i_old->at($i)<$pdl_wave->at(0)) {
	    set($pdl_out,$i,$pdl_spec->at(0));
	} else {
	    if ($pdl_a_i_old->at($i)>$pdl_wave->at($n-2)) {
		set($pdl_out,$i,$pdl_spec->at($n-2));
	    } else {
	        $j=0;
		while ($pdl_a_i_old->at($i)<$pdl_wave->at($j)) {
		    my $val1,$val2;
		    $val1=$pdl_a_i_old->at($i);
		    $val2=$pdl_wave->at($j);
		    print "$i $j $val1 $val2\n";
		    $j++;
		}
		my $val=$pdl_spec->at($j)+(($pdl_spec->at($j+1)-$pdl_spec->at($j))/($pdl_wave->at($j+1)-$pdl_wave->at($j)))*($pdl_a_i_old->at($i)-$pdl_wave->at($j));
	    }
	}
    }
#    print "NN = $n\n"; exit;

    return $pdl_out;
}
