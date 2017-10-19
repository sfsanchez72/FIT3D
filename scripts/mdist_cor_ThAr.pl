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


if ($#ARGV<9) {
    print "USE: mdist_cor_ppak.pl EXTRACTED_SCI.fits EXTRACTED_CAL.fits APERTURE NSIGMA NPOLY OUTPUT.fits DISTORSION.fits NX_min NX_max plot [APERTURE2] [NPOLY2] [NY_CEN] [INIT_SHIFT]\n";
    exit;
}

$infile=$ARGV[0];
$infile_cal=$ARGV[1];
#$dist_txt=$ARGV[1];
$aperture=$ARGV[2];
$nsigma=$ARGV[3];
$npoly=$ARGV[4];
$out_file=$ARGV[5];
$disp_file=$ARGV[6];
$nx_min=$ARGV[7];
$nx_max=$ARGV[8];
$plot=$ARGV[9];
$nsigma2=4;
$ny_cen=0;
$do_med=0;
$aperture2=$aperture;
if ($#ARGV==10) {
    $aperture2=$ARGV[10];
}

$npoly2=3;
if ($#ARGV==11) {
    $aperture2=$ARGV[10];
    $npoly2=$ARGV[11];
}

if ($#ARGV==12) {
    $aperture2=$ARGV[10];
    $npoly2=$ARGV[11];
    $ny_cen=$ARGV[12];
}
$init_shift=0;
if ($#ARGV==13) {
    $aperture2=$ARGV[10];
    $npoly2=$ARGV[11];
    $ny_cen=$ARGV[12];
    $init_shift=$ARGV[13];
}

print "NY_CEN=$ny_cen, $nbox,$nsigma2 $#ARGV\n";
print "Reading file $infile ";
$in_array=rfits($infile);
($nx,$ny)=$in_array->dims();
#$nax=read_naxes($infile);   
#@naxis=@$nax;
#$nx=$naxis[0];
#$ny=$naxis[1];
#
#$pdl=rfits($infile);
#@in_array=list($pdl);
#@in_array=list($pdl);
$in_array_cal=rfits($infile_cal);
#$in_array_cal=list($pdl);
#read_img($infile);
print "DONE\n";

#
$n=0;
$nc=15;
$ns=36;
#open(FH,"<ppak_positions.txt");
while($line=<DATA>) {
    chop($line);
    @data=split(" ",$line);
    $type[$n]=0;
    if (($data[1]>400)&&($data[1]<500))    {
	$type[$n]=1;
#	print "\n";
    } 

    if (($data[1]>500)&&($data[1]<900))    {
	$type[$n]=2;
#	print "\n";
    }

    if ($data[1]<400)    {
#	print "$data[0] ";
    }


    if ($data[1]<900) {
	$n++;
    }

    
}


#
# We construct the checking spectrum
#
if ($type[$ny_cen]!=2) {
    $ny_cen=0;
}
if ($nx_min<0) {
    $nx_min=0;
}
if ($nx_max>$nx) {
    $nx_max=$nx;
}

#print "$nx_min,$nx_max\n";
$min=1e12;
$max=-1e12;
my @spec;
for ($i=$nx_min;$i<$nx_max;$i++) {
    $w[$i]=$i*1.0;
    $w_cal[$i]=$w[$i]+$init_shift;
    $spec[$i]=$in_array->at($i,$ny_cen);
    $spec_cal[$i]=$in_array_cal->at($i,$ny_cen);
    if ($min>$spec[$i]) {
	$min=$spec[$i];
    }
    if ($max<$spec[$i]) {
	$max=$spec[$i];
    }
#    print "$i $ny_cen $spec[$i] $min $max\n";
}
$mean_flux=mean(@spec);
$sigma_flux=sigma(@spec);
print "Looking for peaks\n";
my @xp;
my @yp;
my @wp;
my @mp;
$np=0;
for ($i=$nx_min;$i<$nx_max;$i++) {
#    print "$spec[$i]-$mean_flux $nsigma*$sigma_flux\n";
    if (($spec[$i]-$mean_flux)>($nsigma*$sigma_flux)) {
	$xs=$i;
	$ys=$spec[int($xs)];
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
#	print "$xp[$np] $yp[$np]\n";
	$np++;
    }
}
print "DONE\n";
print "$np peaks found\n";

@borders_ini=(0,$nx,$min,$max);
@borders=(0,$nx,$min,$max);

#
# Main Loop to select the features.
#
$QUIT="Y";

$np=0;
open(FH,"ThAr.id");
while($line=<FH>) {
    @data=split(" ",$line);
    $xp[$np]=$data[0];
    $mp[$np]=1;
    $yp[$np]=$spec[int($data[0])];
    $np++;
}
close(FH);

    print "a => Unzoom\n";
    print "x => X-Zoom\n";
    print "y => Y-Zoom\n";
    print "m => Mark a peak\n";
    print "d => delete a peak\n";
    print "r => Move right\n";
    print "l => Move left\n";
    print "R => Load a file\n";
    print "W => Save a file\n";    
    print "Q => Quit\n";



while($QUIT eq "N") {
    pgbegin(0,"/null",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($borders[0],$borders[1],$borders[2],$borders[3],0,1);
    pgsci(1);
    pgsci(8);
    pgline($nx,\@w_cal,\@spec_cal);
    pgsci(7);
    pgline($nx,\@w,\@spec);
    pgsch(2);
    pgsci(5);
    for ($i=0;$i<$np;$i++) {
	if ($mp[$i]==1) {
	    pgpoint(1,[$xp[$i]],[$yp[$i]],2);
	}
    }
    pgsch(1.2);
    pgsci(1);
    $command="";
    while ($command ne "X") {
#	$x_cen=0.5*($borders[1]-$borders[0])+$borders[0];
#	$y_cen=0.5*($borders[3]-$borders[2])+$borders[2];
	pgband(7,0,50,50,$x,$y,$command);
#	print "$x $y $command\n";
#	pgsci(2);
#	pgpoint(1,[$x],[$y],2);
#	pgsci(1);
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
#		$wp[$np]=$data[1];
		$mp[$np]=1;
		$yp[$np]=$spec[int($data[0])];
		$np++;
	    }
	    close(FH);
	    $command="X";
	}
	if ($command eq "W") {
	    print "Saving Identification file?"; $file=<stdin>;
	    chop($file);
	    open(FH,">$file");
	    for ($kk=0;$kk<$np;$kk++) {
		if ($mp[$kk]==1) {
	#	    print FH "$xp[$kk] $wp[$kk]\n";
		    print FH "$xp[$kk]\n";
		}
	    }
	    close(FH);
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
	    $maximun=-1e12;
	    for ($i=$i_min;$i<$i_max;$i++) {
		if (($spec[$i]>$maximun)&&($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])) {
		    $i_peak=$i;
		    $maximun=$spec[$i];
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
	    $np++;
	    pgsch(2);
	    pgsci(5);
	    for ($i=0;$i<$np;$i++) {
		if ($mp[$i]==1) {
		    pgpoint(1,[$xp[$i]],[$yp[$i]],2);
		}
	    }
	    pgsch(1.2);

#	    print "Wavelength at $x_peak?"; $data=<stdin>;
#	    chop($data);
#	    $wp[$np]=$data;

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
	    pgsch(2);
	    for ($i=0;$i<$np;$i++) {
		if ($mp[$i]==1) {
		    pgsci(5);
		    pgpoint(1,[$xp[$i]],[$yp[$i]],2);
		} else {
		    pgsci(2);
		    pgpoint(1,[$xp[$i]],[$yp[$i]],22);
		}
	    }
	    pgsci(1);
	    pgsch(1.2);
	}


    }

    pgclos();
    pgend();

}
print "DONE\n";
$command="A";
#
# Cleaning the peaks
#

$n_peak=0;
for ($i=0;$i<$np;$i++) {
    if ($mp[$i]==1) {
	$peaks_cen[$n_peak]=$xp[$i];
	$peaks_now[$n_peak]=$xp[$i];
	$peaks_now_cal[$n_peak]=$xp[$i];#-$init_shift;
	$peaks_last[$n_peak]=$xp[$i];
	$n_peak++;
    }
}
@peaks_cen = sort {$a <=> $b} @peaks_cen;
@peaks_now = sort {$a <=> $b} @peaks_now;
@peaks_last = sort {$a <=> $b} @peaks_last;
print "Final number of peaks=$n_peak\n";
#@peaks_now=@peaks_cen;
#
# We determine the offset with respect to the reference spectrum:
#
#@borders=@borders_ini;
$ny_cen=int($ny/2);
$min=1e12;
$max=-1e12;
my @spec;
for ($j=0;$j<$ny;$j++) {
#    print "$j $type[$j]\n";
    if ($type[$j]==2) {
	for ($i=0;$i<$nx;$i++) {
	    $w[$i]=$i*1.0;
	    $w_cal[$i]=$w[$i]+$init_shift;
	    $spec[$i]=$in_array->at($i,$j);
	    $spec_cal[$i]=$in_array_cal->at($i,$j);
	    if ($min>$spec[$i]) {
		$min=$spec[$i];
	    }
	    if ($max<$spec[$i]) {
		$max=$spec[$i];
	    }
	}
	$mean_flux=median(@spec);
	$sigma_flux=sigma(@spec);
#
# We look for the peaks in the SCIENCE frame
#
	for ($k=0;$k<$n_peak;$k++) {	
	    $xs=$peaks_last[$k];#+$shift_all[$j];
	    $i_min=int($xs)-$aperture;
	    $i_max=int($xs)+$aperture;
	    $i_peak=int($xs);
	    $maximun=-1e12;
	    $peak_found=0;
	    for ($i=$i_min;$i<$i_max;$i++) {
		if ($spec[$i]>$maximun) {
		    $i_peak=$i;
		    $maximun=$spec[$i];
		    $peak_found=1;
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
	    if ($peak_found==0) {
		$x_peak=$xs;
	    }
	    $peaks_now[$k]=$x_peak;
	    $dpeaks_now[$k]=$peaks_cen[$k]-$x_peak;
	    $a_peaks[$k][$j]=$x_peak;
    }

#
# We look for the peaks in the Calibration frame
#
	for ($k=0;$k<$n_peak;$k++) {	
	    $xs=$peaks_now[$k]+$init_shift;#+$shift_all[$j];
#	    $i_min=int($xs)-$aperture;
#	    $i_max=int($xs)+$aperture;
	    $i_min=int($xs)-$aperture2;
	    $i_max=int($xs)+$aperture2;
	    $i_peak=int($xs);
	    $maximun=-1e12;
	    $peak_found=0;
	    for ($i=$i_min;$i<$i_max;$i++) {
		if ($spec_cal[$i]>$maximun) {
		    $i_peak=$i;
		    $maximun=$spec_cal[$i];
		    $peak_found=1;
		}
	    }
	    $a=$i_peak-1;
	    $b=$i_peak;
	    $c=$i_peak+1;
	    $fa=-$spec_cal[$a];
	    $fb=-$spec_cal[$b];
	    $fc=-$spec_cal[$c];
	    $den=($fc-2*$fb+$fa);
	    if ($den!=0) {
		$x_peak=$c-($b-$a)*(($fc-$fb)/$den+0.5);
		
	    } else {
		$x_peak=$i_peak;
	    }
	    if ($peak_found==0) {
		$x_peak=$xs;
	    }
	    $peaks_now_cal[$k]=$x_peak;
	    $dpeaks_now_cal[$k]=$peaks_cen[$k]-$x_peak;
	    $a_peaks_cal[$k][$j]=$x_peak;
    }

#    $command="";

#    $command="";
	if (($plot==1)||($plot==2)) {
	    pgbegin(0,"/null",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv($borders[0],$borders[1],$borders[2],$borders[3],0,1);
	    pglabel("X","counts","$j/$ny");
	    pgsci(1);
	    pgsci(8);
	    pgline($nx,\@w_cal,\@spec_cal);
	    pgsci(7);
	    pgline($nx,\@w,\@spec);
	    pgsch(2);
	    for ($i=0;$i<$n_peak;$i++) {
		$xx=$peaks_last[$i];#+$shift_all[$j];
		$yy=$spec[int($peaks_last[$i])];
		pgsci(4);
		pgpoint(1,[$xx],[$yy],2);
		$xx=$peaks_now[$i];
		$yy=$spec[int($peaks_now[$i])];
		pgsci(3);
		pgslw(3);
		pgline(2,[$xx,$xx],[1.2*$yy,1.6*$yy]);
		pgslw(4);
		$xx=$peaks_now_cal[$i];
		$yy=$spec[int($peaks_now_cal[$i])];
		pgsci(2);
		pgslw(3);
		pgline(2,[$xx,$xx],[1.2*$yy,1.6*$yy]);
		pgslw(4);
	    }
	    pgsch(1.2);
	    pgsci(1);
	    pgclos();
	    pgend();
	    if (($command ne "A")&&($command ne "0")) {
		print "Press Enter (A=Automatic, 0=No plot, N=Next plot)";
		$command=<stdin>;
		chop($command);
		if ($command eq "Q") {
		    $plot=0;
		}
		if ($command eq "0") {
		    $plot=0;
		}
		if ($command eq "N") {
		    $plot++;
		}
	    }
	    
	} else {
	    if ($j==100*int($j/100)) {
		print "$j/$ny\n";
	    }
#	print "$j/$ny\n";
	}
	for ($k=0;$k<$n_peak;$k++) {	
	    $peaks_last[$k]=$peaks_now[$k];
	}
    }
}
#
#
# Smoothing process
#
#
#$nbox=5;
#$nsigma2=2;


    $command="";
$command="A";

if (($plot==1)||($plot==2)||($plot==3)) {

    pgbegin(0,"/null",1,1);
    pgask(0);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    for ($k=0;$k<$n_peak;$k++) {	
	$kk=0;
	for ($j=0;$j<$ny;$j++) {	
	    if ($type[$j]==2) {
		$xp[$kk]=$a_peaks[$k][$j];		
		$xp_cal[$kk]=$a_peaks_cal[$k][$j];		
		$delta_xp[$kk]=$a_peaks[$k][$j]-$a_peaks_cal[$k][$j];
		$yp[$kk]=$j;		
		$kk++;
	    }
	    $yyp[$j]=$j;
	}
#	print "PASO\n";
	($s_f,$coeff) = fitpoly1d(pdl(@yp),pdl(@xp),$npoly2);
	for ($j=0;$j<$ny;$j++) {	
	    $a_new_peaks[$k][$j]=0;
	    for ($kp=0;$kp<$npoly2;$kp++) {
		$c=$coeff->at($kp);
		$a_new_peaks[$k][$j]=$a_new_peaks[$k][$j]+$c*($j**$kp);
	    }	    
	    $sxp[$j]=$a_new_peaks[$k][$j];
	}	
#	print "PASO 1\n";
	($s_f,$coeff) = fitpoly1d(pdl(@yp),pdl(@xp_cal),$npoly2);
#	print "PASO 1\n";
	for ($j=0;$j<$ny;$j++) {	
	    $a_new_peaks_cal[$k][$j]=0;
	    for ($kp=0;$kp<$npoly2;$kp++) {
		$c=$coeff->at($kp);
		$a_new_peaks_cal[$k][$j]=$a_new_peaks_cal[$k][$j]+$c*($j**$kp);
	    }	    
	    $sxp_cal[$j]=$a_new_peaks_cal[$k][$j];
	}	
#	print "PASO 2\n";
	($s_f,$coeff) = fitpoly1d(pdl(@yp),pdl(@delta_xp),$npoly2);
	for ($j=0;$j<$ny;$j++) {	
	    $a_delta_peaks[$k][$j]=0;
	    for ($kp=0;$kp<$npoly2;$kp++) {
		$c=$coeff->at($kp);
		$a_delta_peaks[$k][$j]=$a_delta_peaks[$k][$j]+$c*($j**$kp);
	    }	    
	    $delta_sxp[$j]=$a_delta_peaks[$k][$j];
	}	
#	print "PASO 3\n";

	$nx_med=median(@xp);
	$sig_med=sigma(@xp);
	pgenv($nx_med-2*$sig_med,$nx_med+2*$sig_med,0,$ny,0,0);
	pglabel("Y-Axis","X-axis","Distortion");
	pgsci(1);
	pgpoint($kk,\@xp,\@yp,2);
	pgsci(2);
	pgpoint($kk,\@xp_cal,\@yp,2);
	pgsci(3);
	pgline($ny,\@sxp,\@yyp);
	pgsci(8);
	pgline($ny,\@sxp_cal,\@yyp);
	pgsci(1);
	if (($command ne "A")&&($command ne "0")) {
#	if ($command ne "A") {
	    print "Press Enter (A=Automatic, 0=No plot, N=Next plot)";
#	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }	
	    if ($command eq "N") {
		$plot++;
	    }
	}
	$sig_med=sigma(@delta_xp);
	pgenv(-2,+2,0,$ny,0,0);
	pglabel("Y-Axis","Delta X","Distortion");
	pgsci(1);
	pgpoint($kk,\@delta_xp,\@yp,2);
	pgsci(3);
	pgline($ny,\@delta_sxp,\@yyp);
	pgsci(1);
	if (($command ne "A")&&($command ne "0")) {
#	if ($command ne "A") {
	    print "Press Enter (A=Automatic, 0=No plot, N=Next plot)";
#	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }	
	    if ($command eq "N") {
		$plot++;
	    }
	}
    }
    pgclos();
    pgend();
}

$command="A";
$plot=0;

$dnx_max=-1e12;
$dnx_min=1e12;
print "Fitting the distortion\n";
for ($j=0;$j<$ny;$j++) {
    for ($k=0;$k<$n_peak;$k++) {
	$peaks_now[$k]=$a_new_peaks[$k][$j];	
	$peaks_cal[$k]=$a_new_peaks_cal[$k][$j];	
#	$dpeaks_now[$k]=$peaks_cen[$k]-$a_new_peaks[$k][$j];
	$dpeaks_now[$k]=$a_new_peaks_cal[$k][$j]-$a_new_peaks[$k][$j];
#	print "($j/$ny) ($k/$n_peak) $peaks_now[$j] $peaks_cen[$j] \n";
	if ($dnx_max<$dpeaks_now[$k]) {
	    $dnx_max=$dpeaks_now[$k]+0.5;#*$dpeaks_now[$k];
	}
	if ($dnx_min>$dpeaks_now[$k]) {
	    $dnx_min=$dpeaks_now[$k]-0.5;#*$dpeaks_now[$k];
	}
    }
    if ($npoly>0) {
	($s_f,$coeff) = fitpoly1d(pdl(@peaks_cal),pdl(@peaks_now),$npoly);
	for ($i=0;$i<$nx;$i++) {
	    $a_dist[$j][$i]=0;
	    for ($kk=0;$kk<$npoly;$kk++) {
		$c=$coeff->at($kk);
		$a_dist[$j][$i]=$a_dist[$j][$i]+$c*($i**$kk);
	    }	    
	    $a_i_new[$i]=$a_dist[$j][$i];
	    $a_i_old[$i]=$i;
	    $da_i_new[$i]=-$a_i_new[$i]+$a_i_old[$i];
	}	
    } else {
	my $spline=new Math::Spline(\@peaks_cal,\@peaks_now);
	for ($i=0;$i<$nx;$i++) {
	    $a_dist[$j][$i]=$spline->evaluate($i);
	    $a_i_new[$i]=$a_dist[$j][$i];
	    $a_i_old[$i]=$i;
	    $da_i_new[$i]=-$a_i_new[$i]+$a_i_old[$i];

	}
    }
    if (($plot==1)||($plot==2)||($plot==3)||($plot==4)) {
	pgbegin(0,"/null",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
#	pgenv(0,$nx,0,$nx,1,0);
#	pglabel("X old","X new","Distorsion solution ($j/$ny)");
#	pgsci(1);
#	pgpoint($n_peak,\@peaks_now,\@peaks_cen,2);
#	pgsci(2);
#	pgline($nx,\@a_i_old,\@a_i_new);
	pgenv(0,$nx,$dnx_min,$dnx_max,0,0);
	pglabel("X old","\\gD X new","Distorsion solution ($j/$ny)");
	pgsci(1);
	pgpoint($n_peak,\@peaks_now,\@dpeaks_now,2);
	pgsci(2);
	pgline($nx,\@a_i_old,\@da_i_new);
	pgclos();
	pgend();
	if (($command ne "A")&&($command ne "0")) {
#	if ($command ne "A") {
	    print "Press Enter (A=Automatic, 0=No plot, N=Next plot)";
#	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	    if ($command eq "N") {
		$plot++;
	    }
	}
    }  else {
	if ($j==100*int($j/100)) {
	    print "$j/$ny\n";
	}
    }
#    print ".";
}
print "DONE\n";
    $command="";
$command="0";
$plot=0;
#
# We save the dispersion correction
#
print "Writting the distorsion solution at $disp_file\n";
system("rm $disp_file");
write_img($disp_file,$nx,$ny,\@a_dist);
print "DONE\n";


#
# We correct for the new dispersion
#
my @out_array;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$a_i_old[$i]=$a_dist[$j][$i];
	$spec[$i]=$in_array->at($i,$j);
	$wave[$i]=$i;
    }
    my $out_spec_pdl = interpol(pdl(@a_i_old), pdl(@wave), pdl(@spec));
    $min=1e12;
    $max=-1e12;
    for ($i=0;$i<$nx;$i++) {
	$wave_new[$i]=$i;
	$out_spec[$i]=$out_spec_pdl->at($i);
	$out_array[$j][$i]=$out_spec[$i];
	if ($min>$spec[$i]) {
	    $min=$spec[$i];
	}
	if ($max<$spec[$i]) {
	    $max=$spec[$i];
	}
    }
    if (($plot==1)||($plot==2)) {
	pgbegin(0,"/null",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$nx,$min,$max,0,0);
	pglabel("Wavelength","Counts","$j/$ny");
	pgsci(2);
	pgline($nx,\@wave_new,\@spec);
	pgsci(1);
	pgline($nx,\@wave_new,\@out_spec);
	pgclos();
	pgend();
	if (($command ne "A")&&($command ne "0")) {
#	if ($command ne "A") {
	    print "Press Enter (A=Automatic, 0=No plot, N=Next plot)";
#	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;	# 
	    }
	    if ($command eq "N") {
		$plot++;
	    }
	}
    }      else {
	if ($j==100*int($j/100)) {
	    print "$j/$ny\n";
	}
    }
}
    $command="";


print "Writting the $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$ny,\@out_array);
#write_rss($out_file,$n_cut,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;

if ($nbox>0) {
    for ($k=0;$k<$n_peak;$k++) {	
	my @xp;
	for ($j=0;$j<$ny;$j++) {
	    my @tmp;
	    $kk=0;
	    if ($j<$nbox) {
		$j_min=$j;
	    $j_max=$j+2*$nbox;
	    } else {
		if ($j>=($ny-$nbox)) {
		    $j_min=$j-2*$nbox;
		    $j_max=$j;
		} else {
		    $j_min=$j-$nbox;
		    $j_max=$j+$nbox;
		}
	    }
	    for ($i=$j_min;$i<$j_max;$i++) {
		$tmp[$kk]=$a_peaks[$k][$i];
		$kk++;
	    }
	    $me=median(@tmp);
	    $sig=sigma(@tmp);
	    if (abs($xp[$j]-$me)<$nsigma2*$sig) {
		$xp[$kk]=$a_peaks[$k][$j];
		$kk++;
	    } else {
		$xp[$j]=$me;
	    }
	}
	for ($j=0;$j<$ny;$j++) {	
	    $a_new_peaks[$k][$j]=$xp[$j];
	}
    }
} else {
    @a_new_peaks=@a_peaks;
}




__DATA__
      1     501    501       501    
      2     502    502       502    
      3     999    999       999    
      4       1      0.780     0.000
      5       2     -0.878    -0.169
      6     401    401       401    
      7       3      0.000    -0.676
      8       4      0.780    -0.676
      9       5     -0.780    -0.676
     10       6      0.585    -0.338
     11       7     -0.390    -1.013
     12       8      0.390    -1.013
     13       9      0.975    -0.338
     14      10     -0.488    -0.507
     15      11      0.390    -0.676
     16      12      0.000    -1.013
     17     402    402       402    
     18      13     -0.975    -0.338
     19      14      1.267    -0.169
     20      15     -0.390    -0.676
     21      16      0.683    -0.844
     22      17     -0.683    -0.169
     23      18      0.098    -0.844
     24      19     -0.683    -0.844
     25      20      1.170     0.000
     26      21     -1.267    -0.169
     27      22      0.683    -0.507
     28     403    403       403    
     29      23     -0.683    -1.182
     30      24      0.293    -0.844
     31      25      0.878    -0.844
     32      26     -0.683    -0.507
     33     503    503       503    
     34      27      0.878    -0.169
     35      28     -0.293    -0.844
     36      29      0.098    -1.182
     37      30     -1.072    -0.507
     38     404    404       404    
     39      31      1.072    -0.507
     40      32     -0.585    -0.338
     41      33     -0.293    -1.182
     42      34      0.195    -0.676
     43      35      0.488    -1.182
     44      36     -0.195    -0.676
     45      37     -0.878    -0.844
     46      38      0.683    -0.169
     47      39      0.780    -1.013
     48      40     -0.585    -0.676
     49     405    405       405    
     50      41      1.072    -0.169
     51      42     -0.780    -0.338
     52      43      0.585    -0.676
     53      44     -0.195    -1.013
     54      45     -0.975    -0.676
     55      46      0.683    -1.182
     56      47      1.170    -0.338
     57      48     -0.585    -1.013
     58      49      0.488    -0.844
     59      50     -1.072    -0.169
     60     406    406       406    
     61      51     -0.098    -0.844
     62      52      0.293    -1.182
     63      53      0.878    -0.507
     64      54     -0.780    -1.013
     65     504    504       504    
     66      55     -0.878    -0.507
     67      56      0.975     0.000
     68      57     -0.098    -1.182
     69      58      0.585    -1.013
     70     407    407       407    
     71      59     -1.170    -0.338
     72      60     -0.488    -0.844
     73      61      0.488    -0.507
     74      62      1.365     0.000
     75      63      0.195    -1.013
     76      64      0.975    -0.676
     77      65     -0.488    -1.182
     78      66      0.780    -0.338
     79      67     -1.365    -0.338
     80      68      0.390    -1.351
     81     408    408       408    
     82      69      1.365    -0.338
     83      70     -0.975    -1.013
     84      71     -0.488    -1.520
     85      72      1.072    -0.844
     86      73      1.658    -0.169
     87      74     -1.365    -0.676
     88      75      0.000    -1.351
     89      76      0.780    -1.351
     90      77     -0.975    -1.351
     91      78      1.365    -0.676
     92     409    409       409    
     93      79      0.293    -1.520
     94      80     -1.170    -0.676
     95      81      1.560     0.000
     96      82     -0.390    -1.351
     97     505    505       505    
     98      83     -1.755    -0.338
     99      84      1.072    -1.182
    100      85      1.560    -0.338
    101      86     -1.170    -1.013
    102     410    410       410    
    103      87     -0.195    -1.689
    104      88      0.683    -1.520
    105      89      1.267    -0.844
    106      90     -0.780    -1.351
    107      91     -1.462    -0.507
    108      92      1.950     0.000
    109      93     -1.365    -1.013
    110      94      1.072    -1.520
    111      95     -0.195    -1.351
    112      96     -1.658    -0.169
    113     411    411       411    
    114      97     -0.975    -1.689
    115      98      1.267    -0.507
    116      99      0.098    -1.520
    117     100      1.365    -1.013
    118     101     -1.560    -0.676
    119     102      0.878    -1.182
    120     103     -1.072    -1.182
    121     104      0.780    -1.689
    122     105     -0.683    -1.520
    123     106     -1.462    -0.169
    124     412    412       412    
    125     107      1.755     0.000
    126     108      0.390    -1.689
    127     109      1.658    -0.507
    128     110     -1.267    -0.844
    129     506    506       506    
    130     111     -1.170    -1.351
    131     112      0.975    -1.013
    132     113      0.195    -1.351
    133     114     -0.585    -1.689
    134     413    413       413    
    135     115      1.462    -0.169
    136     116     -1.852    -0.169
    137     117     -0.098    -1.520
    138     118      0.975    -1.351
    139     119     -1.267    -0.507
    140     120      1.462    -0.844
    141     121      0.488    -1.520
    142     122     -1.072    -0.844
    143     123     -0.878    -1.520
    144     124      1.170    -0.676
    145     414    414       414    
    146     125      1.755    -0.338
    147     126      1.170    -1.013
    148     127     -1.560    -0.338
    149     128      0.878    -1.520
    150     129     -1.267    -1.182
    151     130     -0.293    -1.520
    152     131      1.462    -0.507
    153     132      1.170    -1.351
    154     133      0.195    -1.689
    155     134     -1.658    -0.507
    156     415    415       415    
    157     135      1.852    -0.169
    158     136      0.585    -1.351
    159     137     -1.072    -1.520
    160     138      1.560    -0.676
    161     507    507       507    
    162     139     -0.390    -1.689
    163     140      1.267    -1.182
    164     141      0.975    -1.689
    165     142     -0.585    -1.351
    166     416    416       416    
    167     143     -1.462    -0.844
    168     144      0.000    -1.689
    169     145     -0.878    -1.182
    170     146     -0.780    -1.689
    171     147      0.585    -1.689
    172     148     -0.585     0.000
    173     149     -0.098    -0.507
    174     150      0.488    -0.169
    175     151     -0.293     0.507
    176     152      0.195     0.338
    177     417    417       417    
    178     153     -0.293    -0.169
    179     154      0.195    -0.338
    180     155      0.488     0.169
    181     156     -0.293    -0.507
    182     157     -0.488     0.169
    183     158      0.293    -0.507
    184     159      0.098     0.507
    185     160     -0.195    -0.338
    186     161      0.585     0.000
    187     162      0.390     0.338
    188     418    418       418    
    189     163     -0.195     0.000
    190     164      0.000     0.000
    191     165      0.195     0.000
    192     508    508       508    
    193     509    509       509    
    194     166     -0.098     0.169
    195     167     -0.098    -0.169
    196     168      0.098     0.169
    197     169      0.098    -0.169
    198     419    419       419    
    199     170     -0.390     0.338
    200     171     -0.390    -0.338
    201     172      0.293     0.507
    202     173      0.293    -0.169
    203     174     -0.390     0.000
    204     175      0.098    -0.507
    205     176      0.000     0.338
    206     177      0.293     0.169
    207     178     -0.293     0.169
    208     179      0.000    -0.338
    209     420    420       420    
    210     180     -0.098     0.507
    211     181      0.390    -0.338
    212     182     -0.488    -0.169
    213     183      0.390     0.000
    214     184     -0.195     0.338
    215     185     -1.950     0.000
    216     186      0.585     1.689
    217     187     -1.170     1.351
    218     188      1.267     0.844
    219     189     -1.365     0.676
    220     421    421       421    
    221     190     -0.195     1.689
    222     191      1.170     1.351
    223     192     -1.658     0.169
    224     193      1.755     0.338
    225     510    510       510    
    226     194     -0.780     1.689
    227     195      0.293     1.520
    228     196     -1.170     1.013
    229     197      0.878     1.520
    230     422    422       422    
    231     198     -1.658     0.507
    232     199      1.462     0.844
    233     200     -1.072     1.520
    234     201     -0.293     1.520
    235     202      1.852     0.169
    236     203     -1.462     0.844
    237     204      1.462     0.507
    238     205     -1.852     0.169
    239     206      0.195     1.689
    240     207      1.072     1.182
    241     423    423       423    
    242     208      0.975     1.689
    243     209     -0.585     1.689
    244     210     -1.072     1.182
    245     211      1.658     0.169
    246     212      1.365     1.013
    247     213     -1.560     0.676
    248     214      1.560     0.676
    249     215      0.975     1.351
    250     216     -1.755     0.000
    251     217      0.000     1.689
    252     424    424       424    
    253     218     -0.780     1.351
    254     219     -1.365     0.338
    255     220      0.585     1.351
    256     221     -1.267     1.182
    257     511    511       511    
    258     222      1.170     0.676
    259     223     -0.195     1.351
    260     224      1.267     1.182
    261     225     -1.560     0.000
    262     425    425       425    
    263     226      0.780     1.689
    264     227      1.462     0.169
    265     228     -1.072     0.844
    266     229     -0.683     1.520
    267     230     -1.755     0.338
    268     231      0.975     1.013
    269     232      0.390     1.689
    270     233      1.658     0.507
    271     234      0.195     1.351
    272     235     -1.365     1.013
    273     426    426       426    
    274     236     -0.975     1.689
    275     237      1.072     1.520
    276     238     -0.390     1.689
    277     239      1.365     0.338
    278     240     -1.462     0.507
    279     241      0.683     1.520
    280     242     -0.975     1.351
    281     243      1.365     0.676
    282     244     -0.390     1.351
    283     245     -1.462     0.169
    284     427    427       427    
    285     246      0.098     1.520
    286     247     -1.170     0.676
    287     248      0.878     1.182
    288     249      1.560     0.338
    289     512    512       512    
    290     250     -0.878     1.182
    291     251      0.390     1.351
    292     252      1.072     0.844
    293     253     -1.560     0.338
    294     428    428       428    
    295     254     -0.488     1.520
    296     255      0.780     1.351
    297     256     -1.267     0.844
    298     257      0.000     1.351
    299     258      1.267     0.507
    300     259     -0.878     1.520
    301     260     -1.267     0.507
    302     261      1.170     1.013
    303     262      0.488     1.520
    304     263     -0.585     1.351
    305     429    429       429    
    306     264     -0.975     1.013
    307     265     -0.098     1.520
    308     266     -1.365     0.000
    309     267      1.072     0.507
    310     268      0.488     1.182
    311     269     -0.293     1.182
    312     270     -0.975     0.676
    313     271      0.878     0.844
    314     272     -1.072     0.169
    315     273      1.267     0.169
    316     430    430       430    
    317     274      0.195     1.013
    318     275     -0.780     1.013
    319     276      0.683     1.182
    320     277     -0.390     1.013
    321     513    513       513    
    322     278      0.878     0.507
    323     279     -1.072     0.507
    324     280      0.098     1.182
    325     281      0.488     0.844
    326     431    431       431    
    327     282     -1.170     0.000
    328     283      1.072     0.169
    329     284     -0.683     1.182
    330     285      0.975     0.676
    331     286     -0.585     0.676
    332     287     -0.098     0.844
    333     288     -0.878     0.169
    334     289      0.780     1.013
    335     290      0.780     0.338
    336     291     -0.878     0.844
    337     432    432       432    
    338     292      0.293     1.182
    339     293     -1.267     0.169
    340     294     -0.488     1.182
    341     295      0.780     0.676
    342     296     -0.878     0.507
    343     297      0.195     0.676
    344     298     -0.098     1.182
    345     299      1.170     0.338
    346     300      0.585     1.013
    347     301     -0.975     0.000
    348     433    433       433    
    349     302     -0.488     0.844
    350     303     -1.170     0.338
    351     304      0.683     0.507
    352     305     -0.585     0.338
    353     514    514       514    
    354     306      0.000     1.013
    355     307      0.683     0.844
    356     308     -0.683     0.844
    357     309      0.390     1.013
    358     434    434       434    
    359     310     -0.195     0.676
    360     311     -0.975     0.338
    361     312      0.585     0.338
    362     313     -0.780     0.676
    363     314      0.390     0.676
    364     315     -0.585     1.013
    365     316     -0.195     1.013
    366     317     -0.780     0.338
    367     318      0.098     0.844
    368     319      0.975     0.338
    369     435    435       435    
    370     320     -0.390     0.676
    371     321      0.488     0.507
    372     322     -0.488     0.507
    373     323      0.293     0.844
    374     324      0.878     0.169
    375     325     -0.780     0.000
    376     326     -0.293     0.844
    377     327      0.000     0.676
    378     328      0.585     0.676
    379     329     -0.683     0.507
    380     436    436       436    
    381     330     -0.683     0.169
    382     331      0.683     0.169
    383     999    999       999    
    384     515    515       515    
