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


if ($#ARGV<8) {
    print "USE: mdist_cor_sp.pl EXTRACTED.fits APERTURE NSIGMA NPOLY OUTPUT.fits DISTORSION.fits NX_min NX_max plot [NSIGMA2] [BOX] [NY_CEN]\n";
    exit;
}

$infile=$ARGV[0];
#$infile_cal=$ARGV[0];
#$dist_txt=$ARGV[1];
$aperture=$ARGV[1];
$nsigma=$ARGV[2];
$npoly=$ARGV[3];
$out_file=$ARGV[4];
$disp_file=$ARGV[5];
$nx_min=$ARGV[6];
$nx_max=$ARGV[7];
$plot=$ARGV[8];
$nsigma2=4;
$ny_cen=0;
$do_med=0;
$nsigma2=3;
if ($#ARGV==9) {
    $nsigma2=$ARGV[9];
}
$nbox=0;
if ($#ARGV==10) {
    $nsigma2=$ARGV[9];
    $nbox=$ARGV[10];
}
if ($#ARGV==11) {
    $nsigma2=$ARGV[9];
    $nbox=$ARGV[10];
    $ny_cen=$ARGV[11];
}

print "NY_CEN=$ny_cen, $nbox,$nsigma2 $#ARGV\n";
print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
#
#$pdl=rfits($infile);
#@in_array=list($pdl);
@in_array=read_img($infile);
print "DONE\n";





#
# We construct the checking spectrum
#
print "$nx_min,$nx_max\n";
$min=1e12;
$max=-1e12;
my @spec;
for ($i=$nx_min;$i<$nx_max;$i++) {
    $w[$i]=$i*1.0;
    $spec[$i]=$in_array[$ny_cen][$i];
    if ($min>$spec[$i]) {
	$min=$spec[$i];
    }
    if ($max<$spec[$i]) {
	$max=$spec[$i];
    }
#    print "$i $spec[$i] $min $max\n";
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
$QUIT="N";

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
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($borders[0],$borders[1],$borders[2],$borders[3],0,1);
    pgsci(1);
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
#
# Cleaning the peaks
#

$n_peak=0;
for ($i=0;$i<$np;$i++) {
    if ($mp[$i]==1) {
	$peaks_cen[$n_peak]=$xp[$i];
	$peaks_now[$n_peak]=$xp[$i];
	$n_peak++;
    }
}
@peaks_cen = sort {$a <=> $b} @peaks_cen;
@peaks_now = sort {$a <=> $b} @peaks_now;
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
    for ($i=0;$i<$nx;$i++) {
	$w[$i]=$i*1.0;
	$spec[$i]=$in_array[$j][$i];
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
# We look for the new peaks
#
#    my @peaks_now;
    for ($k=0;$k<$n_peak;$k++) {	
	$xs=$peaks_cen[$k];#+$shift_all[$j];
	#print "$k $xs\n";
	$i_min=int($xs)-$aperture;
	$i_max=int($xs)+$aperture;
	$i_peak=int($xs);
	$maximun=-1e12;
	$peak_found=0;
	for ($i=$i_min;$i<$i_max;$i++) {
#		if(($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])) {
#	    if (($spec[$i]>$maximun)&&($spec[$i]>$spec[$i-1])&&($spec[$i]>$spec[$i+1])) {
	    if ($spec[$i]>$maximun) {
		$i_peak=$i;
		$maximun=$spec[$i];
		$peak_found=1;
	    }
	}
#	    }
#	}
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
#	print "$j,$k $a_peaks[$k][$j]\n";

#	print "$k $peaks_cen[$k] $peaks_now[$k]\n";

    }

#    $command="";
    if (($plot==1)||($plot==2)) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($borders[0],$borders[1],$borders[2],$borders[3],0,1);
	pglabel("X","counts","$j/$ny");
	pgsci(1);
	pgsci(7);
	pgline($nx,\@w,\@spec);
	pgsch(2);
	for ($i=0;$i<$n_peak;$i++) {
	    $xx=$peaks_cen[$i];#+$shift_all[$j];
#	    $yy=$spec[int($peaks_cen[$i]+$shift_all[$j])];
	    $yy=$spec[int($peaks_cen[$i])];
	    pgsci(4);
	    pgpoint(1,[$xx],[$yy],2);
	    $xx=$peaks_now[$i];
	    $yy=$spec[int($peaks_now[$i])];
	    pgsci(3);
	    pgslw(3);
	    pgline(2,[$xx,$xx],[1.2*$yy,1.6*$yy]);
	    pgslw(4);
#	    print "$i $peaks_cen[$i] $peaks_now[$i]\n";
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
}

#
#
# Smoothing process
#
#
#$nbox=5;
#$nsigma2=2;
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

    $command="";
if (($plot==1)||($plot==2)||($plot==3)) {

    pgbegin(0,"/xs",1,1);
    pgask(0);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    for ($k=0;$k<$n_peak;$k++) {	
	for ($j=0;$j<$ny;$j++) {	
	    $xp[$j]=$a_peaks[$k][$j];
	    $sxp[$j]=$a_new_peaks[$k][$j];
	    $yp[$j]=$j;
	}
	$nx_med=median(@xp);
	$sig_med=sigma(@xp);
	pgenv($nx_med-2*$sig_med,$nx_med+2*$sig_med,0,$ny,0,0);
	pglabel("Y-Axis","X-axis","Distortion");
	pgsci(1);
	pgpoint($ny,\@xp,\@yp,2);
	pgsci(3);
	pgline($ny,\@sxp,\@yp);
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

$dnx_max=-1e12;
$dnx_min=1e12;
print "Fitting the distortion\n";
for ($j=0;$j<$ny;$j++) {
    for ($k=0;$k<$n_peak;$k++) {
	$peaks_now[$k]=$a_new_peaks[$k][$j];	
	$dpeaks_now[$k]=$peaks_cen[$k]-$a_new_peaks[$k][$j];
#	print "($j/$ny) ($k/$n_peak) $peaks_now[$j] $peaks_cen[$j] \n";
	if ($dnx_max<$dpeaks_now[$k]) {
	    $dnx_max=$dpeaks_now[$k]+0.5;#*$dpeaks_now[$k];
	}
	if ($dnx_min>$dpeaks_now[$k]) {
	    $dnx_min=$dpeaks_now[$k]-0.5;#*$dpeaks_now[$k];
	}
    }
    if ($npoly>0) {
	($s_f,$coeff) = fitpoly1d(pdl(@peaks_cen),pdl(@peaks_now),$npoly);
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
	my $spline=new Math::Spline(\@peaks_cen,\@peaks_now);
	for ($i=0;$i<$nx;$i++) {
	    $a_dist[$j][$i]=$spline->evaluate($i);
	    $a_i_new[$i]=$a_dist[$j][$i];
	    $a_i_old[$i]=$i;
	    $da_i_new[$i]=-$a_i_new[$i]+$a_i_old[$i];

	}
    }
    if (($plot==1)||($plot==2)||($plot==3)||($plot==4)) {
	pgbegin(0,"/xs",1,1);
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
	$spec[$i]=$in_array[$j][$i];
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
	pgbegin(0,"/xs",1,1);
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
