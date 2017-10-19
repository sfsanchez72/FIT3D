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
    print "USE: mdist_cor.pl EXTRACTED.fits APERTURE NSIGMA NPOLY OUTPUT.fits DISTORSION.fits NX_min NX_max plot [NPOLY2] [NY_CEN] [MEDIAN=0/1]\n";
    exit;
}

$infile=$ARGV[0];
#$dist_txt=$ARGV[1];
$aperture=$ARGV[1];
$nsigma=$ARGV[2];
$npoly=$ARGV[3];
$out_file=$ARGV[4];
$disp_file=$ARGV[5];
$nx_min=$ARGV[6];
$nx_max=$ARGV[7];
$plot=$ARGV[8];
$npoly2=4;
$ny_cen=0;
$do_med=0;
if ($#ARGV==9) {
    $npoly2=$ARGV[9];
}
if ($#ARGV==10) {
    $npoly2=$ARGV[9];
    $ny_cen=$ARGV[10];
}
if ($#ARGV==11) {
    $npoly2=$ARGV[9];
    $ny_cen=$ARGV[10];
    $do_med=$ARGV[11];
}

print "NY_CEN=$ny_cen $#ARGV\n";
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
	$maximun=0;
	for ($i=$i_min;$i<$i_max;$i++) {
#	    if (($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i-2]>$spec[$i-3])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])&&($spec[$i+2]>$spec[$i+3])) {
#	    if (($spec[$i]-$mean_flux)>($nsigma*$sigma_flux)) {
		if (($spec[$i]>$spec[$i-1])&&($spec[$i-1]>$spec[$i-2])&&($spec[$i]>$spec[$i+1])&&($spec[$i+1]>$spec[$i+2])) {
		    if ($spec[$i]>$maximun) {
			$i_peak=$i;
			$maximun=$spec[$i];
		    }
		}
	    }
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
	$peaks_now[$k]=$x_peak;
#	print "$k $peaks_cen[$k] $peaks_now[$k]\n";

    }


    if ($plot==1) {
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
	    print "Press Enter (A=Automatic, 0=No plot)";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	    if ($command eq "N") {
		$plot=2;
	    }
	}

    }

#
# We determine now the pixels in the new coordinate system:
#
#    for ($i=0;$i<$nx;$i++) {
#	$xpoint[$i]=$i;
#    }
#    my $out_xx_pdl = interpol(pdl(@xpoint), pdl(@peaks_cen), pdl(@peaks_now));
#    for ($i=0;$i<$nx;$i++) {
#	$a_i_new[$i]=$i;
#	$a_i_old[$i]=$out_xx_pdl->at($i);
#	$a_dist[$j][$i]=$a_i_old[$i];
#    }


    $pdl_x_old=pdl(@peaks_now);
    $pdl_x_new=pdl(@peaks_cen);
    ($s_y,$coeff) = fitpoly1d $pdl_x_new,$pdl_x_old,$npoly;
    for ($i=0;$i<$npoly;$i++) {
	$c[$i]=$coeff->at($i);
	$a_c[$j][$i]=$c[$i];	
    }
    $i_old_min=int($c[0]);
    $i_old_max=0;
    for ($kk=0;$kk<$npoly;$kk++) {
	$i_old_max=$i_old_max+$c[$kk]*($nx**$kk);

    }
    for ($i=0;$i<$nx;$i++) {
	$a_i_new[$i]=$i;
	$a_i_old[$i]=0;
	for ($kk=0;$kk<$npoly;$kk++) {
	    $a_i_old[$i]=$a_i_old[$i]+$c[$kk]*($i**$kk);
	}
	$a_dist[$j][$i]=$a_i_old[$i];
#	print "$a_i_new[$i] $a_i_old[$i]\n";
    }



    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$nx,0,$nx,1,0);
	pglabel("X old","X new","Dispersion solution");
	pgsci(1);
	pgpoint($n_peak,\@peaks_now,\@peaks_cen,2);
	pgsci(2);
	pgline($nx,\@a_i_old,\@a_i_new);
	pgclos();
	pgend();
	if (($command ne "A")&&($command ne "0")) {
#	if ($command ne "A") {
	    print "Press Enter (A=Automatic, 0=No plot)";
#	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	}
    }

    
}

#for ($k=0;$k<$npoly;$k++) {
#    my @array;	
#    for ($j=0;$j<$ny;$j++) {
#	$array[$j]=$a_c[$j][$k];
#    }
#    $median=median(@array);
#    $sigma=sigma(@array);
#    for ($j=0;$j<$ny;$j++) {
#	if (abs($array[$j]-$median)>(3*$sigma)) {
#	    $a_c[$j][$k]=$median;
#	}
#    }
#}

$nbox=5;
for ($k=0;$k<$npoly;$k++) {
#    my @array2;
    for ($j=$nbox;$j<($ny-$nbox);$j++) {
	my @array;	
	$nn=0;
	for ($kk=-$nbox;$kk<$nbox;$kk++) {
	    $array[$nn]=$a_c[$j+$kk][$k];
	    $nn++;
	}
	$median=median(@array);
	$sigma=sigma(@array);
#	$array2[$j]=$median;
#    $sigma=sigma(@array);
	if (abs($a_c[$j][$k]-$median)>3*$sigma) {
	    $a_c[$j][$k]=$median;
	}
    }

#    for ($j=$nbox;$j<($ny-$nbox);$j++) {
#    }
}

#$nbox=3;
#for ($k=0;$k<$npoly;$k++) {
#    my @array2;
#   for ($j=$nbox;$j<($ny-$nbox);$j++) {
#	my @array;	
#	$nn=0;
#	for ($kk=-$nbox;$kk<$nbox;$kk++) {
#	    $array[$nn]=$a_c[$j+$kk][$k];
#	    $nn++;
#	}
#	$median=median(@array);
#	$array2[$j]=$median;
#    $sigma=sigma(@array);
#    }
#
#    for ($j=$nbox;$j<($ny-$nbox);$j++) {
#	$a_c[$j][$k]=$array2[$j];
#    }
#}
#
# Smoothing the coefficients;
# 
$command="P";
for ($k=0;$k<$npoly;$k++) {
    my @a_c_now;	
    $min=1e12;
    $max=-1e12;
    for ($j=0;$j<$ny;$j++) {
	$a_y_now[$j]=$j;
#	$a_c_now[$j]=$smoothed->at($k,$j);#$a_c[$j][$k];
	$a_c_now_in[$j]=$a_c[$j][$k];
	$a_c_now[$j]=$a_c[$j][$k];
#	print "$a_y_now[$j] $a_c_now[$j]\n";
	if ($min>$a_c_now[$j]) {
	    $min=$a_c_now[$j];
	}
	if ($max<$a_c_now[$j]) {
	    $max=$a_c_now[$j];
	}
    }
    @a_y_med=median_box($do_med,\@a_y_now);

    $mean=median(@a_c_now);
    $sigma=sigma(@a_c_now);


    if ($sigma>(2*abs($mean))) {
#	$min=$mean-2*abs($mean);
#	$max=$mean+2*abs($mean);
#    } else {
	$min=$mean-1.5*$sigma;
	$max=$mean+1.5*$sigma;
    }

    if ($do_med!=0) {
	print "Doing median filter of coefficient $k\n";
	@a_c_now=median_filter($do_med,\@a_c_now_in);
	print "DONE\n";
    }

    if ($npoly2==0) {
#	my @y2=Derivative2(\@a_y_now,\@a_c_now);
#	my $index=binsearch(\@a_y_now,\@a_y_now);
#	my $index=linsearch(\@a_y_now,\@a_y_now,$index);
#	my $y_interp=spline(\@a_y_now,\@a_c_now,\@y2,$index,\@a_y_now);
	my @a_c_now2=median_box($do_med,\@a_c_now);
	my $spline=new Math::Spline(\@a_y_med,\@a_c_now2);
	for ($j=0;$j<$ny;$j++) {
	    $a_c_new[$j]=$spline->evaluate($a_y_now[$j]);
	}
    } else {
	my $pdl_y_now=pdl(@a_y_now);
	my $pdl_c_now=pdl(@a_c_now);
	my $new_a = fitpoly1d $pdl_y_now,$pdl_c_now,$npoly2;
	for ($j=0;$j<$ny;$j++) {
	    $a_c_new[$j]=$new_a->at($j);
	}
    }

    print "PASO ($npoly2)\n";
    for ($j=0;$j<$ny;$j++) {
	$a_c[$j][$k]=$a_c_new[$j];
    }

    

    if (($plot==1)||($plot==2)) {
	pgbegin(0,"/xs",0,0);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv(0,$ny,$min,$max,0,0);
	pglabel("Y Coord","Pol.Coeff ($k/$npoly)","");
	pgsci(8);
	pgpoint($ny,\@a_y_now,\@a_c_now,4);
	pgsci(1);
	pgpoint($ny,\@a_y_now,\@a_c_now_in,1);
	pgsci(3);
	pgline($ny,\@a_y_now,\@a_c_new);
	pgsci(1);
	pgclos();
	pgend();
	if ($command ne "A") {
	    print "Press Enter";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	}
    }
    print "$k $plot\n";
}

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$a_i_old[$i]=0;
	for ($kk=0;$kk<$npoly;$kk++) {
	    $a_i_old[$i]=$a_i_old[$i]+$a_c[$j][$kk]*($i**$kk);
	}
	$a_dist[$j][$i]=$a_i_old[$i];
    }
}
print "DONE\n";


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
	print "Press Enter";
	$command=<stdin>;
	chop($command);
	if ($command eq "A") {
	    $plot=0;
	}
    }    
}


print "Writting the $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$ny,\@out_array);
#write_rss($out_file,$n_cut,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;


