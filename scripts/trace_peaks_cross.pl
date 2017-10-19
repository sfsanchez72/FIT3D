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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<9) {
    print "USE: trace_peaks_cross.pl RAW.fits Spectral_axis[0/1] PEAKS_FILE x_width y_width plot nplot nsearch TRACE.fits y_shift_limit [shit_ini]\n";
    exit;
}

$y_shift=0;
$infile=$ARGV[0];
$spec_axis=$ARGV[1];
$peaks_file=$ARGV[2];
$width=$ARGV[3];
if ($width==0) {
    $width=1;
}
$y_width=$ARGV[4];
$plot=$ARGV[5];
$nplot=$ARGV[6];
$nsearch=$ARGV[7];
$trace_file=$ARGV[8];
$y_shift_limit=abs($ARGV[9]);
$y_shift=0;
if ($#ARGV==10) {
    $y_shift=$ARGV[10];
}
print "Y_shift=$y_shift_lim $y_shift\n";
#$dmin=$ARGV[7];
#$imin=$ARGV[8];

print "Reading $peaks_file...";
$npeaks=0;
open(FH,"<$peaks_file");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$id[$npeaks]=$data[0];
	$peak_y_pixel[$npeaks]=$data[1]+$y_shift;
	$peak_y_pixel_ini[$npeaks]=$data[1]+$y_shift;
	$peak_y_max[$npeaks]=$data[2]+$y_shift;
	$npeaks++;
    }
}
close(FH);
print "DONE\n";
print "NPEAKS=$npeaks\n";


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

#
# We select the central spectral regions
#
$mean_point=int($nx/2);
#
# We cross-correlate to get the $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #d shift
#
print "Starting Cross-Correlation\n";
$y_shift_tot=$y_shift;
do {
$mm=$mean_point;
my @sec_array;
    $min=1e12;
    $max=-1e12;

for ($j=0;$j<$ny;$j++) {
    $sec_array[$j]=0;
    $peaks_array[$j]=0;
    if ($width>0) {
	for ($i=($mm-$width);$i<($mm+$width);$i++) {
	    $sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	}
    } else {
	for ($i=$mm;$i<($mm+1);$i++) {
	    $sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	}
    }
    if ($min>$sec_array[$j]) {
	$min=$sec_array[$j];
    }
    if ($max<$sec_array[$j]) {
	$max=$sec_array[$j];
    }	        
    $cut[$j]=$j;
}

$sigma=($y_width/(2*2.345));
for ($k=0;$k<$npeaks;$k++) {
    $centroid=$peak_y_max[$k];
    for ($j=0;$j<$ny;$j++) {
	$peaks_array[$j]=$peaks_array[$j]+$sec_array[int($centroid)]*exp(-0.5*(($j-$centroid)/$sigma)**2);	
    }
#    print "$cut[$j] $sec_array[$j] $min $max $peaks_array[$j]\n";
}

if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsubp(1,$nplot);
	for ($i=0;$i<$nplot;$i++) {
	    pgenv(($ny/$nplot)*$i,($ny/$nplot)*($i+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","");
	    pgline($ny,\@cut,\@sec_array);    
	    pgsci(3);
	    pgline($ny,\@cut,\@peaks_array);    
	    pgsch(4.0);           # Set character height
	    for ($k=0;$k<$npeaks;$k++) {
		pgsci(7);
		$x=$peak_y_pixel[$k];
		$y=0.8*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],1);
		pgsci(2);
		$x=$peak_y_pixel_new[$k][$mm];
		$y=0.7*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],5);
		pgsci(4);
		$x=$peak_y_max[$k];
		$y=0.6*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],2);
	    }
	    pgsci(3);
	    pgline(2,[0,$ny],[$imin*$max,$imin*$max]);
	    pgsci(1);
	    pgsch(1.2);           # Set character height
	}
	pgsci(1);
	pgclose;
	pgend;
	if ($command ne "A") {
	    print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	    chop($command);
	}
	if ($command eq "0") {
	    $plot=0;
	}
    }

$k=20;
$n_max=2**$k;
while ($n_max>$ny) {
    $k--;
    $n_max=2**$k;
#    print "$k $n_max $ny\n";
}
#$n_max=512;
for ($j=0;$j<$n_max;$j++) {
    $c_sec_array[$j]=$sec_array[$j];
    $c_peaks_array[$j]=$peaks_array[$j];
}

my $fft1=new Math::FFT(\@c_sec_array);
my $fft2=new Math::FFT(\@c_peaks_array);
my $corr= $fft1->correl($fft2);
@a_tmp=@$corr;
$max=-10000;
$min=1e17;
$sum=0;
for ($j=0;$j<$n_max;$j++) {    
    if ($j<($n_max/2)) {
	$a_corr[$j]=abs($a_tmp[$j+int($n_max/2)]);
    } else {
	$a_corr[$j]=abs($a_tmp[$j-int($n_max/2)]);
    }
#    $a_corr[$j]=$a_tmp[$j];
	$a_xx[$j]=$j;
#	print "$a_corr[$j] $a_xx[$j]\n";
	if ($a_corr[$j]<$min) {
	    $min=$a_corr[$j];
	}
	if ($a_corr[$j]>$max) {
	    $max=$a_corr[$j];
	}
}
$max=-1e12;
for ($j=($n_max/2-$y_width*2);$j<($n_max/2+$y_width*2);$j++) {
    if ($max<$a_corr[$j]) {
	$max=$a_corr[$j];
	$peak_point=$j;
    }   	
    $a=$peak_point-1;
    $b=$peak_point;
    $c=$peak_point+1;
    $fa=-$a_corr[$a];
    $fb=-$a_corr[$b];
    $fc=-$a_corr[$c];
    $den=($fc-2*$fb+$fa);
    $real_peak_point=$c-($b-$a)*(($fc-$fb)/$den+0.5);
}
$y_shift=$real_peak_point-int($n_max/2);
print "Y_shift=$y_shift ($y_shift_limit)\n";

if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($n_max/2-$y_width*2,$n_max/2+$y_width*2,$min,$max,0,0);
    pgline($n_max,\@a_xx,\@a_corr);    
    pgsci(1);
    pgclose;
    pgend;
    if ($command ne "A") {
	print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	chop($command);
    }
    if ($command eq "0") {
	$plot=0;
    }
}

for ($k=0;$k<$npeaks;$k++) {
    $peak_y_pixel_ini[$k]=$peak_y_pixel_ini[$k]+$y_shift;
    $peak_y_max[$k]=$peak_y_max[$k]+$y_shift;
}
$y_shift_tot=$y_shift_tot+$y_shift;
} while (abs($y_shift)>$y_shift_limit);
print "Total Shift=$y_shift_tot\n";
print "DONE\n";

print "We start to trace\n";
#
# 1st Half
#
@peak_y_pixel=@peak_y_pixel_ini;
for ($mm=$mean_point;$mm<($nx-$width);$mm++) {
#    print "COLUMN=$mm\n";
    $min=1e12;
    $max=-1e12;
    my @sec_array;
    for ($j=0;$j<$ny;$j++) {
	$sec_array[$j]=0;
	if ($width>0) {
	    for ($i=($mm-$width);$i<($mm+$width);$i++) {
		$sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	    }
	} else {
	    for ($i=$mm;$i<($mm+1);$i++) {
		$sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	    }
	}
	if ($min>$sec_array[$j]) {
	    $min=$sec_array[$j];
	}
	if ($max<$sec_array[$j]) {
	    $max=$sec_array[$j];
	}	    
	$cut[$j]=$j;
    }
#
# We look for the peaks around the original ones...
#
#print "Looking for peaks with in $nsearch pixels...";

    for ($k=0;$k<$npeaks;$k++) {
	$existe_peak=0;
	for ($j=($peak_y_pixel[$k]-$y_width);$j<($peak_y_pixel[$k]+$y_width);$j++) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($sec_array[$j-$i]<$sec_array[$j-$i-1]) {
		    $peak=0;
		}
		if ($sec_array[$j+$i]<$sec_array[$j+$i+1]) {
		    $peak=0;
		}	       
	    }
	    if ($peak==1) {
		$peak_y_pixel_new[$k][$mm]=$j*1.0;		
		$existe_peak=1;
	    } 

	}
	if ($existe_peak==0) {
	    $peak_y_pixel_new[$k][$mm]=$peak_y_pixel[$k]*1.0;		
	    $peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	} else {
	    $j=$peak_y_pixel_new[$k][$mm];		
	    $a=$j-1;
	    $b=$j;
	    $c=$j+1;
	    $fa=-$sec_array[$a];
	    $fb=-$sec_array[$b];
	    $fc=-$sec_array[$c];
	    $den=($fc-2*$fb+$fa);
	    if ($den!=0) {
#		$peak_y_max_new[$k][$mm]=$a-0.5*((($b-$a)**2)*($fa-$fc)-(($b-$c)**2)*($fb-$fa))/(($b-$a)*($fb-$fc)-($b-$c)*($fb-$fa));	    
		$peak_y_max_new[$k][$mm]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	    } else {
		$peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	    }
	}
	$peak_y_pixel[$k]=$peak_y_pixel_new[$k][$mm];
    }
#
# We plot the section
#
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgsubp(1,$nplot);
	for ($i=0;$i<$nplot;$i++) {
	    pgenv(($ny/$nplot)*$i,($ny/$nplot)*($i+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","");
	    pgline($ny,\@cut,\@sec_array);    
	    pgsch(4.0);           # Set character height
	    for ($k=0;$k<$npeaks;$k++) {
		pgsci(7);
		$x=$peak_y_pixel[$k];
		$y=0.8*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],1);
		pgsci(2);
		$x=$peak_y_pixel_new[$k][$mm];
		$y=0.7*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],5);
		pgsci(4);
		$x=$peak_y_max_new[$k][$mm];
		$y=0.6*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],2);
	    }
	    pgsci(3);
	    pgline(2,[0,$ny],[$imin*$max,$imin*$max]);
	    pgsci(1);
	    pgsch(1.2);           # Set character height
	}

	pgsci(1);
	pgclose;
	pgend;
	if ($command ne "A") {
	    print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	    chop($command);
	}
	if ($command eq "0") {
	    $plot=0;
	}
    }


}

#
# 2nd half
#
@peak_y_pixel=@peak_y_pixel_ini;
for ($mm=$mean_point;$mm>$width;$mm--) {
#    print "COLUMN=$mm\n";
    $min=1e12;
    $max=-1e12;
    my @sec_array;
    for ($j=0;$j<$ny;$j++) {
	$sec_array[$j]=0;
	if ($width>0) {
	    for ($i=($mm-$width);$i<($mm+$width);$i++) {
		$sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	    }
	} else {
	    for ($i=$mm;$i<($mm+1);$i++) {
		$sec_array[$j]=$sec_array[$j]+$in_array[$j][$i];
	    }
	}
	if ($min>$sec_array[$j]) {
	    $min=$sec_array[$j];
	}
	if ($max<$sec_array[$j]) {
	    $max=$sec_array[$j];
	}	    
	$cut[$j]=$j;
    }
#
# We look for the peaks around the original ones...
#
#print "Looking for peaks with in $nsearch pixels...";

    for ($k=0;$k<$npeaks;$k++) {
	$existe_peak=0;
	for ($j=($peak_y_pixel[$k]-$y_width);$j<($peak_y_pixel[$k]+$y_width);$j++) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($sec_array[$j-$i]<$sec_array[$j-$i-1]) {
		    $peak=0;
		}
		if ($sec_array[$j+$i]<$sec_array[$j+$i+1]) {
		    $peak=0;
		}	       
	    }
	    if ($peak==1) {
		$peak_y_pixel_new[$k][$mm]=$j*1.0;		
		$existe_peak=1;
	    } 

	}
	if ($existe_peak==0) {
	    $peak_y_pixel_new[$k][$mm]=$peak_y_pixel[$k]*1.0;		
	    $peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	} else {
	    $j=$peak_y_pixel_new[$k][$mm];		
	    $a=$j-1;
	    $b=$j;
	    $c=$j+1;
	    $fa=-$sec_array[$a];
	    $fb=-$sec_array[$b];
	    $fc=-$sec_array[$c];
	    $den=($fc-2*$fb+$fa);
	    if ($den!=0) {
#		$peak_y_max_new[$k][$mm]=$a-0.5*((($b-$a)**2)*($fa-$fc)-(($b-$c)**2)*($fb-$fa))/(($b-$a)*($fb-$fc)-($b-$c)*($fb-$fa));	    
		$peak_y_max_new[$k][$mm]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
	    } else {
		$peak_y_max_new[$k][$mm]=$peak_y_max[$k];		
	    }
	}
	$peak_y_pixel[$k]=$peak_y_pixel_new[$k][$mm];
    }
#
# We plot the section
#
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgsubp(1,$nplot);
	for ($i=0;$i<$nplot;$i++) {
	    pgenv(($ny/$nplot)*$i,($ny/$nplot)*($i+1),$min,$max,0,0);
	    pglabel("Y-axis","Counts","");
	    pgline($ny,\@cut,\@sec_array);    
	    pgsch(4.0);           # Set character height
	    for ($k=0;$k<$npeaks;$k++) {
		pgsci(7);
		$x=$peak_y_pixel[$k];
		$y=0.8*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],1);
		pgsci(2);
		$x=$peak_y_pixel_new[$k][$mm];
		$y=0.7*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],5);
		pgsci(4);
		$x=$peak_y_max_new[$k][$mm];
		$y=0.6*($max-$min)+$min;		
		pgpoint(1,[$x],[$y],2);
	    }
	    pgsci(3);
	    pgline(2,[0,$ny],[$imin*$max,$imin*$max]);
	    pgsci(1);
	    pgsch(1.2);           # Set character height
	}
	pgsci(1);
	pgclose;
	pgend;
	if ($command ne "A") {
	    print "Press Enter (A for automatic, 0 for no plot):"; 
	    $command=<stdin>;
	    chop($command);
	}
	if ($command eq "0") {
	    $plot=0;
	}
    }
}

for ($mm=0;$mm<$width;$mm++) {
    for ($k=0;$k<$npeaks;$k++) {
	$peak_y_pixel_new[$k][$mm]=$peak_y_pixel_new[$k][$width];
	$peak_y_max_new[$k][$mm]=$peak_y_max_new[$k][$width];
    }
}
for ($mm=($nx-$width);$mm<$nx;$mm++) {
    for ($k=0;$k<$npeaks;$k++) {
	$peak_y_pixel_new[$k][$mm]=$peak_y_pixel_new[$k][$width];
	$peak_y_max_new[$k][$mm]=$peak_y_max_new[$k][$width];
    }
}



#
# We plot the section
#
#$plot=1;
if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    for ($k=0;$k<$npeaks;$k++) {
	pgenv(0,$nx,$peak_y_pixel[$k]-$y_width,$peak_y_pixel[$k]+$y_width,0,0);
	pglabel("X-axis","Y-Axis","$k/$npeaks");
	for ($i=0;$i<$nx;$i++) {
	    pgsci(7);
	    $x=$i;
	    $y=$peak_y_pixel[$k];
	    pgpoint(1,[$x],[$y],1);
	    pgsci(2);
	    $x=$i;
	    $y=$peak_y_pixel_new[$k][$i];
	    pgpoint(1,[$x],[$y],1);
	    pgsci(4);
	    $x=$i;
	    $y=$peak_y_max_new[$k][$i];
	    pgpoint(1,[$x],[$y],1);
	}
	pgsci(1);
	pgclose;
    }
    pgend;
    print "Press Enter (A for automatic, 0 for no plot):"; 
    $command=<stdin>;
}





print "DONE\n";
print "Writting the tracing solution $trace_file\n";
system("rm $trace_file");
#system("cp $input_cube $output_cube");
#system("chmod 755 $output_cube");
#print "Writting the file $output_cube\n";
write_fits($trace_file,[$nx,$npeaks],2,\@peak_y_max_new);
#update_fits($output_cube,$nax,2,\@out_cube);
#write_fits($output_cube,$nax,3,\@out_cube);
print "DONE\n";



exit;



