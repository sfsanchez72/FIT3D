#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#



use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<4) {
    print "USE: rss2cube_pt_mask.pl input_RSS.fits pos_table.txt mask.txt SMOOTH[0/1] output_cube.fits [NBOX] [SHIFTS_FILE]\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$mask_table=$ARGV[2];
$flag=$ARGV[3];
$output_cube=$ARGV[4];
$nbox=3;
if ($#ARGV==5) {
    $nbox=$ARGV[5];
}

$shift_flag=0;
if ($#ARGV==6) {
    $nbox=$ARGV[5];
    $shift_flag=1;
    $shift_file=$ARGV[6];
}


$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$dpix=$data[1];
$nx_min=1e12;
$ny_min=1e12;
$nx_max=-1e12;
$ny_max=-1e12;

while($line=<PT>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $xx[$n]=$data[1]/$dpix;
    $yy[$n]=$data[2]/$dpix;
    if ($xx[$n]>$nx_max) {
	$nx_max=$xx[$n];
    }
    if ($yy[$n]>$ny_max) {
	$ny_max=$yy[$n];
    }
    if ($xx[$n]<$nx_min) {
	$nx_min=$xx[$n];
    }
    if ($yy[$n]<$ny_min) {
	$ny_min=$yy[$n];
    }
#    print "$xx[$n] $yy[$n] $nx_min,$nx_max $ny_min,$ny_max\n";
    $n++;
}
close(PT);



$nm=0;
open(MASK,"<$mask_table");
while($line=<MASK>) {
    chop($line);
    @data=split(" ",$line);
    $idm[$nm]=$data[0];
    $a_mask[$nm]=$data[1];
    $nm++;
}
close(MASK);

$ns=0;
if ($shift_flag==1) {
    open(SHIFT,"<$shift_file");
    $ns=0;
    while ($line=<SHIFT>) {
	chop($line);
	@data=split(" ",$line);
	$start_shift[$ns]=$data[0];
	$end_shift[$ns]=$data[1];
	$val_shift[$ns]=$data[2];
	$ns++;
    }
    close(SHIFT);
}
print "$ns shifts\n";


for ($i=0;$i<$n;$i++) {
    if ($shift_flag==1) {
	for ($j=0;$j<$ns;$j++) {
	    if (($i>=$start_shift[$j])&&($i<=$end_shift[$j])) {
		$a_shift[$i]=$val_shift[$j];
	    } else {
		$a_shift[$i]=0;
	    }
	}
    }
}



$nx=$nx_max-$nx_min+1;
$ny=$ny_max-$ny_min+1;
for ($j=0;$j<$n;$j++) {
    $xx[$j]=$xx[$j]-$nx_min;
    $yy[$j]=$yy[$j]-$ny_min;
    $ii=intval($xx[$j]);
    $jj=intval($yy[$j]);
#    print "$id[$j] $xx[$j] $yy[$j] $ii,$jj\n";
}

#exit;
#$nx=$ARGV[2];
#$ny=$ARGV[3];


$nax=read_naxes($input_rss);   
@naxis=@$nax;
@a_rss=read_img($input_rss);
$check=$nx*$ny;
($start_w,$delta_w)=read_img_headers($input_rss,["CRVAL1","CDELT1"]);
if ($delta_w==0) {
    $delta_w=1;
}
if ($check!=$naxis[1]) {
    print "Error: The Y-axis of the RSS has $naxis[1] values\n";
    print "It does not correspond with $nx X $ny\n";
    exit;
}


#    for ($j=0;$j<$ny;$j++) {
#	for ($i=0;$i<$nx;$i++) {
#	    $a_cube[$k][$j][$i]=$a_rss[$j+($nx-1-$i)*$ny][$k];	    
for ($nn=0;$nn<$n;$nn++) {
    $i=intval($xx[$nn]);
    $j=intval($yy[$nn]);
    for ($k=0;$k<$naxis[0];$k++) {    
	$new_nn=$nn+$a_shift[$nn];
	if ($new_nn<0) {
	    $new_nn=0;
	}
	if ($new_nn>=$n) {
	    $new_nn=$n-1;
	}
#	$a_cube[$k][$j][$i]=$a_rss[$nn][$k]*$a_mask[$nn];	    	    
	$a_cube[$k][$j][$i]=$a_rss[$new_nn][$k]*$a_mask[$nn];	    	    
	$a_done[$j][$i]=$a_mask[$nn];
    }
}

if ($flag==1) {
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
	    if ($a_done[$j][$i]==0) {
#		print "$i,$j->$a_done[$j][$i]\n";
		for ($k=0;$k<$naxis[0];$k++) {    
		    $j_min=$j-$nbox;
		    $j_max=$j+$nbox;
		    $i_min=$i-$nbox;
		    $i_max=$i+$nbox;
		    if ($j_min<0) { 
			$j_min=0;
		    }
		    if ($j_max>=$ny) { 
			$j_max=$ny;
		    }
		    if ($i_min<0) { 
			$i_min=0;
		    }
		    if ($i_max>=$nx) { 
			$i_max=$nx;
		    }
		    $sum=0;
		    $val=0;
		    for ($jj=$j_min;$jj<$j_max;$jj++) {
			for ($ii=$i_min;$ii<$i_max;$ii++) {
			    $val=$val+$a_cube[$k][$jj][$ii]*$a_done[$jj][$ii];
			    $sum=$sum+$a_done[$jj][$ii];			
			}
		    }
		    if ($sum>0) {
			$a_cube[$k][$j][$i]=$val/$sum;
		    }		    
		}	    
	    }
	}
    }
}



#$nax=[$nx,$ny,$naxis[0]];
$nax=[$nx,$ny,$naxis[0]];
#system("rm $output_cube");
#print "PASO\n";
$pdl=pdl(@a_cube);
$pdl->hdr->{CRPIX1}=1;
$pdl->hdr->{CRVAL1}=$nx_min;
$pdl->hdr->{CDELT1}=$dpix;
$pdl->hdr->{CRPIX2}=1;
$pdl->hdr->{CRVAL2}=$ny_min;
$pdl->hdr->{CDELT2}=$dpix;
$pdl->hdr->{CRPIX3}=1;
$pdl->hdr->{CRVAL3}=$start_w;
$pdl->hdr->{CDELT3}=$delta_w;


$pdl->wfits($output_cube,-32);
#write_fits($output_cube,$nax,3,\@a_cube);
#write_crval_cube($output_cube,[1,0,1,1,0,1,1,$start_w,$delta_w]);
exit;


sub intval {
    my $a=@_[0];
    my $ia=int($a);
    my $d=$a-$ia;
    if ($d>0.5) {
	$ia=$ia+1;
    }
    return $ia;
}
