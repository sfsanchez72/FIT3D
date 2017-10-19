#!/usr/bin/perl
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<3) {
    print "USE: get_slice_filter_rss.pl INPUT_RSS.fits POS_TABLE.txt filter.dat output_slice.txt [MASK_LIST] [PLOT]\n";
    exit;
}

$input_rss=$ARGV[0];
$slice=$ARGV[1];
$filter=$ARGV[2];
$out_file=$ARGV[3];
$mask_list=$ARGV[4];
$plot=$ARGV[5];



if ($mask_list eq "") {
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
}

print "$nmask regions\n";




$ns=0;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$ns++;
    } else {
	$header=$line;
    }
}
close(FH);


#
#  Reading filter transmission curve
#

$nf=0;
$sum_f=0;
open(FH,"<$filter");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
        $flux_f[$nf]=$data[2];
        $wave_f[$nf]=$data[1];
    } else {
        $flux_f[$nf]=$data[1];
        $wave_f[$nf]=$data[0];
    }
    $sum_f=$sum_f+$flux_f[$nf];
    $nf++;
}
close(FH);
print "NF=$nf\n";

$b=rfits("$input_rss");
$h=$b->gethdr;
$crval1=$b->hdr->{CRVAL1};
$cdelt1=$b->hdr->{CDELT1};
$crpix1=$b->hdr->{CRPIX1};
($nx,$ny)=$b->dims;


for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval1+$cdelt1*($i-($crpix1-1));    
    $masked[$i]=1;
    for ($j=0;$j<$nmask;$j++) {
	if (($wave[$i]>$start_mask[$j])&&($wave[$i]<$end_mask[$j])) {
	    $masked[$i]=0;	    
	}
    }
}
$pdl_mask=pdl(@masked);
#
# Fraction included
#
$sum_f_in=0;
for ($j=0;$j<$nf;$j++) {
    if (($wave_f[$j]>$wave[0])&&($wave_f[$j]<$wave[$nx-1])) {
	$sum_f_in=$sum_f_in+$flux_f[$j];
    }
}

$fraction_in=$sum_f_in/$sum_f;

my $pdl_flux = interpol(pdl(@wave), pdl(@wave_f), pdl(@flux_f));
$pdl_flux=$pdl_flux*$pdl_mask;
$sum=sum($pdl_flux);
$pdl_flux=($pdl_flux/$sum)*$fraction_in;



print "Fin=$fraction_in\n";
open(OUT,">$out_file");
print OUT "$header\n";
for ($j=0;$j<$ny;$j++) {
#    print "$pdl_flux\n";
    $a=$b->slice(",($j)");
    $c=$a*$pdl_flux;
    $val=sum($c);
    @list_a=list($a);
    @list_flux=list($pdl_flux*$sum);
    @list_c=list($c*$sum);
    ($min,$max)=minmax(@list_a);
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wave[0],$wave[$nx-1],$min,$max,0,0);
	pglabel("Wavelength","Flux","");
	pgsci(1);
	pgline($nx,\@wave,\@list_a);
	pgsci(8);
	pgline($nx,\@wave,\@list_flux);
	pgsci(2);
	pgline($nx,\@wave,\@list_c);
	pgclose;
	pgend;
    }


#    print "$a $pdl_flux\n";
    print OUT "$id[$j] $x[$j] $y[$j] $val\n";
}
close(OUT);


exit;






