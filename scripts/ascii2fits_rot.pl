#!/usr/bin/perl
#
#
#

use PGPLOT;
use Statistics::OLS;
use Math::Stat;

use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

$ENV{'PGPLOT_ENVOPT'}="IV";
if ($#ARGV<1) {
    print "USE: ascii2img.pl ASCII FITS\n";
    exit;
}
$ascii=$ARGV[0];
$outfile=$ARGV[1];

system("rm $outfile");

#$pf=$ARGV[2];
#$dpix=$ARGV[3];
#$crval1=$ARGV[4];
#$cd1_1=$ARGV[5];
#$min_y=$ARGV[6];
#$max_y=$ARGV[7];


open(FH,"<$ascii");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($line =~ "#") {
	$npix=$data[1]-1;
	$nb_spec=$data[2];
#	$nb_spec=$data[1];
#	$npix=$data[2];
	$crval1=$data[3];
	$cd1_1=$data[4];
	$cd1_2=$data[5];
	$nx=0;
	$ny=0;
    } else {
#	if ($ny>=$npix) {
#	if ($ny>=$np) {
	if ($ny>=$nb_spec) {
#	if ($ny>=$nb_spec) {
#	    $nx=0;
#	    $ny++;
	    $ny=0;
	    $nx++;
	}	
	$array[$ny][$nx]=$line;
#	$array[$nx][$ny]=$line;
#	$nx++;
	$ny++;
    }
}
print "$npix $nb_spec $crval1 $cd1_1 $nx $ny\n";
#$npix=$nx;
#$n=10;

# We build de image
#
#$naxes={$nb_spec,$npix};
my $status;

fits_clear_errmsg();
$fptr=Astro::FITS::CFITSIO::create_file($outfile,$status);
print "ffinit create new file ($outfile) status = $status\n";
$simple=1;
$bitpix=-32;
$naxis=2;
$naxes2=[$npix,$nb_spec];
$fptr->write_imghdr($bitpix,$naxis,$naxes2,$status);
$fptr->write_key_flt('CRVAL1',$crval1,5,'CRVAL',$status) and print "CRVAL status = $status\n";
$fptr->write_key_flt('CRPIX1',1,5,'CRpix1',$status) and print "CRPIX1 status = $status\n";
$fptr->write_key_flt('CD1_1',$cd1_1,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
$fptr->write_key_flt('CDELT1',$cd1_1,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
$fptr->write_key_flt('CD1_2',$cd1_2,5,'Cd1_2',$status) and print "CD1_2 status = $status\n";


$fptr->write_img(TFLOAT, 1, $nb_spec*$npix+1, \@array, $status);
print "STATUS=$status\n";

$fptr->close_file($status);
print "STATUS=$status\n";


exit;

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub read_naxes() {
    my $file=@_[0];
#    print "FILE=$file";
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
#    print "\nSTATUS=$status\n";
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    ($naxis1,$naxis2) = @$naxes;
    $fptr->close_file($status);
#    print "done\n";
    return $naxes;
}

sub read_img() {
    my $file=@_[0];
#    print "FILE=$file";
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
#    print "\nSTATUS=$status\n";
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    ($naxis1,$naxis2) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_pixnull(Astro::FITS::CFITSIO::TFLOAT(), [1,1], $naxis1*$naxis2, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
#    print "done\n";
#    print "$array[0][0]\n";
    return @array;
}
