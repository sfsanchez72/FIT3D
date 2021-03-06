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
if ($#ARGV<4) {
    print "USE: maps2fits.pl ASCII FITS NMAP NX NY\n";
    exit;
}
$ascii=$ARGV[0];
$outfile=$ARGV[1];
$n_map=$ARGV[2];
$nx_tot=$ARGV[3];
$ny_tot=$ARGV[4];

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
    if ($line !~ "#") {	
	@data=split(" ",$line);
	if ($nx>=$nx_tot) {
	    $nx=0;
	    $ny++;
	}	
	$array[$ny][$nx]=$line[$n_map];
	$nx++;
    }
}
#print "$npix $nb_spec $crval1 $cd1_1 $nx $ny\n";
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
$naxes2=[$nx_tot,$ny_tot];
$fptr->write_imghdr($bitpix,$naxis,$naxes2,$status);
#$fptr->write_key_flt('CRVAL1',$crval1,5,'CRVAL',$status) and print "CRVAL status = $status\n";
#$fptr->write_key_flt('CRPIX1',1,5,'CRpix1',$status) and print "CRPIX1 status = $status\n";
#$fptr->write_key_flt('CD1_1',$cd1_1,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
#$fptr->write_key_flt('CD1_2',1,5,'Cd1_2',$status) and print "CD1_2 status = $status\n";


$fptr->write_img(TFLOAT, 1, $nx_tot*$ny_tot+1, \@array, $status);
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
