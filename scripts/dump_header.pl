#!/usr/bin/perl
#!/usr/bin/perl
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
use Astro::FITS::Header;
#use PDL::Matrix;
$R3DPATH=$ENV{"R3DPATH"};
$mypath=$R3DPATH."/my.pl";
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "$mypath";


@days=(31,28,31,30,31,30,31,31,30,31,30,31);

if ($#ARGV<0) {
    print "USE: read_header.pl fitsfile [add]\n";
    exit;
}

($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);

$row=$ARGV[0];
if ($#ARGV==1) {
    $add=$ARGV[1];
}

$nk=0;
#$pdl->rfits($row);
$h=rfits($row,{data=>0});
open(HDR,">dump_header.tmp");
while( my ($k, $v) = each %$h ) {
    if ($#ARGV==1) {
	$field_now=$add."_".$k;
    } else {
	$field_now=$k;
    }
    $field_now =~ s/\=/\|/g;
    $v =~ s/\=/\|/g;
    $v =~ s/\'//g;
    $field_now =~ s/\=/\|/g;
    print HDR "$field_now | $v \n";
}
close(HDR);
open(HDR,"<dump_header.tmp");
while($line=<HDR>) {
    chop($line);
    if (($line =~ "|")&&($line ne "")) {
	@V=split("/",$line);
	print "$V[0]\n";
    }
}
close(HDR);
exit;
