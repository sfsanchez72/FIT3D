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


if ($#ARGV<0) {
    print "USE: read_Euro3D.pl FILE.e3d \n";
    exit;
}


$infile=$ARGV[0];
my $tmp=$infile."[1]";
my $e3d=rfits($tmp);
$hdr=rfits($tmp,{data=>0});
$crvals=$hdr->{CRVALS};
$cdelts=$hdr->{CDELTS};

@nspax=list($e3d->{NSPAX});
$nb_spec=$#nspax+1;
@spec_len=list($e3d->{SPEC_LEN});
@spec_sta=list($e3d->{SPEC_STA});
print "NB_SPEC=$nb_spec, CRVALS=$crvals, CDELTS=$cdelts\n";
#
# Reading the spectra, and building the RSS
#
($ns_min,$ns_max)=minmax(list($e3d->{SPEC_STA}));
($ne_min,$ne_max)=(0,0);
for ($i=0;$i<$nb_spec;$i++) {
    $n_lim[$i]=$spec_len[$i]+$spec_sta[$i];
    if ($ne_max<$n_lim[$i]) {
	$ne_max=$n_lim[$i];
    }
}
@spax_id=list($e3d->{SPAX_ID});
for ($j=0;$j<$nb_spec;$j++) {
    print "$j/$nb_spec\n";
    my @data_spe=list($e3d->{DATA_SPE});
    for ($i=0;$i<$ne_max;$i++) {
	$rss[$j][$i]=0;
	if (($i>=($spec_sta[$j]-1))&&($i<($n_lim[$i]))) {
	    $rss[$j][$i]=$data_spec[$i-($spec_sta[$j]-1)];
	}
    }
}
$rss->wfits("test.fits");


exit;
