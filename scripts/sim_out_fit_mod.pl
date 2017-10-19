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


$CRVAL=5061.000000;
$CDELT=2;
$n=45;

open(FH,">spec_base.txt");
for ($i=0;$i<$n;$i++) {
    $wave[$i]=$CRVAL+$CDELT*$i;
    $spec_base[$i]=10;    
    $K=$i+1;
    print FH "$K $wave[$i] $spec_base[$i]\n";
}
close(FH);



$NS=100;
$F=1;
#$CONF="OIII_V500_now.config";
#$call="rm -f sim_out_fit_mod.out";
#system($call);
@FS=(0.1,0.25,0.5,1,1.25,1.5,1.75,2,2.5,3.5,5,7.5,10,15,20);
#@FS=(1,1.25,1.5,1.75,2,2.5,3.5,5,7.5,10,15,20);
open(FOUT,">sim_out.out");
print FOUT "# (1) S/N\n";
print FOUT "# (2) mean\n";
print FOUT "# (3) stddev\n";
print FOUT "# (4) median\n";
print FOUT "# (5) min\n";
print FOUT "# (6) max\n";
foreach $F (@FS) {

    $outfile="sim_out.".$F.".out";
#open(COUT,">sim_out_fit_mod.out");
    open(COUT,">$outfile");
for ($IS=0;$IS<$NS;$IS++) {
    my @flux_now;
    $R=rand(50);
    $S=(-1)**(int(rand(1.5)));
    $vr=5600+$R*$S;
    open(FH,">fix.conf");
    print FH "0 2 0.1 0.1\n";
    print FH "eline\n";
    print FH "5007     0       0       0       -1\n";
    print FH "100      0       -0.1    1e10    -1\n";
    print FH "2.8      0       2.3     3.2     -1\n";
    print FH "$vr      0       4400     7600   -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "poly1d\n";
    print FH "10       0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    close(FH);

    open(FH,">fit.conf");
    print FH "0 2 0.1 0.1\n";
    print FH "eline\n";
    print FH "5007     0       0       0       -1\n";
    print FH "10000    1       -0.1    1e10    -1\n";
    print FH "3        1       2.3     3.2     -1\n";
    print FH "5600     1       4400     7600   -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "poly1d\n";
    print FH "10       1       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       -1e13   1e13    -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    print FH "0        0       0       0       -1\n";
    close(FH);

    $call="fit_spec_back.pl spec_base.txt none 0 10000 none 0 none fix.conf /null > junk.txt";
    system($call);
    open(FH,"<out_mod_res.fit_spectra");    
    open(OUT,">spec_now.txt");
    $id=0;
    while($line=<FH>) {
	chop($line);
        ($w,$f0,$f)=split(" ",$line);
	$id++;
	$R=$F*rand(1);
	$S=(-1)**(int(rand(1.5)));
	$flux_now=$f+$S*$R;
	print OUT "$id $w $flux_now\n";
    }
    close(OUT);
    close(FH);
#    $call="spec_plot.pl spec_now.txt 22/xs 0 100 4500 5300";
#    system($call);
#    <stdin>;
    $call="fit_spec_back.pl spec_now.txt none 0 10000 none 0 none fit.conf /xs > junk.txt";
#    $call="fit_spec_back.pl spec_now.txt none 0 10000 none 0 none fit.conf /null > junk.txt";
    system($call);
#    <stdin>;
    $call="cat out.fit_spectra | grep 5007.0000 > sim_out.now";
    system($call);
    open(OUT,"<sim_out.now");
    $line=<OUT>; chop($line);
    @data=split(" ",$line);
    $delta=$data[7]-$vr;
    close(OUT);
    print COUT "$IS $delta $data[7] $vr\n";
#    print  "$IS $delta $data[7] $vr\n";

}
    
close(COUT);
    $call="table_hist.pl ".$outfile." 1 ".$outfile.".ps/CPS 10 > out_now.txt";
    system($call);
    $SN=100/(2.8*2.345*($F+0.1));
    open(FH,"<out_now.txt");
    $line=<FH>;
    close(FH);
    print FOUT "$SN $line";

#    $call="table_hist.pl ".$outfile." 1 33/xs 10";
#    system($call);
}
    close(FOUT);
exit;
