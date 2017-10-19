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
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;
use PDL::Image2D;
use PDL::ImageRGB;
use PDL::Graphics::PGPLOT;
use PDL::Ufunc;

#$ENV{PGPLOT_FOREGROUND} = "black";
#$ENV{PGPLOT_BACKGROUND} = "white";

#$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<13) {
    print "USE: maps_orion.pl map.fits min max bright contrast Label factor dev VAL_TO_WHITE MAP_CONT NLEVELS MIN STEP NWHITE[49] [TABLE REVERSE] [ROT]\n";
    exit;
}
$mapfiles=$ARGV[0];
($Rfile,$Gfile,$Bfile)=split(/\,/,$mapfiles);
#print "$Rfile | $Gfile | $Bfile \n"; exit;
$min=$ARGV[1];
$max=$ARGV[2];
$bright=$ARGV[3];
$contrast=$ARGV[4];
$label=$ARGV[5];
$factor=$ARGV[6];
$dev=$ARGV[7];
$vw=$ARGV[8];
$map_cont=$ARGV[9];
$nlevels=$ARGV[10];
$min_cont=$ARGV[11];
$steps=$ARGV[12];
$NW=$ARGV[13];
if ($dev !~ "TPNG") { 
    $ENV{PGPLOT_FOREGROUND} = "black";
    $ENV{PGPLOT_BACKGROUND} = "white";
} else {
    $ENV{PGPLOT_FOREGROUND} = "white";
    $ENV{PGPLOT_BACKGROUND} = "white";
}

$table="idl5"; $reverse=1; 
if ($#ARGV==15){
    $table=$ARGV[14];
    $reverse=$ARGV[15];
}
$rot=0;
if ($#ARGV==16){
    $table=$ARGV[14];
    $reverse=$ARGV[15];
    $rot=$ARGV[16];
}




#$table="smooth2"; $reverse=0; 
#$table="real"; $reverse=0; 
#$bright=0.7; 
#$contra=0.5; 
$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");
#$pdl_mask=rfits("mask_sc.fits");


$pdl_R_tmp=rfits($Rfile);
$pdl_G_tmp=rfits($Gfile);
$pdl_B_tmp=rfits($Bfile);
$pdl_cont_tmp=rfits($map_cont);

$is_zero=sum(abs($pdl_R_tmp));


if ($rot==0) {
    $pdl_R=$pdl_R_tmp;
    $pdl_G=$pdl_G_tmp;
    $pdl_B=$pdl_B_tmp;
    $pdl_cont=$pdl_cont_tmp;
} else {
#    $pdl_R=$pdl_R_tmp->rot2d($rot,0,1);
#    $pdl_cont=$pdl_cont_tmp->rot2d($rot,0,1);
    $pdl_R=transpose($pdl_R_tmp);
    $pdl_G=transpose($pdl_G_tmp);
    $pdl_B=transpose($pdl_B_tmp);
    $pdl_cont=transpose($pdl_cont_tmp);
}


#$crval1=$pdl_R->hdr->{CRVAL1};
#$cdelt1=$pdl_R->hdr->{CDELT1};
#$crval2=$pdl_R->hdr->{CRVAL2};
#$cdelt2=$pdl_R->hdr->{CDELT2};
#print "$crval1 $cdelt1 $crval2 $cdelt2\n";
($nx,$ny)=$pdl_R->dims;
$h=$pdl_R->gethdr;
$crval1=-138.839996337891;
$cdelt1=1;
$crval2=-210.726806640625;
$cdelt2=1;


$pdl_R=$pdl_R*$factor;


$NNx=$nx-1;
$NNy=$ny-1;
$rpdl_R=$pdl_R->slice("$NNx:0,$NNy:0");
#@mapR=list($rpdl_R);
$rpdl_G=$pdl_G->slice("$NNx:0,$NNy:0");
#@mapG=list($rpdl_G);
$rpdl_B=$pdl_B->slice("$NNx:0,$NNy:0");
#@mapB=list($rpdl_B);
$rpdl_cont=$pdl_cont->slice("$NNx:0,$NNy:0");
@cont=list($rpdl_cont);
$crval1=-30;
$crval2=-40;
$x_max=$crval1;
$x_min=$crval1+$cdelt1*$nx;
$y_min=$crval2;
$y_max=$crval2+$cdelt2*$ny;

@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);


#dev($dev);

dev($dev);
#pgbeg(0,$dev,1,1);
pgscf(2.0);
pgsch(1.2);
pgenv($x_min,$x_max,$y_min,$y_max,1,0);
#env($x_min,$x_max,$y_min,$y_max);
#env(0,$nx,0,$ny);
if ($dev !~ "TPNG") { 
    pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
}


my ($lo,$hi);
my ($oldlo,$oldhi);
my ($nlevelm1,$imag,$nimag);
my ($out,$lut);
my ($r,$g,$b);
my ($levels);
my ($i,@cr,@cg,@cb);
pgqcir($lo,$hi);

$nlevelm1=($hi-$lo)*1.;

$imag=cat($rpdl_R,$rpdl_G,$rpdl_B)->xchg(0,2)->xchg(1,2);
#my $test=$imag->slice("(0),:,(0)");
#print "$test\n";
$min_imag=$min;
$max_imag=$max;
#($min_imag,$max_imag)=minmax($imag);
#print "$min_imag,$max_imag $nlevelm1\n";
#$nimag=byte(PDL::copy((((float($imag)-$min_imag)/($max_imag-$min_imag))*($nlevelm1))));
#$nimag=byte((((float($imag)-$min_imag)/($max_imag-$min_imag))*($nlevelm1)));

$fimag=(((float($imag)-$min_imag)/($max_imag-$min_imag))*(255*1000));
$nimag=byte($fimag);

#my $test=$nimag->slice("(0),:,:");
#$test->wfits("test.fits");

#print "$test\n $nlevelm1\n";
#$nimag=byte(PDL::copy((($imag-min($imag))/(max($imag)-min($imag))*($nlevelm1))));

#$nimag_float=(PDL::copy((($imag-min($imag))/(max($imag)-min($imag))*(255))));
#$nimag_float->wfits("test.fits");
#$nimag=byte($nimag_float);


($out,$lut)=cquant($nimag,$nlevelm1+1);
#print "nlevelsm1=$nlevelm1 $hi $lo\n";
#exit;

#$out->wfits("test.fits");

$out=float($out)/$nlevelm1;
$lut=float($lut)/$nlevelm1;
$r=$lut->slice('(0),:');
$g=$lut->slice('(1),:');
$b=$lut->slice('(2),:');
$levels=sequence(float,($nlevelm1+1))/$nlevelm1;

ctab($levels,$r,$g,$b);
@al=list($levels);
@ar=list($r);
@ag=list($g);
@ab=list($b);

for ($i=0;$i<$#ar+1;$i++) {
    print "$i $al[$i] $ar[$i] $ag[$i] $ab[$i]\n";
}

#print "@al\n";
#pgctab(\@al,\@ar,\@ag,\@ab,$nlevelm1,$bright,$contrast);
#imag($out,0,1,$tr);


#@a_out=list($out);
#pgimag(\@a_out,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);
#pgimag(\@a_out,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);
#$out->wfits("test.fits");
imag($out,0,1,pdl(@tr));

pgsci(0);
pgsch(1);
#pgsch(1.1);
@tmp;
$nt=0;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$mapR[$i+$j*$nx];
#	if ($mapR[$i+$j*$nx]<=0) {
	if ($val<=0) {
	    $X=$crval1+$cdelt1*$i+1;
	    $Y=$crval2+$cdelt2*$j+1;
	#    print "$i,$j $X,$Y MAPR=$mapR[$i+$j*$nx]\n";
#	    pgpoint(1,[$X],[$Y],16);
	} else {
	    if ($val>6) {
		$tmp[$nt]=$val;
#		print "$tmp[$nt] $nt\n";
		$nt++;
	    }
	}
    }
}
$mean=mean(@tmp);
$sigma=sigma(@tmp);
$sum=$nt*$mean;
print "FLUX=$sum\n";
print "MEAN=$mean+-$sigma\n";


pgsch(3);
pgsci(3);
pgslw(3);
#pgpoint1,[3.65],[-3.13],3);
#pgpoint(1,0,0,3);
pgsch(1.2);
pgslw(1);

pgsci(1);
pgsch(1.3);
if ($dev !~ "TPNG") { 
    pgsci(1);
    pgsch(1.3);
    if ($rot==0) {
	pgptxt(-15,28,0,0.5,$label);
    } else {
	pgptxt(-15,24,0,0.5,$label);
    }
#pgptxt(-13,28,0,0.5,$label);
}

#pgptxt(-13,28,0,0.5,$label);
pgsci(1);
pgsch(1);

pgsci(1);
pgslw(3);
#$nlevels=10;
#$min_cont=$ARGV[11];
#$steps=$ARGV[12];
for ($i=0;$i<$nlevels;$i++) {    
    $levels[$i]=$min_cont+$steps*($i**1.5);
#    print "$i $levels[$i]\n";
}

pgsci(1);
pgslw(3);
if ($dev !~ "TPNG") { 
    pgsci(1);
    pgslw(3);
} else { 
    if ($reverse==0) {
	pgsci(4);
    } else {
	pgsci(4);
    }
    pgslw(8);
}



if ($is_zero!=0) {
    pgcont(\@cont,$nx,$ny,1,$nx,1,$ny,\@levels,$nlevels,\@tr);
}
pgsci(1);
pgslw(1);

if ($dev !~ "TPNG") { 
#pgptxt(-85,160,0,0,$label);
#pgptxt(-110,160,0,0,$label);
pgwedg("RI",0.0,4.5,$min,$max,"");
#pgwedg("RI",0,1.0,1,500,"");
pgbox("SBC",0,0,"SBC",0,0);
}




pgclos();
pgend();



exit;


sub mean {
    local(@data)=@_;
    my $sum=0;
    my $j;
    for ($j=0;$j<$#data;$j++) {
	$sum=$sum+$data[$j];
    }
    my $mean = $sum/$#data;
    return $mean;
}

sub sigma {
    local(@data)=@_;
    my $mean = mean(@data);
    my $stddev = 0;
    my $j;
    my $sum=0;
    for ($j=0;$j<$#data;$j++) {
	$sum=$sum+($data[$j]-$mean)**2;
    }
    if ($#data>1) {
	$sum=$sum/($#data-1);
    }
    $stddev=sqrt($sum);
    return $stddev;
}

