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
#use PDL::Graphics::TriD;
#use PDL::Graphics::TriD::Image;
use PDL::Fit::Gaussian;
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";

$R3DPATH=$ENV{"R3DPATH"};
$mypath=$R3DPATH."/my.pl"; 
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "$mypath";


if ($#ARGV<4) {
    print "USE: radial_sum_cube.pl CUBE.fits Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]\n";
    print "It will use images coordinates: 0,nx 0,ny\n";
    exit;
}

$input=$ARGV[0];
$Dr=$ARGV[1];
$x_c=$ARGV[2];
$y_c=$ARGV[3];
$output=$ARGV[4];
$e_output="e_".$output;
$weight_output="weight.".$output;
$plot=0;
if ($#ARGV==5) {
    $plot=$ARGV[5];
}

$r_max=0;
$n=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;

$pdl_cube=rfits("$input");
($NX,$NY,$NZ)=$pdl_cube->dims();
$y_c=$NY-$y_c;


$h=$pdl_cube->gethdr;
$crval=$pdl_cube->hdr->{"CRVAL3"};
$cdelt=$pdl_cube->hdr->{"CDELT3"};
$crpix=$pdl_cube->hdr->{"CRPIX3"};

$pdl_mask=ones($NX,$NY,$NZ);                                                                                                                             
$extend=$h->{"EXTEND"};
if ($extend==1) {
    $mask_ext=$input."[3]"; 
    my $b_mask=rfits($mask_ext);
    $pdl_mask=$pdl_mask-$b_mask;
    $e_ext=$input."[1]"; 
    $pdl_e=rfits($e_ext);
    $w_ext=$input."[2]"; 
    $pdl_w=rfits($w_ext);
}
    
$pdl_cube=$pdl_cube*$pdl_mask;

$nz_med=int($NZ/2);
$start_i=$nz_med-50;
$end_i=$nz_med+50;
#$npix=$end_i-$start_i+1;
my $a=$pdl_cube->slice(",,$start_i:$end_i");
my $c=average($a->xchg(0,2));
my $d=$c->xchg(0,1);
my $a_mask=$pdl_mask->slice(",,$start_i:$end_i");
my $c_mask=average($a_mask->xchg(0,2));
my $d_mask=$c_mask->xchg(0,1);
$pdl_weight=$d/$d_mask;



for ($i=0;$i<$NZ;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
}
$wmin=$wave[0];
$wmax=$wave[$NZ-1];

$shape="R";
$size=1;

$nx=$NZ;
$ny=$NX*$NY;

$pdl_in=zeroes($nx,$ny);
$pdl_in_e=zeroes($nx,$ny);
$pdl_in_mask=zeroes($nx,$ny);
$pdl_in_weight=zeroes($ny);
$r_max=0;
$n=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;
for ($j=0;$j<$NY;$j++) {
    for ($i=0;$i<$NX;$i++) {
	$px=$i+($NY-1-$j)*$NX;
#	print "($NX,$NY,$NZ) ($nx,$ny) $i,$j $px $n\n";
        $t = $pdl_in->slice(":,($px)");
	$t .= $pdl_cube->slice("($i),($j),:");
        $t_e = $pdl_in_e->slice(":,($px)");
	$t_e .= $pdl_e->slice("($i),($j),:");
        $t_m = $pdl_in_mask->slice(":,($px)");
	$t_m .= $pdl_mask->slice("($i),($j),:");
	$t_w = $pdl_weight->at($i,$j);
	set($pdl_in_weight,$px,$t_w);
	$id[$n]=$i+$j*$nx+1;
	$x[$n]=$i;
	$y[$n]=$j;
	$r[$n]=sqrt(($x[$n]-$x_c)**2+($y[$n]-$y_c)**2);    
	$nk=$r[$n]/$Dr;
	$c[$n]=2+int($nk/2);
	$s[$n]=3+int($nk/2);
	$S[$n]=2.5-$nk/15;

	if ($r_max<$r[$n]) {
	    $r_max=$r[$n];	
	}
	if ($x_min>$x[$n]) {
	    $x_min=$x[$n];
	}
	if ($y_min>$y[$n]) {
	    $y_min=$y[$n];
	}
	if ($x_max<$x[$n]) {
	    $x_max=$x[$n];
	}
	if ($y_max<$y[$n]) {
	    $y_max=$y[$n];
	}
	

	$n++;

	
    }
}
$nr=int($r_max/$Dr)+1;
$pdl_out=zeroes($nx,$nr);
$pdl_out_e=zeroes($nx,$nr);
$pdl_out_mask=zeroes($nx,$nr);




for ($i=0;$i<$nr;$i++) {
    $r_min=$Dr*$i;
    $r_max=$Dr*($i+1);
    $nsum=0;
    $sum_all=0;
    $t=$pdl_out->slice(":,($i)");
    $t_e=$pdl_out_e->slice(":,($i)");
    $t_m=$pdl_out_mask->slice(":,($i)");
    my @spaxels;
    $nspaxels=0;
#    print "$r_min $r_max\n";


    $flux_max=-1e20;
    $flux_min=1e20;
    $nx1=int(0.25*$nx);
    $nx2=int(0.6*$nx);
    $slice=$pdl_in->slice("$nx1:$nx2,");
    $image = sumover $slice;#->xchg(0,1);
    for ($ns=0;$ns<$ny;$ns++) {
	$flux[$ns]=$image->at($ns);
	if ($flux[$ns] eq "nan") {
	    $flux[$ns]=0;
	}
	if ($flux_min>$flux[$ns]) {
	    $flux_min=$flux[$ns];
	}
	if ($flux_max<$flux[$ns]) {
	    $flux_max=$flux[$ns];
	}
    }
    if ($flux_min eq "-inf") {
	$flux_min=0;
    }

    for ($j=0;$j<$ny;$j++) {
	if ($r[$j]<$r_max) {
#	if (($r[$j]>=$r_min)&&($r[$j]<$r_max)) {
	    $spaxels[$nspaxels]=$j;
	    $nspaxels++;
	    $slice=$pdl_in->slice(":,($j)");
	    $slice_e=$pdl_in_e->slice(":,($j)");
	    $slice_mask=$pdl_in_mask->slice(":,($j)");
	    $weight=$pdl_in_weight->at($j);
	    if ($flux[$j]>0) {
		$t.=$t+$slice;
		$t_e.=$t_e+$slice_e;
		$t_m.=$t_m+$slice_mask*$weight;
		$nsum=$nsum+$weight;
	    }
	}
    }
    $NN[$i]=$nsum;

#    $t=$pdl_out->slice(":,($i)");
#    $t_m=$pdl_out_mask->slice(":,($i)");
    if ($NN[$i]>0) {
	$t .= ($t/$t_m)*$NN[$i];
	$t_e .= ($t_e/$t_m)*sqrt(abs($NN[$i]));
	$t_m .= $t_m/$NN[$i];
    }

#    $t.=$t/$nspaxels;
#    print "$nspaxels\n";
    if ($plot>0) {
	#$size=1;
	if ($plot==1) {
	    $dev="/xs";
	} else {
	    $dev="radial_sum_cube_".$i.".ps/CPS";
	}

#	$flux_min=-1;
#	$flux_max=3*$flux_max;
#	print "$flux_min,$flux_max\n";
	$bright=0.55; 
	$contrast=0.55;
	$table="idl5"; 
	$reverse=1; 
	pgbeg(0,$dev,1,1);
	pgscf(2.0);
	pgsch(1.2);
	pgpap(14.0,0.4);
	pgsubp(2,1);

	($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
	$nc=$pl->getdim(0);
	for ($jj=0;$jj<$nc;$jj++) {
	    $l[$jj] = $pl->slice($jj)->sclr;
	    $rr[$jj] = $pr->slice($jj)->sclr;
	    $g[$jj] = $pg->slice($jj)->sclr;
	    $b[$jj] = $pb->slice($jj)->sclr;
	}
    
    
	#$rr[0]=1;  
	#$g[0]=1;  
	#$b[0]=1;  
	
	pgscir(50,100);
	pgctab(\@l,\@rr,\@g,\@b,$nc,$bright,$contrast);
	
	pgenv($x_max,$x_min,$y_min,$y_max,1,0);
	pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
	

	for ($ii=0;$ii<$ns;$ii++) {
	$color[$ii]=50+ceil((($flux[$ii]-$flux_min)/($flux_max-$flux_min))*50);

	if ($color[$ii] eq "nan") {
	    $color[$ii]=0;
	}

	pgsci($color[$ii]);

	if ($shape eq "C") {
	    pgsfs(1);
	    pgcirc($x[$ii],$y[$ii],$size);
	    pgsfs(2);
	    pgsci(1);
	    
	    pgcirc($x[$ii],$y[$ii],$size);

	}
	if ($shape eq "R") {

	    pgsfs(1);
	    pgrect($x[$ii]-$size/2,$x[$ii]+$size/2,$y[$ii]-$size/2,$y[$ii]+$size/2);
	    pgsfs(2);
	    pgsci(1);
	    
	    pgrect($x[$ii]-$size/2,$x[$ii]+$size/2,$y[$ii]-$size/2,$y[$ii]+$size/2);
	}

	}

    pgsci(3);
	pgsch(1.2);
	$sys=3;
    for ($ii=0;$ii<$nspaxels;$ii++) {
	$jj=$spaxels[$ii];
	pgpoint(1,[$x[$jj]],[$y[$jj]],$sys);
    }

	pgsci(1);
	pgwedg("RI", 0.1, 5, $flux_min, $flux_max, "");
	
	@spec=list($t);
	if ($i==0) {
	    $med=median(@spec);
	    $y_max_s=3*$med;
	    $y_min_s=-0.5*$med;
	} else {
	    $med=median(@spec);
	    $y_max_s_now=3*$med;
	    $y_min_s_now=-0.5*$med;
	    if ($y_max_s_now>$y_max_s) {
		$y_max_s=$y_max_s_now;
	    }
	    if ($y_min_s_now<$y_min_s) {
		$y_min_s=$y_min_s_now;
	    }
	}
#	print "$y_min_s,$y_max_s $y_min_s_now,$y_max_s_now\n";
	pgenv($wmin,$wmax,$y_min_s,$y_max_s,0,0);
	pglabel("Wavelength","Flux","");		   
	for ($ii=0;$ii<$i+1;$ii++) {
	    pgsci($ii+1);
	    $t_now=$pdl_out->slice(":,($ii)");
	    @spec=list($t_now);
	    pgline($nx,\@wave,\@spec);	
	}
	








    pgclos();
    pgend();

	if ($plot==1) {
	    print "Press Enter\n"; <stdin>;
	}
}





}

#$$h{"NAXIS2"}=$nr;
#$pdl_out->sethdr($h);
#$h = {NAXIS=>2, NAXIS1=>$nx, NAXIS=>$nr, COMMENT=>"Sample FITS-style header"};

#for ($i=0;$i<$nr;$i++) {
#    $t=$pdl_out->slice(":,($i)");
#    $t_m=$pdl_out_mask->slice(":,($i)");
#    if ($NN[$i]>0) {
#	$t .= ($t/$t_m)*$NN[$i];
#    }
#    print "$t \n $t_m\n $NN[$i]\n";
#}
#$pdl_out=$pdl_out/$pdl_out_mask;

$pdl_out->hdr->{CRPIX1}=$crpix; 
$pdl_out->hdr->{CRVAL1}=$crval; 
$pdl_out->hdr->{CDELT1}=$cdelt; 
$pdl_out->wfits($output);

$pdl_out_e->hdr->{CRPIX1}=$crpix; 
$pdl_out_e->hdr->{CRVAL1}=$crval; 
$pdl_out_e->hdr->{CDELT1}=$cdelt; 
$pdl_out_e->wfits($e_output);

#$e_output="e_".$output;

$pdl_out_mask->hdr->{CRPIX1}=$crpix; 
$pdl_out_mask->hdr->{CRVAL1}=$crval; 
$pdl_out_mask->hdr->{CDELT1}=$cdelt; 
$pdl_out_mask->wfits($weight_output);


#$weight_output="weight.".$output;

#
# We plot
#


exit;
