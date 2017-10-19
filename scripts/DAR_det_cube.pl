#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
use PDL::Transform;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<5) {
    print "USE: DAR_det_cube.pl input_cube.fits X Y NP_X NP_Y output_cube [N_START] [N_END] [Aperture] [SCALE] [DELTA_X DELTA_Y] [RECENTER=0/1] [Y_WIDTH]\n";
    exit;
}


$input_cube=$ARGV[0];
$x0=$ARGV[1];
$y0=$ARGV[2];
$npoly_x=$ARGV[3];
$npoly_y=$ARGV[4];
$output_cube=$ARGV[5];
$ps="DAR_".$input_cube;
$ps =~ s/.fits/.ps/;
$dev=$ps."/CPS";
$n_start=0;
$n_end=0;
if ($#ARGV==7) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
}
$aperture=5;
if ($#ARGV==8) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
    $aperture=$ARGV[8];
}

$scale=1;
if ($#ARGV==9) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
    $aperture=$ARGV[8];
    $scale=$ARGV[9];
}
$delta_x=0;
$delta_y=0;
if ($#ARGV==11) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
    $aperture=$ARGV[8];
    $scale=$ARGV[9];
    $delta_x=$ARGV[10];
    $delta_y=$ARGV[11];
}


$recenter=1;
if ($#ARGV==12) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
    $aperture=$ARGV[8];
    $scale=$ARGV[9];
    $delta_x=$ARGV[10];
    $delta_y=$ARGV[11];
    $recenter=$ARGV[12];
}

$yw=1;
if ($#ARGV==13) {
    $n_start=$ARGV[6];
    $n_end=$ARGV[7];
    $aperture=$ARGV[8];
    $scale=$ARGV[9];
    $delta_x=$ARGV[10];
    $delta_y=$ARGV[11];
    $recenter=$ARGV[12];
    $yw=$ARGV[13];
    if ($yw<1) {
	$yw=1;
    }
}

#print "$n_start,$n_end $#ARGV\n";

$pdl_in=rfits($input_cube);
($nx,$ny,$nz)=$pdl_in->dims;
$H=$pdl_in->gethdr;

$nx_0=int($nx*$scale);
$ny_0=int($ny*$scale);

print "$scale $nx_0 $ny_0 $nz\n";
$pdl_out=zeroes($nx_0,$ny_0,$nz);

$x_c=zeroes($nz);
$y_c=zeroes($nz);
$flux=zeroes($nz);

if ($n_start<0) {
    $n_start=0;
}
if ($n_end>$nz) {
    $n_end=$nz;
}
#print "$n_start,$n_end $nz\n";
$kmed=int($nz/2);
for ($k=$kmed;$k<$nz;$k++) {

#    if ($k==$kmed) {
	$x1=$x0;
	$y1=$y0;
#    }     else {
#	if (($x_val>0)&&($y_val>0)) {
#	    $x1=$x_val;
#	    $y1=$y_val;
#	}
#    }
    #
    $i_min=int($x1-$aperture);
    $i_max=int($x1+$aperture);
    $j_min=int($y1-$aperture);
    $j_max=int($y1+$aperture);
    


    if ($j_min<0) {
	$j_min=0;
    }
    if ($j_max<0) {
	$j_max=0;
    }
    if ($j_max>=$ny) {
	$j_max=$ny-1;
    }
    if ($j_min>=$ny) {
	$j_min=$ny-1;
    }
    if ($i_min<0) {
	$i_min=0;
    }
    if ($i_max<0) {
	$i_max=0;
    }
    if ($i_max>=$nx) {
	$i_max=$nx-1;
    }
    if ($i_min>=$nx) {
	$i_min=$nx-1;
    }
    $x_val=0;
    $y_val=0;
    $sum=0;
    for ($i=$i_min;$i<$i_max+1;$i++) {
	for ($j=$j_min;$j<$j_max+1;$j++) {
	    $k_min=$k-$yw;
	    $k_max=$k+$yw;
	    if ($k_min<0) {
		$k_min=0;
		$k_max=2*$yw;
	    }
	    if ($k_max>$nz-1) {
		$k_min=$nz-1-2*$yw;
		$k_max=$nz-1;
	    }

	    for ($kk=$k_min;$kk<$k_max;$kk++) {
		$val=$pdl_in->at($i,$j,$kk);
		$val=($val**4);
		if ($val>0) {
		    $x_val=$x_val+$i*$val;
		    $y_val=$y_val+$j*$val;
		    $sum=$sum+$val;	 
		}
	    }
   
	}


    }
    if ($sum>0) {
	$x_val=$x_val/$sum;
	$y_val=$y_val/$sum;
    }
    set($x_c,$k,$x_val);
    set($y_c,$k,$y_val);
    set($flux,$k,$sum);
    $id[$k]=$k;
    
#print "$kmed $k $x0,$y0 $x_val,$y_val\n";
    $x_old=$x_val;
    $y_old=$y_val;
}

for ($k=$kmed;$k>-1;$k--) {
#    if ($k==$kmed) {
	$x1=$x0;
	$y1=$y0;
#    } 
#    else {
#	if (($x_val>0)&&($y_val>0)) {
#	    $x1=$x_val;
#	    $y1=$y_val;
#	}
#    }
#    if ($k==$kmed) {
#	$x1=$x0;
#	$y1=$y0;
#    } else {
#	$x1=$x_val;
#	$y1=$y_val;
#    }

    #  my $img_in=$pdl_in->slice(":,:,$k");
    $i_min=int($x1-$aperture);
    $i_max=int($x1+$aperture);
    $j_min=int($y1-$aperture);
    $j_max=int($y1+$aperture);

    if ($j_min<0) {
	$j_min=0;
    }
    if ($j_max<0) {
	$j_max=0;
    }
    if ($j_max>=$ny) {
	$j_max=$ny-1;
    }
    if ($j_min>=$ny) {
	$j_min=$ny-1;
    }
    if ($i_min<0) {
	$i_min=0;
    }
    if ($i_max<0) {
	$i_max=0;
    }
    if ($i_max>=$nx) {
	$i_max=$nx-1;
    }
    if ($i_min>=$nx) {
	$i_min=$nx-1;
    } 
    $x_val=0;
    $y_val=0;
    $sum=0;
    for ($i=$i_min;$i<$i_max+1;$i++) {
	for ($j=$j_min;$j<$j_max+1;$j++) {
#	    print "$k $i,$j,$k $nx,$ny,$nz $x0 $y0 [$i_min:$i_max,$j_min:$j_max] $x_val $y_val\n";

	    $k_min=$k-$yw;
	    $k_max=$k+$yw;
	    if ($k_min<0) {
		$k_min=0;
		$k_max=2*$yw;
	    }
	    if ($k_max>$nz-1) {
		$k_min=$nz-1-2*$yw;
		$k_max=$nz-1;
	    }

	    for ($kk=$k_min;$kk<$k_max;$kk++) {
		$val=$pdl_in->at($i,$j,$kk);
		$val=($val**4);
		if ($val>0) {
		    $x_val=$x_val+$i*$val;
		    $y_val=$y_val+$j*$val;
		    $sum=$sum+$val;	 
		}
	    }

	}


    }
    if ($sum>0) {
	$x_val=$x_val/$sum;
	$y_val=$y_val/$sum;
    }
    set($x_c,$k,$x_val);
    set($y_c,$k,$y_val);
    set($flux,$k,$sum);
    $id[$k]=$k;
#    print "$kmed $k $x0,$y0 $x_val,$y_val\n";
}

@X=list($x_c);
@Y=list($y_c);


#
# We mask the wrong values 
#
#print "K=$k $n_start $n_end\n";
$NZ=0;
for ($k=$n_start;$k<$n_end;$k++) {
    $point[$NZ]=$k;
    $y[$NZ]=$Y[$k];
    $x[$NZ]=$X[$k];
    $flux_now[$NZ]=$flux->at($k);
#    print "$point[$NZ] $y[$NZ] $x[$NZ]\n";
    $NZ++;
}



$weights=1/(0.01+abs(pdl(@flux_now)));

my $point_pdl=pdl(@point);
my $y_pdl = pdl(@y);
($s_y,$coeff) = fitpoly1d $point_pdl,$y_pdl,$npoly_y,{Weights => $weights};

print "Fitting...\n";


for ($j=0;$j<$nz;$j++) {
    $y_out[$j]=0;
    for ($i=0;$i<$npoly_y;$i++) {
	$C=$coeff->at($i);
	$y_out[$j]=$y_out[$j]+$C*($id[$j]**$i);
    }
}

$x_pdl = pdl(@x);
($s_x,$coeff) = fitpoly1d $point_pdl,$x_pdl,$npoly_x;

for ($j=0;$j<$nz;$j++) {
    $x_out[$j]=0;
    for ($i=0;$i<$npoly_x;$i++) {
	$C=$coeff->at($i);
	$x_out[$j]=$x_out[$j]+$C*($id[$j]**$i);
    }
}


($i_min,$i_max)=minmax(@x_out);
($j_min,$j_max)=minmax(@y_out);
if ($i_min<$j_min) {
    $min=$i_min;
} else {
    $min=$j_min;
}

if ($i_max>$j_max) {
    $max=$i_max;
} else {
    $max=$j_max;
}

$max=$max+0.5;
$min=$min-0.5;



pgbegin(0,"666/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
pgenv(1,$nz,$min-1,$max+1,0,0);
pgsci(2);
pgpoint($nz,\@point,\@x,3);
pgsci(4);
pgpoint($nz,\@point,\@y,3);
pgsci(1);
pgline($nz,\@id,\@x_out);
pgline($nz,\@id,\@y_out);
pgsci(3);
@F=list($flux);
pgline($nz,\@id,\@F);
pgsci(1);
pglabel("spec pix","Xc,Yc","$input_cube");
pgclose;
pgend;

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
pgenv(1,$nz,$min-1,$max+1,0,0);
pgsci(2);
pgpoint($nz,\@point,\@x,3);
pgsci(4);
pgpoint($nz,\@point,\@y,3);
pgsci(1);
pgline($nz,\@id,\@x_out);
pgline($nz,\@id,\@y_out);
pgsci(3);
@F=list($flux);
pgline($nz,\@id,\@F);
pgsci(1);
pglabel("spec pix","Xc,Yc","$input_cube");
pgclose;
pgend;


#
# We create the new cube
#

print "Shifting...\n";

if ($recenter==0) {
    $x0=mean(@x_out);
    $y0=mean(@y_out);
}

for ($k=0;$k<$nz-1;$k++) {
    $DX=($x0-$x_out[$k]);
    $DY=($y0-$y_out[$k]);
#    print "$DX $DY\n";
    $pdl_in_sec=$pdl_in->slice(",,($k)");
    $shift=pdl($DX,$DY);
    $tr=t_offset($shift);
    my $a=$pdl_in_sec->map($tr,{pix=>1});
    my $t = $pdl_out->slice(",,($k)");
    $t .= $a;
#    print "$t\n";
    

#    @dims=$t->dims;
}

# This was the only difference I found

$nz_med=int($nz/2);

for ($i=0;$i<$nx_0;$i++) {
    for ($j=0;$j<$ny_0;$j++) {
	$val=$pdl_out->at($i,$j,$nz_med);
	if (abs($val)>1e100) {
	    my $t = $pdl_out->slice("($i),($j),");
	    my $zero = zeroes($nz);
	    $t .= $zero;
	}
    }
}



#$H->{$NAXIS1}=$nx_0;
#$H->{$NAXIS2}=$ny_0;
#
$CDELT1=$H->{CDELT1};
$CDELT2=$H->{CDELT2};
$CDELT3=$H->{CDELT3};
$CRPIX1=$H->{CRPIX1};
$CRPIX2=$H->{CRPIX2};
$CRPIX3=$H->{CRPIX3};
$CRVAL1=$H->{CRVAL1};
$CRVAL2=$H->{CRVAL2};
$CRVAL3=$H->{CRVAL3};

#$h = {NAXIS=>3, NAXIS1=>$nx_0, NAXIS=>$ny_0, NAXIS3=>$nz, CDELT1=>$CDELT1,  CDELT2=>$CDELT2,  CDELT3=>$CDELT3,  CRPIX1=>$CRPIX1,  CRPIX2=>$CRPIX2,  CRPIX3=>$CRPIX3,  CRVAL1=>$CRVAL1,  CRVAL2=>$CRVAL2,  CRVAL3=>$CRVAL4,COMMENT=>"Sample FITS-style header"};
$pdl_out=$pdl_out/($scale*$scale);
$pdl_out->hdr->{CDELT1}=$CDELT1/$scale;
$pdl_out->hdr->{CDELT2}=$CDELT2/$scale;
$pdl_out->hdr->{CDELT3}=$CDELT3;
$pdl_out->hdr->{CRPIX1}=$CRPIX1;
$pdl_out->hdr->{CRPIX2}=$CRPIX2;
$pdl_out->hdr->{CRPIX3}=$CRPIX3;
$pdl_out->hdr->{CRVAL1}=$CRVAL1;
$pdl_out->hdr->{CRVAL2}=$CRVAL2;
$pdl_out->hdr->{CRVAL3}=$CRVAL3;


#$pdl_out->sethdr($h);


$pdl_out->wfits($output_cube);

exit;

