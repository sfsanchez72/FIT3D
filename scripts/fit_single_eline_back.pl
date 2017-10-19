#!/usr/bin/perl
#
#
#

#use PGPLOT;

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");


if ($#ARGV<12) {
    print "USE:fit_single_eline_back.pl SPECTRUM.TXT START_WAVELENGTH END_WAVELENGTH BACK_LIST NBACK MASK_LIST WHICH_BACK[=-1 ALL] REDSHIFT VEL_DISP FILE_OUT ELINE_CONFIG VISUAL=[0/1] AV(DUST)\n";    
    exit;
}

$spec_file=$ARGV[0];
$w_start=$ARGV[1];
$w_end=$ARGV[2];
$back_list=$ARGV[3];
$nback=$ARGV[4];
$mask_list=$ARGV[5];
$which_back=$ARGV[6];
$redshift=$ARGV[7];
$vel_disp=$ARGV[8];
$file_out=$ARGV[9];
$eline_config=$ARGV[10];
$visual=$ARGV[11];
$AV=$ARGV[12];
$no_config=0;

$nf=0;
open(FH,"<$back_list");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$file[$nf]=$line;
	$nf++;
    }
}
close(FH);

for ($j=0;$j<$nf;$j++) {
    open(FH,"<$file[$j]");
    $nx=0;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($j==0) {
	    if ($nx==0) {
		$crval=$data[1];
	    }
	    if ($nx==1) {
		$cdelt=$data[1]-$crval;
	    }
	}
	$id[$nx]=$data[0];
	$w[$nx]=$data[1];
	$f[$nx]=$data[2];
	$nx++;
    }
    close(FH);
    if ($j==0) {
	$pdl_out=zeroes($nx,$nf);
    }
    for ($i=0;$i<$nx;$i++) {
	set($pdl_out,$i,$j,$f[$i]);
    }
}
$$h{"CRPIX1"}=1;
$$h{"CRVAL1"}=$crval;
$$h{"CDELT1"}=$cdelt;

$pdl_out->sethdr( $h );
$pdl_out->wfits("BACK_JUNK.fits");

$ne_lines=0;
$ne=0;
open(ELINE,"<$eline_config");
$line=<ELINE>;
chop($line);
@data=split(" ",$line);
$ne=$data[1];
for ($i=0;$i<$ne;$i++) {
    for ($j=0;$j<10;$j++) {
	$line=<ELINE>;
	chop($line);
	$eline[$i][$j]=$line;
    }
}
close(ELINE);



if ($mask_list eq "none") {
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

#print "PASO\n";

$w_med=($w_start+$w_end)/2;
print "$spec_file $w_start $w_end\n";
$nw=0;
$y_min=10000;
$y_max=-10000;
$min_dist=1e12;
open(FH,"<$spec_file");
open(FHOUT,">fit_spectra.input");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id=$data[0];
    $w_test=$data[1];
    $f_test=$data[2];
    
#    print "$w_spec $w_start $w_end\n";
    if (($w_test>$w_start)&&($w_test<$w_end)) {
	$masked=0;
	for ($j=0;$j<$nmask;$j++) {
	    if (($w_test>$start_mask[$j])&&($w_test<$end_mask[$j])) {
		$masked=1;
		
	    }
	}
	if ($f_test<=0) {
	    $masked=1;
	}
	if ($masked==0) {
	    print FHOUT "$id $w_test $f_test\n";
	    if (abs($w_test-$w_med)<$min_dist) {
		$min_dist=abs($w_test-$w_med);
		$f_med=$f_test;
	    }
	    $w_spec[$nw]=$data[1];
	    $f_spec[$nw]=$data[2];
	    if ($y_min>$f_spec[$nw]) {
		$y_min=$f_spec[$nw];
	    }
	    if ($y_max<$f_spec[$nw]) {
	    $y_max=$f_spec[$nw];
	}
	    $nw++;
	}
    }
}
close(FHOUT);
close(FH);
#print "F_MED=$f_med\n";
#exit;
if ($which_back>-1) {
    $nnback=0;
    $min_dist=1e12;
    open(FH,"<$back_list");
    while($line=<FH>) {
	if ($nnback==$which_back) {
	    print "FILE_BACK = '$line'\n";
	    open(SPEC,"<$line");
	    while($line_spec=<SPEC>) {
		chop($line_spec);
		@data=split(" ",$line_spec);
		$id_back=$data[0];
		$w_back=$data[1]*(1+$redshift);
		$f_back=$data[2];	
#		print "$line\n";
#		print "$which_back $nnback $w_back $f_back $f_med $ratio\n";
		if (abs($w_back-$w_med)<$min_dist) {
		    $min_dist=abs($w_back-$w_med);
		    $ratio=$f_med/$f_back;
		}
	    }
	    close(SPEC);
	}
	$nnback++;
    }
    close(FH);

    $vel=$redshift*299792.458;
    print "Systemic Velocity=$vel km/s\n";
    $flex_w=0;
    $dz=$vel_disp;    
    open(FH,">tmp.config");
    $nmod=$ne+1;
    if ($dz>0) {
	print FH "0 $nmod 0.5 0.1\n";
    } else {
	print FH "0 $nmod 0.002 0.00001\n";
    }
    for ($i=0;$i<$ne;$i++) {
	for ($j=0;$j<10;$j++) {
	    print FH "$eline[$i][$j]\n";
	}
    }

    print FH "back\n";
    print FH "$which_back\t 0\t 0\t 0\t -1\n"; 
    if ($dz>0) {
	print FH "$ratio\t 1\t 0.001\t 1e12\t -1\n";
    } else {
	print FH "$ratio\t 1\t 0.001\t 1e12\t -1\n";
    }
    if ($dz>0) {
	print FH "$dz\t  0\t 0\t 0\t -1\n";
    } else {
	$dza=abs($dz);
	$dza_min=0.2*$dza;
	$dza_max=2.0*$dza;
	print FH "$dza\t  1\t $dza_min\t $dza_max\t -1\n";
    }
    if ($vel>0) {
	print FH "$vel\t 0\t 0\t 0\t -1\n"; 
    } else {
	$vela=abs($vel);
	$vela_min=0.2*$vela;
	$vela_max=2.0*$vela;
	print FH "$vela\t  1\t $vela_min\t $vela_max\t -1\n";
    }
    print FH "$AV\t 1\t 0.0\t 50\t -1\n";
    for ($k=0;$k<5;$k++) {
	print FH "0\t 0\t 0\t 0\t -1\n";
    }
    close(FH);
    
    system("rm out_config.fit_spectra");
    system("rm out.fit_spectra");
#    $call="fit_spec_back fit_spectra.input tmp.config ".$back_list." ".$nback;
    $call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback;
    system($call);
    print "$call\n";
    if ($visual==1) {
	$call="plot_out_fit.pl out_mod_res.fit_spectra /xs";
	system($call);
    }

    $n=0;
    $chisq=0;
    open(FH,"<out_mod_res.fit_spectra");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$res[$n]=$data[3];
	if ($data[1]!=0) {
	    $chisq=$chisq+(($res[$n])**2)/abs($data[1]);
	}
	$n++;
    }
    $chisq=$chisq/($n-1);
    close(FH);
    print "CHISQ $chisq\n";
} else {
    print "PASO\n";
    $nn=0;
    open(LIST,"<$back_list");
    while($line=<LIST>) {
	chop($line);
	@data=split(/\//,$line);
	$all_name[$nn]=$line;
	$name[$nn]=$data[$#data];
	$nn++;
    }
    close(FH);

    close(LIST);
    open(OUTFIT,">$file_out");
    for ($i=0;$i<$nback;$i++) {
	$min_dist=1e12;
	print "FILE_BACK = '$name[$i]'\n";
	open(SPEC,"<$all_name[$i]");
	while($line_spec=<SPEC>) {
	    chop($line_spec);
	    @data=split(" ",$line_spec);
	    $id_back=$data[0];
	    $w_back=$data[1]*(1+$redshift);
	    $f_back=$data[2];	    
#	    print "$line_spec\n";
#	    print "$min_dist $w_back $w_med $f_med $ratio\n";
	    if (abs($w_back-$w_med)<$min_dist) {
		$min_dist=abs($w_back-$w_med);
		$ratio=$f_med/$f_back;
	    }
	}
	close(SPEC);
#	<stdin>;




	$vel=$redshift*299792.458;
	print "Systemic Velocity=$vel km/s\n";
	$flex_w=0;
	$dz=$vel_disp;
	open(FH,">tmp.config");
	$nmod=$ne+1;
	if ($dz>0) {
	    print FH "0 $nmod 0.5 0.1\n";
	} else {
	    print FH "0 $nmod 0.002 0.00001\n";
	}
	for ($ii=0;$ii<$ne;$ii++) {
	    for ($j=0;$j<10;$j++) {
		print FH "$eline[$ii][$j]\n";
	    }
	}
	print FH "back\n";
	print FH "$i\t 0\t 0\t 0\t -1\n"; 
	if ($dz>0) {
	    print FH "$ratio\t 1\t 0.001\t 1e12\t -1\n";
	} else {
	    print FH "$ratio\t 1\t 0.001\t 1e12\t -1\n";
	}
	if ($dz>0) {
	    print FH "$dz\t  0\t 0\t 0\t -1\n";
	} else {
	    $dza=abs($dz);
	    $dza_min=0.2*$dza;
	    $dza_max=2.0*$dza;
	    print FH "$dza\t  1\t $dza_min\t $dza_max\t -1\n";
	}
	if ($vel>0) {
	    print FH "$vel\t 0\t 0\t 0\t -1\n"; 
	} else {
	    $vela=abs($vel);
	    $vela_min=0.2*$vela;
	    $vela_max=2.0*$vela;
	    print FH "$vela\t  1\t $vela_min\t $vela_max\t -1\n";
	}
	print FH "$AV\t 1\t 0\t 50\t -1\n";
	for ($k=0;$k<4;$k++) {
	    print FH "0\t 0\t 0\t 0\t -1\n";
	}
	close(FH);
	
	system("rm out_config.fit_spectra");
	system("rm out.fit_spectra");
#	$call="fit_spec_back fit_spectra.input tmp.config ".$back_list." ".$nback;
	$call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback;
	system($call);
	print "$call\n";
	if ($visual==1) {
	    $call="plot_out_fit.pl out_mod_res.fit_spectra /xs";
	    system($call);
	}
	$n=0;
	$chisq=0;
	open(FH,"<out_mod_res.fit_spectra");
	while($line=<FH>) {
	    chop($line);
	    @data=split(" ",$line);
	    $res[$n]=$data[3];
	    if ($data[1]!=0) {
		$chisq=$chisq+(($res[$n])**2)/abs($data[1]);
	    }
	    $n++;
	}
	$chisq=$chisq/($n-1);
	close(FH);
	open(FH,"<out.fit_spectra");
	$tmp=<FH>;
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$flux=$data[3];
	$eflux=$data[4];
	$disp=$data[5];
	$edisp=$data[6];
	$vel=$data[7];
	$evel=$data[8];
	close(FH);

	print "$name[$i] $chisq \n";
	print OUTFIT "$name[$i] $chisq $flux $eflux $disp $edisp $vel $evel\n";		
    }
    close(OUTFIT);
}




exit;

sub mean {
  local(@data)=@_;
  my $sum=0; 
  my $j;
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+$data[$j];
  }
  my $mean;
  if ($#data>0) {
     $mean = $sum/($#data);
  } else {
     $mean=-666;
  }
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
  if ($#data>0) {
      if ($#data>1) {
	  $sum=$sum/($#data-1);
      } else {
	  $sum=$sum/($#data);
      } 
      $stddev=sqrt($sum);
  } else {
      $stddev=-666;
  }

  return $stddev;
}

