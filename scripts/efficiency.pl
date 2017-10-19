#!/usr/bin/perl

#
# Read FITS image, unpacking data into Perl array.
# Display image with PGPLOT
#

use PGPLOT;
use Carp;
use Statistics::OLS;
use Math::Stat;

use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

Astro::FITS::CFITSIO::PerlyUnpacking(0);

if ($#ARGV<4) {
    print "USE: efficiency.pl REAL_FLUX.SPEC OBSERVED_FLUX.SPEC N_SPAXELS NB_SPEC OUTFILE\n";
    exit;
}

$file_real=$ARGV[0];
$file_obs=$ARGV[1];
$n_spax=$ARGV[2];
$nb_spec=$ARGV[3];
$outfile=$ARGV[4];

$nr=0;
open(FH,"<$file_real");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$wr[$nr]=$data[0];
	$fr[$nr]=$data[1];
	$nr++;
    }
}
close(FH);
print "#$nr data on file $file_real\n";


$spline=new Math::Spline(\@wr,\@fr);
@y2=Derivative2(\@wr,\@fr);



$rat_max=-10000000000;
$rat_min=10000000000;



$mean=0;
$no=0;
$n_mean=0;
open(FH,"<$file_obs");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$wo[$no]=$data[1];
	$fo[$no]=$data[2]*$n_spax;
#
# We interpolate the real value for the
# proper wavelength
#
#	$fr_wo[$no]=$spline->evaluate($wo[$no]);


	if ($wo[$no]<$wr[$nr-1]) {
	    $index=binsearch(\@wr,$wo[$no]);
	    $index=linsearch(\@wr,$wo[$no],$index);
	    $fr_wo[$no]=spline(\@wr,\@fr,\@y2,$index,$wo[$no]);
	    $fr_last=$fr_wo[$no];
	} else {
	    $fr_wo[$no]=$fr_last;
	}

#
# We determine the ratio
#
#	print "$wo[$no] $fo[$no] $fr_wo[$no] $ratio[$no]\n";
	if ($fo[$no]) {
	    $ratio[$no]=$fr_wo[$no]/$fo[$no];
	} else {
	    $ratio[$no]=0;
	}
	if ($rat_min>$ratio[$no]) {
	    $rat_min=$ratio[$no];
	}
	if ($rat_max<$ratio[$no]) {
	    $rat_max=$ratio[$no];
	}
	if (($wo[$no]>3800)&&($wo[$no]<9000)) {
	    $val[$n_mean]=$ratio[$no];
	    $n_mean++;
	}
	$no++;
    }
}
close(FH);
print "#$no data on file $file_obs\n";
$stat = Math::Stat->new(\@val);
$mean = $stat->average();
$sig = $stat->stddev();
#$mean=$mean/$n_mean;
#$sig=0.5*$mean;






$rat_min=$mean-5*$sig;
if ($rat_min<0) {
    $rat_min=0;
}
$rat_max=$mean+5*$sig;
print "#$mean $rat_min $rat_max\n";



#$a = new Math::Approx (\&poly, 10, %y);
#$a->print;
#$a->plot("math-approx-demo.out");
#print "Fit: ", $a->fit, "\n";      




pgbegin(0,$dev,1,1);  # Open plot device
pgscf(2.0);             # Set character font
pgslw(1.2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wo[0],$wo[$no-1],$rat_min,$rat_max,0,0);
pglabel("Wavelength (AA)","F\\dreal/\\uF\\dobs",$file_obs);
#pgsci(1);
pgpoint($no,\@wo,\@ratio,1);
#pgpoint($no,\@wo,\@ratio,1);
#pgsci(2);
#


#
# Cleaning process
#
#@y2_new=Derivative2(\@wo,\@ratio);

$smooth=30;

$n=0;
for ($k=$smooth;$k<$no-$smooth;$k=$k+$smooth) {
    my @val;
    my @w_val;
    my @val2;
    for ($j=0;$j<2*$smooth;$j++) {
	$val[$j]=$ratio[$k+$j-$smooth];
#	$w_val[$j]=$wo[$k+$j-$smooth];
    }
    $stat = Math::Stat->new(\@val);
    $mean = $stat->median();
    $sig = $stat->stddev();
    $nn=0;
    for ($j=0;$j<2*$smooth;$j++) {	
	$val2[$nn]=$ratio[$k+$j-$smooth];
	$w_val[$nn]=$wo[$k+$j-$smooth];
	if (abs($val2[$nn]-$mean)<$sig) {
	    $nn++;
	}
    }
#    if ($nn>2) {    
    $stat = Math::Stat->new(\@val2);
    $y[$n] = $stat->median();
    $stat = Math::Stat->new(\@w_val);
    $x[$n] = $stat->median();
    $n++;
#    }
#    $sig = $stat->stddev();
}

$spline=new Math::Spline(\@x,\@y);
@y2=Derivative2(\@x,\@y);
for ($j=0;$j<$no;$j++){
    $index=binsearch(\@x,$wo[$j]);
    $index=linsearch(\@x,$wo[$j],$index);
    $sp_ratio[$j]=spline(\@x,\@y,\@y2,$index,$wo[$j]);
#    $sp_ratio[$j]=median_boxspline(\@x,\@y,\@y2,$index,$wo[$j]);

}

@sp_ratio2=median_box(3,\@ratio);


#for ($k=0;$k<$smooth;$k++) {
#    $sp_ratio[$k]=$ratio[$k];
#    $val_x[$j]=$wo[$k];
#}

#for ($k=$smooth;$k<$no;$k++) {
#    my @val;
#    my @val2;
#    my @val3;
#    for ($j=0;$j<2*$smooth;$j++) {
#	$val[$j]=$ratio[$k+$j-$smooth];
#    }
#    $stat = Math::Stat->new(\@val);
#    $mean = $stat->median();
#    $sig = $stat->stddev();
#    $n=0;
#    for ($j=0;$j<2*$smooth;$j++) {	
#	$val2[$n]=$ratio[$k+$j-$smooth];
#	if (abs($val2[$n]-$mean)<$sig) {
#	    $n++;
#	}
#    }
#    $stat = Math::Stat->new(\@val2);
#    $mean = $stat->median();
#    $sig = $stat->stddev();
#    $n=0;
#    for ($j=0;$j<2*$smooth;$j++) {	
#	$val3[$n]=$ratio[$k+$j-$smooth];
#	if (abs($val3[$n]-$mean)<0.5*$sig) {
#	    $n++;
#	}
#    }
#    $stat = Math::Stat->new(\@val3);
#    $mean = $stat->median();
#    $sig = $stat->stddev();
#
#    $sp_ratio[$k]=$mean;
#
#spline(\@wo,\@ratio,\@y2_new,$index,$wo[$k]);
#}
#@y2_new=Derivative2(\@wo,\@ratio);
#for ($k=0;$k<$no;$k++) {
#    $index=binsearch(\@wo,$wo[$k]);
#    $index=linsearch(\@wo,$wo[$k],$index);
#    $sp_ratio[$k]=spline(\@wo,\@ratio,\@y2_new,$index,$wo[$k]);
#}

pgsci(2);
pgline($no,\@wo,\@sp_ratio);
pgsci(3);
#pgline($no,\@wo,\@sp_ratio2);
pgclos();
pgend();



#
# Now we write the output
#
system("rm $outfile");
fits_clear_errmsg();
$fptr=Astro::FITS::CFITSIO::create_file($outfile,$status);
print "ffinit create new file ($outfile) status = $status\n";
$simple=1;
$bitpix=-32;
$naxis=2;
$naxes=[$no,$nb_spec];
#$npixels=20;
#$pcount=0;
#$gcount=1;
#$extend=1;

############################
#  write single keywords   #
############################

$fptr->write_imghdr($bitpix,$naxis,$naxes,$status);

#print "ffinit create new file ($outfile) status = $status\n";

for ($j=0;$j<$nb_spec;$j++) {
    for ($i=0;$i<$no;$i++) {
#	$float_array[$j][$i] = $sp_ratio[$i];
	$float_array[$j][$i] = $ratio[$i];
    }
}
#    $fptr->write_img(TFLOAT, $j, $no, $float_pdl->slice('0:10')->get_dataref, $status);
$fptr->write_img(TFLOAT, 1, $no*$nb_spec+1, \@float_array, $status);
#    $fptr->write_img(TFLOAT, $j, $no,\@wo, $status);
print "STATUS=$status\n";

$fptr->close_file($status);
print "STATUS=$status\n";
exit;

sub poly {
    my($n,$x) = @_;
    return $x ** $n;
}



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
  $sum=$sum/($#data-1);
  $stddev=sqrt($sum);
  return $stddev;
}


sub median_box {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    @in_val=@$array;
    my $i,$j,$k;
    $k=0;
    for ($i=$box;$i<($#in_val+1-$box);$i=$i+2*$box) {
	my @tmp;
	for ($j=0;$j<2*$box;$j++) {
	    $tmp[$j]=@$array[$i-$box+$j];
	}
	$out_val[$k]=median(@tmp);
#	$out_sval[$k]=sigma(@tmp);
	$k++;
    }
    return @out_val;
}

sub median {
  my @data=@_;
# sort numerically ascending
  my @out_data = sort {$a <=> $b} @data;
  my $n_median=int(($#out_data+1)*0.5);
  my $median = ($out_data[$n_median-1]+$out_data[$n_median]+$out_data[$n_median+1])/3;
  return $median;
}

