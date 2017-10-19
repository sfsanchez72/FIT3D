#!/usr/bin/perl
#
# S.F.Sanchez March 2003
# Modified on 31.07.03
#
# This program read an input_v.txt (read README.txt) file
# and run GALFIT over the data to extract the information
#



use DBI;
use PGPLOT;  # Load PGPLOT module
use Statistics::OLS;
use Math::Stat;

$galfit="/work1/ssa/galfit/galfit";

print "This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s a 'input_v.txt' file\n";
print "containing the files to fit.\n";
print "It also need a template for the galfit (default=galfit.template)\n";
print "It also need the 'flag' to identify the output-files:\n";
print "e.g.: if 'vauc' selected, then for an input file named\n";
print "qso_z05A_00784.fits\n";
print "the output will be:\n";
print "qso_z05A_00784_gal_vauc.fits\n";

if ($#ARGV<3) {
    print "USE:\n";
    print "analyze_3c120.pl inputfile galfit.template flag FIX\n";
    exit;
}
$input_file=$ARGV[0];
$galfit_template=$ARGV[1];
$flag="GALFIT_".$ARGV[2].".fits";
$flag_nuc="GALFIT_NUC_".$ARGV[2].".fits";
$fix=$ARGV[3];

#
# We open the directory to read the information on it
#

$n_files=0;
open(FH,"<$input_file");
while($line=<FH>) {
  chop($line);
  if ($line !~ "#") {
      @data=split(" ",$line);
      $files[$n_files]=$data[0];
      $psf_in[$n_files]=$data[1];
      $x_cen_in[$n_files]=$data[2];
      $y_cen_in[$n_files]=$data[3];
      $magzero_in[$n_files]=$data[4];
      $mag_in[$n_files]=$data[5];
      $r_in[$n_files]=$data[6];
      $z_in[$n_files]=$data[7];
      $R_in[$n_files]=$data[8];
      $ab_in[$n_files]=$data[9];
      $pang_in[$n_files]=$data[10];
      $n_files++;
  }
}
close(FH);
print "$n_files to fit\n";



$out_file="out_".$ARGV[2]."_".$input_file;
open(FHOUT,">$out_file");
print FHOUT "# 1  filename\n";
print FOOUT "# 2  x_c\n";
print FHOUT "# 3  y_c\n";
print FHOUT "# 4  z \n";
print FHOUT "# 5  R \n";
print FHOUT "# 6  mag_out\n";
print FHOUT "# 7  mag_nuc_out\n";
print FHOUT "# 8  mag_hg_out\n";
print FHOUT "# 9 r_hn_out\n";
print FHOUT "# 10 r_out\n";
print FHOUT "# 11 A/B\n";
print FHOUT "# 12 PA\n";
print FHOUT "# 13 Chisq\n";
print FHOUT "# 14 N-Sersic\n";
close(FHOUT);

$que_file=0;
for ($i=0;$i<$n_files;$i++) {
    $file=$files[$i];
  ###############################################################
  #We read the GALFIT template
  ###############################################################
  $out=$file;
  $out =~ s/fits/$flag/;
  $out_nuc=$file;
  $out_nuc =~ s/fits/$flag_nuc/;
  $var_file=$file;
  $var_file =~ s/.fits/_s.fits/;
# OLD modified on 20.08.03
#  $cnx=$x_cen_in[$i]+0.53;
#  $cny=$y_cen_in[$i]+0.53;
  $cnx=$x_cen_in[$i]+0.5;
  $cny=$y_cen_in[$i]+0.5;
  $mag=$mag_in[$i]-0.5;
  $mag_d=$mag-0.5;
  $mag_e=$mag-0.5;
  $r_e=$r_in[$i];
  $r_d=$r_in[$i];
  open(FH,"<$galfit_template");
    
#    if( !(-e "$dir/ID.txt") ) {


  open(FHOUT,">galfit.tmp");
  while($line=<FH>) {
    chop($line);
    $line =~ s/INPUTFILE/$file/g;
    $line =~ s/OUTPUTFILE/$out/g;
    $line =~ s/NOISEFILE/$var_file/g;
    $line =~ s/PSFFILE/$psf_in[$i]/g;
    $line =~ s/CNX/$cnx/g;
    $line =~ s/CNY/$cny/g;
    $line =~ s/NX/51/g;
    $line =~ s/NY/51/g;
    $line =~ s/COX/51/g;
    $line =~ s/COY/51/g;
    $line =~ s/MAGZERO/$magzero_in[$i]/g;
    $line =~ s/SCALE/0.03/g;
    $line =~ s/CAG_D/$mag_d/g;
    $line =~ s/CAG_E/$mag_e/g;
    $line =~ s/R_D/$r_d/g;
    $line =~ s/R_E/$r_e/g;
    $line =~ s/MAG/$mag/g;
    $line =~ s/NSERSIC/2.5/g;
    $line =~ s/FIX/$fix/g;
    $line =~ s/AB/$ab_in[$i]/g;
    $line =~ s/PANG/$pang_in[$i]/g;
    print FHOUT "$line\n";
  }
  close(FHOUT);
  close(FH);



#    exit;
  ###############################################################
  #We use galfit to analyze this simulation
  ###############################################################
  $command=$galfit." galfit.tmp";
  system($command);


  ###############################################################
  # Output Parameters
  ###############################################################
  open(FH,"<fit.log");
  $error="none";
    $type="none";
  while($line=<FH>) {
    chop($line);

    
    ($p1,$p2)=split(/\)/,$line);
    ($p3,$p4)=split(/\(/,$p1);
    ($x_t,$y_t)=split(",",$p4);

    if ($error eq "gaussian") {
      $portion=substr($line,29);
      ($e_mag_nuc_out,$e_r_nuc_out,$e_ab_nuc_out,$e_pa_nuc_out)=split(" ",$portion);
      ($e_x_nuc,$e_y_nuc)=($x_t,$y_t);
      $e_nser_out=0.0;
      $type=$error;
      $error="none";
    }
    if ($error eq "expdisk") {
      $portion=substr($line,29);
      ($e_mag_d_out,$e_r_d_out,$e_ab_d_out,$e_pa_d_out)=split(" ",$portion);
      ($e_x_d,$e_y_d)=($x_t,$y_t);
      $e_nser_out=0.0;
      $type=$error;
      $error="none";
    }
    if ($error eq "devauc") {
      $portion=substr($line,29);
      ($e_mag_e_out,$e_r_e_out,$e_ab_e_out,$e_pa_e_out)=split(" ",$portion);
      ($e_x_e,$e_y_e)=($x_t,$y_t);
      $e_nser_out=0.0;
      $type=$error;
      $error="none";
    }
    if ($error eq "sersic") {
#	print "PASO Error!!!";
      $portion=substr($line,29);
      ($e_mag_s_out,$e_r_s_out,$e_nser_out,$e_ab_s_out,$e_pa_s_out)=split(" ",$portion);
      ($e_x_s,$e_y_s)=($x_t,$y_t);
      $type=$error;
      $error="none";
    }
    if ($line =~ "gaussian") {
      $portion=substr($line,29);
      ($mag_nuc_out,$r_nuc_out,$ab_nuc_out,$pa_nuc_out)=split(" ",$portion);
      ($x_nuc,$y_nuc)=($x_t,$y_t);
      $nser_out=0;
      $error="gaussian";
    }
    if ($line =~ "expdisk") {
      $portion=substr($line,29);
      ($mag_d_out,$r_d_out,$ab_d_out,$pa_d_out)=split(" ",$portion);
      ($x_d,$y_d)=($x_t,$y_t);
      $nser_out=1;
      $error="expdisk";
    }
    if ($line =~ "devauc") {
      $portion=substr($line,29);
      ($mag_e_out,$r_e_out,$ab_e_out,$pa_e_out)=split(" ",$portion);
      ($x_e,$y_e)=($x_t,$y_t);
      $nser_out=4;
      $error="devauc";
    }
    if ($line =~ "sersic") {
#	print "PASO!!!";
      $portion=substr($line,29);
      ($mag_s_out,$r_s_out,$nser_out,$ab_s_out,$pa_s_out)=split(" ",$portion);
      ($x_s,$y_s)=($x_t,$y_t);
      $error="sersic";
    }


#    print "$line\n $type|$error\n";

    if ($line =~ "Chi") {
       ($code,$chi_sq)=split("=",$line);
    }
  }
  close(FH);
#    exit;
  system("/bin/rm fit.log");
  system("/bin/rm galfit.??");

    print "TIPO=$type\n";
    if ($type eq "expdisk") {
	$mag_hg_out=$mag_d_out;
	$e_mag_hg_out=$e_mag_d_out;
	$r_hg_out=$r_d_out;
	$e_r_hg_out=$e_r_d_out;
	$ab_hg_out=$ab_d_out;
	$e_ab_hg_out=$e_ab_d_out;
	$pa_hg_out=$pa_d_out;
	$e_pa_hg_out=$e_pa_d_out;
    } 
    if ($type eq "devauc") {
	$mag_hg_out=$mag_e_out;
	$e_mag_hg_out=$e_mag_e_out;
	$r_hg_out=$r_e_out;
	$e_r_hg_out=$e_r_e_out;
	$ab_hg_out=$ab_e_out;
	$e_ab_hg_out=$e_ab_e_out;
	$pa_hg_out=$pa_e_out;
	$e_pa_hg_out=$e_pa_e_out;
    }
    if ($type eq "sersic") {
	$mag_hg_out=$mag_s_out;
	$e_mag_hg_out=$e_mag_s_out;
	$r_hg_out=$r_s_out;
	$e_r_hg_out=$e_r_s_out;
	$ab_hg_out=$ab_s_out;
	$e_ab_hg_out=$e_ab_s_out;
	$pa_hg_out=$pa_s_out;
	$e_pa_hg_out=$e_pa_s_out;
    }
    if ($type eq "gaussian") {
	$mag_hg_out=99.99;
	$e_mag_hg_out=99.99;
	$r_hg_out=99.99;
	$e_r_hg_out=99.99;
	$ab_hg_out=99.99;
	$e_ab_hg_out=99.99;
	$pa_hg_out=99.99;
	$e_pa_hg_out=99.99;
    }
  $mag_out=mag_suma($mag_nuc_out,$mag_hg_out);
  $e_mag_out=sqrt($e_mag_nuc_out**2+$e_mag_hg_out**2);
  $r_hn_out=ratio_mag($mag_nuc_out,$mag_hg_out);
  $r_bd_out=0.0;
  $mag_out=apr($mag_out);
  $e_mag_out=apr($e_mag_out);
  $r_hn_out=apr($r_hn_out);
  $r_hg_out=$r_hg_out;
  $e_r_hg_out=apr($e_r_hg_out);

  print "MAG = $mag_in $mag_out \n";
  print "MAG_NUC = $mag_nuc_in $mag_nuc_out \n";
  print "MAG_HG = $mag_hg_in $mag_hg_out \n";
  print "MAG_D = $mag_d_in $mag_d_out \n";
  print "R_HN = $r_hn_in $r_hn_out \n";
  print "R_D = $r_d_in $r_d_out \n";
  print "Chisq/Nu = $chi_sq \n";




  #
  # We save the results!
  #
  
    $x_nuc=$x_nuc-0.53;
    $y_nuc=$y_nuc-0.53;

  print "OUT=$out_file\n";
  open(FHOUT,">>$out_file");
  print FHOUT "$file|$x_nuc|$y_nuc|$z_in[$i]|$R_in[$i]|$mag_out|$mag_nuc_out|$mag_hg_out|$r_hn_out|$r_hg_out|$ab_hg_out|$pa_hg_out|$chi_sq|$nser_out\n";
  print FHOUT "#ERRORS|$e_x_nuc|$e_y_nuc|-|-|$e_mag_out|$e_mag_nuc_out|$e_mag_hg_out|$e_r_hn_out|$e_r_hg_out|$e_ab_hg_out|$e_pa_hg_out|----|$e_nser_out\n";
  close(FHOUT);

#    print "WAITING....\n";
#    <stdin>;

#Only the first!
#  exit;

  
#    print "MAG=$mag\n";
#    exit;
#  $que_file++;
  print "******************************\n";
  print "******************************\n";
  print "*******   FIT  $que_file *****\n";
  print "*******     DONE         *****\n";
  print "******************************\n";
  print "******************************\n";
    print "\n\nWe create the nucleus model\n";
  open(FH,"<galfit.nucleus");
  open(FHOUT,">galfit.tmp");
  while($line=<FH>) {
    chop($line);
    $line =~ s/INPUTFILE/$file/g;
    $line =~ s/OUTPUTFILE/$out_nuc/g;
    $line =~ s/NOISEFILE/$var_file/g;
    $line =~ s/PSFFILE/$psf_in[$i]/g;
    $line =~ s/CNX/$cnx/g;
    $line =~ s/CNY/$cny/g;
    $line =~ s/NX/51/g;
    $line =~ s/NY/51/g;
    $line =~ s/COX/51/g;
    $line =~ s/COY/51/g;
    $line =~ s/MAGZERO/$magzero_in[$i]/g;
    $line =~ s/SCALE/0.03/g;
    $line =~ s/MAG/$mag_nuc_out/g;
    $line =~ s/FWHM/$r_nuc_out/g;
    print FHOUT "$line\n";
  }
  close(FHOUT);
  close(FH);
  ###############################################################
  #We use galfit to analyze this simulation
  ###############################################################
  $command=$galfit." galfit.tmp";
  system($command);

#  system("/bin/rm fit.log");
#  system("/bin/rm galfit.??");

}

exit;


sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub mag_suma {
  local($mag1,$mag2)=@_;
  my $mag=-2.5*log10(10**(-0.4*$mag1)+10**(-0.4*$mag2));
  return $mag;
}

sub ratio_mag {
  local($mag1,$mag2)=@_;
  my $ratio=10**(-0.4*($mag2-$mag1));
  return $ratio;
}

sub apr {
    my $z=@_[0];
    my @s=split(/\./,$z);
    my $last=substr($s[1],0,-length($s[1])+2);
    my $zz=$s[0].".".$last;
    return $zz;
}

