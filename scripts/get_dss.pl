#!/usr/bin/perl
#
#
#

#http://archive.stsci.edu/cgi-bin/dss_search?v=1&r=00+08+54.72&d=%2B23+49+02.2&e=J2000&h=3.0&w=3.0&f=gif&c=none&fov=NONE&v3=



use Astro::FITS::CFITSIO;
use PGPLOT;
use Carp;


if ($ARGV<0) {
    print "USE: get_dss.pl FILE( HH MM SS +DD MM SS)\n";
}
#
# open FITS file
#
my $status = 0;


$xsize=5.0;
$ysize=5.0;
#$equ="B1950";
$equ="J1200";
$file=$ARGV[0];
open(FH,"<$file");
while($line=<FH>) {
  chop($line);
  ($name,$hh,$mm,$ss,$dd,$dm,$ds,$junk)=split(" ",$line);
  $ask="http://archive.stsci.edu/cgi-bin/dss_search?v=1&r=".$hh."+".$mm."+".$ss."&d=".$dd."+".$dm."+".$ds."&e=".$equ."&h=".$ysize."&w=".$xsize."&f=fits&c=none&fov=NONE&v3=";
#  print "$ask\n";
  $name_gif=$name.".fits";
  $name_ps=$name.".ps";
  system("/usr/bin/wget '$ask' -O $name_gif");
  
  my $naxes;
  my $fptr = Astro::FITS::CFITSIO::open_file($name_gif,Astro::FITS::CFITSIO::READONLY(),$status);
  $fptr->get_img_parm(undef,undef,$naxes,$status);
  my ($naxis1,$naxis2) = @$naxes;

  print "Reading ${naxis2}x${naxis1} image...";
  my ($array, $nullarray, $anynull);
  $fptr->read_pixnull(Astro::FITS::CFITSIO::TLONG(), [1,1], $naxis1*$naxis2, $array, $nullarray, $anynull ,$status);
  print "done\n";
  
  $fptr->close_file($status);
  
  $file_out=$name_ps."/CPS";
  @tr=(-$xsize/2,$xsize/$naxis2,0,-$ysize/2,0,$ysize/$naxis1);
  pgbeg(0,'/xs',1,1);
#  pgenv(-1.5,1.5,-1.5,1.5,1,0);
  pgenv(-$xsize/2,$xsize/2,-$ysize/2,$ysize/2,1,0);
  pgimag($array,$naxis1,$naxis2,1,$naxis1,1,$naxis2,0,10000,\@tr);
  pglabel("\\gDRA (\')","\\gDDEC (\')",$name);
  pgclos();
  pgend();

  pgbeg(0,$file_out,1,1);
#  pgenv(-1.5,1.5,-1.5,1.5,1,0);
  pgenv(-$xsize/2,$xsize/2,-$ysize/2,$ysize/2,1,0);
  for ($i=0;$i<20;$i++) {
    $flux[$i]=4500+10**(0.2*$i)*200;
  }
#  pgcons($array,$naxis1,$naxis2,1,$naxis1,1,$naxis2,\@flux,$i,[-1.5,3/$naxis2,0,-1.5,0,3/$naxis1]);
#  pgcons($array,$naxis1,$naxis2,1,$naxis1,1,$naxis2,\@flux,$i,[-1.5,3/$naxis2,0,-1.5,0,3/$naxis1]);
  pggray($array,$naxis1,$naxis2,1,$naxis1,1,$naxis2,10000,0,\@tr);
  pglabel("\\gDRA (\')","\\gDDEC (\')",$name);
  pgclos();
  pgend();




}
close(FH);

exit;
