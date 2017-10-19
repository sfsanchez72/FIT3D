#!/usr/bin/perl
#
#
#

#http://archive.stsci.edu/cgi-bin/dss_search?v=1&r=00+08+54.72&d=%2B23+49+02.2&e=J2000&h=3.0&w=3.0&f=gif&c=none&fov=NONE&v3=



use Astro::FITS::CFITSIO;
use PGPLOT;
use Carp;
#
# open FITS file
#
my $status = 0;


$xsize=1.5*60;
$ysize=1.5*60;
#$equ="B1950";
$equ="J2000";
if ($#ARGV<1) {
    print "mos_pattern.pl Mos_config DEVICE\n";
    exit;
}
$dither_file=$ARGV[0];
$file_out=$ARGV[1];
$ratio=$ARGV[2];

open(FH,"<$dither_file");
$nd=0;
$xd[$nd]=0;
$yd[$nd]=0;
$nd++;
$xmin=0;
$xmax=0;
$ymin=0;
$ymax=0;
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $xd[$nd]=$xd[$nd-1]+$data[0];
    $yd[$nd]=$yd[$nd-1]+$data[1];
    if ($xd[$nd]<$xmin) {
	$xmin=$xd[$nd];
    }
    if ($yd[$nd]<$ymin) {
	$ymin=$yd[$nd];
    }
    if ($xd[$nd]>$xmax) {
	$xmax=$xd[$nd];
    }
    if ($yd[$nd]>$ymax) {
	$ymax=$xd[$nd];
    }
    $nd++;
}
close(FH);
$xmin=$xmin-35;
$ymin=$ymin-35;
$ymax=$ymax+35;
$xmax=$xmax+35;
#$infofile=$ARGV[0];
#$what_ask=$ARGV[1];
#  $file_out="/XS";
#  $file_out="dither.ps/CPS";
#  $file_out2="test2.ps/CPS";
  pgbeg(0,$file_out,1,1);
pgsvp(0.15,0.85,0,1);
#pgswin(-$xsize/2,$xsize/2,-$ysize/2,$ysize/2);
#pgswin(-$xsize/2,$xsize/2,-$ysize/2,$ysize/2);
pgswin($xmax,$xmin,$ymin,$ymax);
 for ($i=0;$i<$nd;$i++) {
   pgsci(1+$i);
#   pattern_dith($xd[$i],$yd[$i]);
   pattern_dith_file($xd[$i],$yd[$i]);
 }
  pgsci(1);
  pgclos();
  pgend();

exit;

sub pattern {

    $x0=-15/60;
    $y0=-27/60;
    $x00=-18/60;
    $y00=-27/60;
    $n=11;
    $n_tot=0;
    $dist=3*1.2;
    pgsfs(2);
    for ($i=0;$i<21;$i++) {
	$j=0;	
	while($j<$n) {
	    $x=$x0+$dist*$j/60;
	    pgcirc($x,$y0,$dist/($ratio*120));	    
	    $n_tot++;
#	    print "$x $y0\n";
	    $j++;
	}
	if ($i<10) {
	    $n=$n+1;
	    $x0=$x0-$dist/120;
	} else {
	    $n=$n-1;
	    $x0=$x0+$dist/120;
	}
	$y0=$y0+($dist/60)*0.87;
    }
    print "N.PPAK=$n_tot\n";
}

sub pattern_dith {
    my $xc=@_[0];
    my $yc=@_[1];
    my $x0=-15/60+$xc/60;
    my $y0=-27/60+$yc/60;
    my $x00=-18/60;
    my $y00=-27/60;
    my $n=11;
    my $n_tot=0;
    my $dist=3*1.2;
#    my $ratio=1.1111;
    my $i,$j;
    pgsfs(1);
    for ($i=0;$i<21;$i++) {
	$j=0;	
	while($j<$n) {
	    $x=$x0+$dist*$j/60;
	    pgcirc($x,$y0,$dist/($ratio*120));	    
	    $n_tot++;
#	    print "$x $y0\n";
	    $j++;
	}
	if ($i<10) {
	    $n=$n+1;
	    $x0=$x0-$dist/120;
	} else {
	    $n=$n-1;
	    $x0=$x0+$dist/120;
	}
	$y0=$y0+($dist/60)*0.87;
    }
    print "N.PPAK=$n_tot\n";
}

sub pattern_dith_file {
    my $xc=@_[0];
    my $yc=@_[1];
    my $x0=-15/60+$xc/60;
    my $y0=-27/60+$yc/60;
    my $x00=-18/60;
    my $y00=-27/60;
    my $n=11;
    my $n_tot=0;
    my $dist=3*1.2;
#    my $ratio=1.1111;
    my $i,$j;
    open(FH,"</home/sanchez/sda2/code/lyon/v3d/data/ppak_pt_arc.txt");
    $line=<FH>;
    chop($line);
    @data=split(" ",$line);
    $size=$data[1];
    $size=$size*$ratio;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$x=$xc+$data[1];
	$y=$yc+$data[2];
#	print "$x $y $size\n";
	pgcirc($x,$y,$size);	    
    }
    close(FH);
    print "N.PPAK=$n_tot\n";
}


sub pattern_pmas_8 {
    $x0=5/60;
    $y0=2/60;
    $x00=-18/60;
    $y00=-27/60;
    $n=11;
    pgsfs(2);
    for ($i=0;$i<16;$i++) {
	for ($j=0;$j<16;$j++) {
	    $x1=$x0+0.5*$i/60;
	    $y1=$y0+0.5*$j/60;
	    $x2=$x0+0.5*$i/60+0.5/60;
	    $y2=$y0+0.5*$j/60+0.5/60;
	    pgrect($x1,$x2,$y1,$y2);
	}
    }
}


sub pattern_pmas_16 {
    $x0=-12/60;
    $y0=-15/60;
    $x00=-18/60;
    $y00=-27/60;
    $n=11;
    pgsfs(2);
    for ($i=0;$i<16;$i++) {
	for ($j=0;$j<16;$j++) {
	    $x1=$x0+1.0*$i/60;
	    $y1=$y0+1.0*$j/60;
	    $x2=$x0+1.0*$i/60+1/60;
	    $y2=$y0+1.0*$j/60+1/60;
	    pgrect($x1,$x2,$y1,$y2);
	}
    }
}



sub pattern_pmas_32 {
    $x0=-20/60;
    $y0=12/60;
    $x00=-18/60;
    $y00=-27/60;
    $n=11;
    pgsfs(2);
    for ($i=0;$i<32;$i++) {
	for ($j=0;$j<32;$j++) {
	    $x1=$x0+1.0*$i/60;
	    $y1=$y0+1.0*$j/60;
	    $x2=$x0+1.0*$i/60+1/60;
	    $y2=$y0+1.0*$j/60+1/60;
	    pgrect($x1,$x2,$y1,$y2);
	}
    }
}

