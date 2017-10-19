#!/star/Perl/bin/perl





use Astro::Telescope;
use Astro::Coords;


if ($#ARGV<4) {
    print "USE: get_heliocentric_vel.pl OBSERVATORY[CAHA]  RA  DEC  EQUINOX[J2000] OBS_DATE_UTC[2010-10-10T22:00:00]\n";
    exit;
}
$OBS=$ARGV[0];
$RA=$ARGV[1];
$DEC=$ARGV[2];
$EQUI=$ARGV[3];
$DATE=$ARGV[4];
($date,$time)=split("T",$DATE);
($year,$month,$day)=split(/\-/,$date);
($hour,$minute,$second)=split(/\:/,$time);
#print "($year,$month,$day) ($hour,$minute,$second)\n";
#exit;


$long=0.64966917648148148146;
$lat=0.04444768074074074071;
$alt=2168;
$tel = new Astro::Telescope(Name => 'CAHA', Long => $long, Lat => $lat, Alt => $alt );
$latitude = $tel->lat;
$longitude = $tel->long;
$altitude = $tel->alt;
 %limits = $tel->limits;
#print "$latitude $longitude $altitude $limits\n";

#$c = new Astro::Coords( name => "TARGET", ra   => '05:22:56', dec  => '-26:20:40.4', type => 'J2000',units=> 'sexagesimal');
#$c->telescope( new Astro::Telescope( 'CAHA' ));

$c = new Astro::Coords( name => "TARGET", ra   => $RA, dec  => $DEC, type => $EQUI ,units=> 'sexagesimal');
$c->telescope( new Astro::Telescope( $OBS ));

my $dt = DateTime->new(
      year       => $year,
      month      => $month,
      day        => $day,
      hour       => $hour,
      minute     => $minute,
      second     => $second,
      nanosecond => 00,
      time_zone  => 'UTC',
  );
$epoch=2000;

$date=$dt;
#$date = DateTime->from_epoch( epoch => $epoch, time_zone => 'UTC' );
#$date="2010-06-10T22:00:10";
#$date="1970-01-01T00:33:20";
#$date="2010-10-10T22:00:00";

$c->datetime( $date );


$vd = $c->vdiff( 'HEL', 'TOP' );

print "$vd\n";

exit;
