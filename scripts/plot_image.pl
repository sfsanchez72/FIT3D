#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
    eval 'exec perl -S $0 "$@"'
        if 0;

use strict;
$|++;

use PDL::Doc::Perldl;
use File::Basename;

use vars qw( $VERSION );
$VERSION = '0.2';

use PDL::Config;
my $bvalflag = $PDL::Config{WITH_BADVAL} || 0;

my %options = 
    ( a => \&apropos, 
      h => \&help, s => \&sig, u => \&usage );

my $name = basename( $0 );
my $usage = <<"EOH";
Usage: $name [-a] [-h] [-s] [-u] <string>

This program provides command-line access to the PDL documentation.
If no flag is specified, -h is assumed.

  -a (apropos) searches the documentation for the string
  -h (help)    prints the help for the function/module/document
  -s (sig)     prints the signature of the function
  -u (usage)   gives usage information on the function

EOH

my $oflag = $#ARGV > -1 ? substr($ARGV[0],0,1) eq "-" : 0;
die $usage unless ($#ARGV == 0 and not $oflag) or ($#ARGV == 1 and $oflag);

my $option = "h";
if ( $oflag ) {
    $option = substr($ARGV[0],1,1);
    die $usage unless exists $options{$option};
    shift @ARGV;
}

&{$options{$option}}( $ARGV[0] );

exit;

__END__

=head1 NAME

pdldoc - shell interface to PDL documentation

=head1 SYNOPSIS

B<pdldoc> <text>

B<pdldoc> [B<-a>] [B<-h>] [B<-s>] [B<-u>] <text>

=head1 DESCRIPTION

The aim of B<pdldoc> is to provide the same functionality
as the C<apropos>, C<help>, C<sig>, 
and C<usage> commands available in the L<perldl|PDL::perldl> shell.

Think of it as the PDL equivalent of C<perldoc -f>.

=head1 OPTIONS

=over 5

=item B<-h> help

print documentation about a PDL function or module or show a PDL manual.
This is the default option.

=item B<-a> apropos

Regex search PDL documentation database.

=item B<-s> sig

prints signature of PDL function.

=item B<-u> usage

Prints usage information for a PDL function.

=back

=head1 VERSION

This is pdldoc v0.2.

=head1 AUTHOR

Doug Burke <burke@ifa.hawaii.edu>.

=cut

