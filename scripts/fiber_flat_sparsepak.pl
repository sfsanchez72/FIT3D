#!/usr/bin/perl

use v5.10;
use strict;
use warnings;
use PDL;
use PDL::IO::FITS;
use File::Basename qw(basename);

# We need PDL >=2.4.10 due to bug fixes contained therein
die "Please upgrade to at least version 2.4.10 of PDL." unless "$PDL::VERSION" ge '2.4.10';

if (@ARGV < 1)
{
	print 'USE: ', basename($0), " SKYFLAT.fits FIBERFLAT.fits [LIMIT MASK]\n";
	exit -1;
}

my (
	$skyflatFile,			# The skyflat file
	$fiberFlatFile,			# The flattened version of $infile
	$limit,					# The lower cutoff as a percent of the median sky value (default=0.2)
	$maskFile				# The name of the file where the masked pixels are stored (default='SKYFLAT.mask.fiberflat.fits')
) = @ARGV;

$limit //= 0.2;

# Create a unique name for the mask file unless the user specified one
unless (defined $maskFile)
{
	$maskFile = (split(quotemeta '.', $fiberFlatFile))[0] . '.mask.fiberflat.fits';
}

my $skyFlux = rfits($skyflatFile);

# compute the median flux at each wavelength
my $medSky = $skyFlux->xchg(0,1)->medover();

# divide each spectrum by its median flux
my $fiberFlat = $skyFlux / $medSky;

# set the mask values to 1
my $mask = zeroes(long, $skyFlux->dims());
$mask->where($fiberFlat > $limit) .= 1;

$fiberFlat->sethdr($skyFlux->hdr);
$fiberFlat->wfits($fiberFlatFile);

$mask->sethdr($skyFlux->hdr);
$mask->wfits($maskFile);
