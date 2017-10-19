#!/usr/bin/perl

use strict;
use warnings;
use File::Basename qw(basename);
use PDL;
use PDL::IO::FITS;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/usr/local/lib/R3D/my.pl";

# We need PDL >=2.4.10 due to bug fixes contained therein
die "Please upgrade to at least version 2.4.10 of PDL." unless "$PDL::VERSION" ge '2.4.10';

if (@ARGV < 4)
{
	print 'USE: ', basename($0), " INPUT.fits OBJECT.fits SKY.fits PLOT\n";
	exit -1;
}

my (
	$inputfile,    # name of the RSS file to split
	$objfile,      # name of the file to store the science data
	$skyfile,      # name of the file to store the sky data
	$spy           # flag to plot
) = @ARGV;

# It is **very** important to keep these indices in order!
my $skyFiberIndices = pdl(1,15,21,36,53,69,79);

my $input = rfits($inputfile);

# Grab the sky fiber data
my $sky = $input->dice_axis(1,$skyFiberIndices);

# Grab the science data
# The indices of the science fibers are at the complementary locations of the sky fibers
my $tmp = sequence(82);
my $data = $input->dice_axis(1, $tmp->where($tmp->in($skyFiberIndices)->not));

# Plot the sky data, if asked
if ($spy)
{
	use PDL::Graphics::PGPLOT;
	
	for my $i (0..$sky->getdim(1)-1)
	{
		line($sky->slice(":,($i)"), {'charsize' => 1.6, 'linewidth' => 2, 'title' => "Sky $i", 'xtitle' => 'ID', 'ytitle' => 'Counts'});
		<STDIN>;
	}
}

$sky->wfits($skyfile);
$data->wfits($objfile);
