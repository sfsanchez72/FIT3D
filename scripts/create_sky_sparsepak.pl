#!/usr/bin/perl

use v5.10;
use strict;
use warnings;
use PDL;
use PDL::Image2D;
use File::Basename qw(basename);

# We need PDL >=2.4.10 due to bug fixes contained therein
die "Please upgrade to at least version 2.4.10 of PDL." unless "$PDL::VERSION" ge '2.4.10';

if (@ARGV < 4)
{
	print 'USE: ', basename($0), " INPUT.fits OBJECT.fits SKY.fits PLOT\n";
	exit -1;
}

my (
	$inputfile,    # name of RSS input
	$objfile,      # name of file to store object (science) fiber data
	$skyfile,      # name of file to store sky fiber data
	$spy           # flag to graph results
) = @ARGV;

# It is **very** important to keep these indices in order!
my $skyFiberIndices = pdl(1,15,21,36,53,69,79);

my $input = rfits($inputfile);

# Grab the sky fiber data
my $sky = $input->dice_axis(1,$skyFiberIndices);

# Grab the science data
# The indices of the science fibers are at the complementary locations of the sky fibers
my $tmp = sequence(82);
my $object = $input->dice_axis(1, $tmp->where($tmp->in($skyFiberIndices)->inplace->not));

# Remove any (nearly) zero values from the sky data by averaging with their nearest neighbors (essentially a mean dithering)
$sky .= $sky->patch2d(approx($sky,0,1e-7));

# compute the median and standard deviation of the sky flux at each wavelength
my (undef, undef, $med, undef, undef, undef, $stdDev) = $sky->xchg(0,1)->statsover;

my $sky_old;
if ($spy == 1)
{
	$sky_old = $sky->copy();
}

# Perform 2-sigma rejection, replacing with the wavelength-specific median value
my $mask = abs($sky - $med) > 2 * $stdDev;
$sky->where($mask) .= $med->dummy(1)->mult($mask,0)->where($mask);

# Interpolate the sky data over all 82 fibers to create a sky image
# the same size as the data image which can be easily subtracted from the object data
my ($new_sky, undef) = interpolate(sequence(82)->dummy(0), $skyFiberIndices, $sky->xchg(0,1));

if ($spy == 1)
{
	use PDL::Graphics::PGPLOT;
	
	my ($column, $command, $maxColumn) = (0, '');
	$maxColumn = $sky_old->getdim(0);
	
	COLUMN_LOOP: while(1)
	{
		last if $column == $maxColumn;
		
		env(0, 82,
			$sky_old->slice("($column),:")->min,
			$sky_old->slice("($column),:")->max,
			{'linewidth' => 2, 'charsize' => 1.6, 'xtitle' => 'SpecID', 'ytitle' => 'Counts', 'title' => "Column $column"}
		);
		points($skyFiberIndices, $sky_old->slice("($column),:"), {'symbol' => 'plus', 'charsize' => 2.2, 'color' => 'blue'});
		points($skyFiberIndices, $sky->slice("($column),:"), {'symbol' => 'dot', 'charsize' => 2.2, 'color' => 'green'});
		line($tmp, $new_sky->slice("($column),:"), {'charsize' => 2.2, 'color' => 'red'});
		
		while(1)
		{
			print 'Press enter to advance to the next column, q to quit, or enter a SpecID [min=0, max=', $maxColumn-1, ']: ';
			chomp($command = <STDIN>);
			
			# Show the next column if just enter was pressed 
			if ($command eq '')
			{
				$column++;
				next COLUMN_LOOP;
			}
			
			# stop everything if 'q[uit]' was entered
			last COLUMN_LOOP if $command =~ m/^\s*q(?:uit)?\s*$/i;
			
			# Show the specified column, if a valid digit sequence was entered
			if ($command =~ m/^\s*(\d+)\s*$/)
			{
				# keep this inside to prevent Boolean short-circuiting
				if ($1 >= 0 && $1 <= $maxColumn)
				{
					$column = $1;
					next COLUMN_LOOP;
				}
			}
		}
	}
}

# Save the object fiber data [n=75]
$object->wfits($objfile);

# Save the sky data that was interpolated into the object fibers [n=75]
$new_sky->dice_axis(1, $tmp->where($tmp->in($skyFiberIndices)->not))->wfits($skyfile);
