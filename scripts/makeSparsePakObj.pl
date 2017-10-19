#!/usr/bin/perl
use v5.10;
use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Long;
use Carp qw(croak);
use Cwd qw(cwd);
use File::Copy qw(copy);
use PDL;
use PDL::IO::FITS;

# We need PDL >=2.4.10 due to bug fixes contained therein
die "Please upgrade to at least version 2.4.10 of PDL." unless "$PDL::VERSION" ge '2.4.10';

my ($darkFileName, $biasFileName, $makeMasterSubtrahend, $objectName, $contSuffix, $arcSuffix);

&parseCommandLineArguments();
&makeMasterSubtrahend(), exit if $makeMasterSubtrahend;

for my $type ('', $arcSuffix, $contSuffix)
{
	croak "$type doesn't exist or is empty.\n" unless -s "${objectName}${type}.list";
	&execute('imcombine.pl', "${objectName}${type}.list ${objectName}${type}_raw.fits");
	&fixupObj("${objectName}${type}_raw.fits", "${objectName}${type}.fits");
}

sub execute
{
	my $prog   = shift;
	my $params = shift;
	system("$prog $params");
	croak "\n\nError executing \n\t'$prog $params'\n\n" if (($? >> 8) != 0 || $? == -1 || ($? & 127) != 0);
}

sub fixupObj
{
	my $inName  = shift;
	my $outName = shift;
	
	$outName //= $inName;

	my $pdlIn = rfits($inName, {'hdrcpy' => 1});
	$pdlIn *= $pdlIn->hdr->{'GAIN'};				# convert from adu to electrons

	# Flip and rotate so that the blue end is on the left and fiber 1 at the bottom
	$pdlIn
	->slice('-1:0:-1')								# flip vertically
	->xchg(0,1)->slice('-1:0:-1')					# rotate and flip vertically
	->wfits($outName);
}

sub makeMasterSubtrahend()
{
	open my $fdOut, '>', "$objectName.list" or die "Unable to open $objectName.list: $!\n";
	
	my $rootName;
	
	for my $file ($darkFileName, $biasFileName)
	{
		if ($file =~ m/(.+)\.list$/)
		{
			$rootName = $1;
			croak "$file doesn't exist or is empty.\n" unless -s $file;
			&execute('imcombine.pl', "$file $rootName.fits");
			$file = "$rootName.fits";
		}
		
		# if the user didn't supply a list file, then it must be a fits file
		die "$file must be a fits file.\n" unless $file =~ m/\.fits$/;
		print $fdOut "$file\n";
	}
	
	close $fdOut;

	&execute('imcombine.pl', "$objectName.list ${objectName}_raw.fits");
	&fixupObj("${objectName}_raw.fits", "$objectName.fits");
}

sub parseCommandLineArguments()
{
	my $help;
	
	my $ret = GetOptions(
		'dark=s'		=> \$darkFileName,
		'bias|zero=s'	=> \$biasFileName,
		'master'		=> \$makeMasterSubtrahend,
		'help|?|h'		=> \$help,
		'contsuffix'	=> \$contSuffix,
		'arcsuffix'		=> \$arcSuffix
	);
	
	&showHelp(), exit if $help;
	
	unless (@ARGV == 1 && $ret == 1)
	{
		&showHelp();
		exit -1;
	}

	$objectName = $ARGV[0];
	
	# remove the fits extension, if present
	$objectName =~ s/\.fits$//;
	
	$darkFileName //= 'dark.fits';
	$biasFileName //= 'bias.fits';
	$contSuffix //= '_cont';
	$arcSuffix //= '_arc';
}

sub showHelp()
{
	my $basename = basename($0);
	print<<__END__;

Usage: $basename \[options\] objectName
	
	--dark              Name of the dark file                   (default='dark.fits')
	--bias              Name of the bias file                   (default='bias.fits')
	--zero              Synonym for --bias
	--master            Compute the master subtrahend file      (default=false)
	--contsuffix        The suffix for the continuum file       (default='_cont')
	--arcsuffix         The suffix for the arc file             (default='_arc')
	--help              Show this help

	If the name given in --dark or --bias ends in '.list', then it is assumed to be a list file
	containing the names of the files to be median combined into the respective image. As a courtesy,
	imcombine will be called against it to create the corresponding fits image.
	
	The --master option is used for creating a single file which consists of a median-combination
	of the bias and dark files. This master subtrahend file can then be used with imarith to simultaneously
	remove the dark current and bias from an image.

Examples
	Create a master subtrahend file called 'master_sub.fits'
	
		$basename --master --dark=dark.list --bias=bias.list master_sub
		
	Create the science, continuum, and arc images for a target object 'NGC224'.
	
		$basename NGC224
		
		This creates NGC224.fits NGC224_cont.fits NGC224_arc.fits assuming there exist the list files
		NGC224.list, NGC224_cont.list, and NGC224_arc.list.
	
	Create the science, continuum, and arc images for a target object 'NGC224' using a non-default prefix
	
		$basename --contsuffix='.cont' --arcsuffix='.arc' NGC224
		
		This creates NGC224.fits NGC224.cont.fits NGC224.arc.fits assuming there exist the list files
		NGC224.list, NGC224.cont.list, and NGC224.arc.list.
__END__
}

