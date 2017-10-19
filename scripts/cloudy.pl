#!/usr/bin/perl -w

# This is an example of how to write a perl program to drive cloudy.
# This method is an alternative to writing a new fortran main program
# to run a grid of models (as described in the manual).  With a little
# practice in perl, it's easy to write a program to generate a grid of
# models, and then another program to extract the desired information
# from the output files that are generated. 

# To run this program, save it to a file called runcloudy, then make
# it executable with "chmod +x runcloudy".  Change lines 1 and 72 to
# give the correct path names to perl and cloudy on your system, then
# type runcloudy to start the grid running.

# The first line announces that this is a perl script rather than a
# shell script.   Change the path name to point to where perl lives on
# your network.  The -w option tells perl to look for problems or
# errors in the script and print warning messages.

# Open a log file to save comments about the results of each cloudy
# run.  The ">>" means to open the file for appending, so if the file
# already exists it won't be clobbered.  Then, to append text to the
# logfile, we use the command:
# print LOGFILE "text to append to logfile."
# Note that normally the output of a print command goes to standard
# output, unless a filehandle such as LOGFILE is specified.

open (LOGFILE, ">>cloudy.log");

# Main program loop.  In this example, we will calculate a 2x2 grid by
# varying the ionization parameter and the density.  This is easily
# done by nesting two "for" loops.  In perl, variable names begin with
# the $ sign.
#
# This loop structure can easily be extended to grids of any
# dimensionality, and any parameters can be looped over, such as
# blackbody temperature, continuum luminosity, abundances, alpha_ox
# (with the "AGN" command), etc.  The loop increment command $i++ can 
# also be written as $i=$i+1, and this can be generalized to any
# desired increment, for example $i=$i+0.2 for a finer grid.

for ($i = 3; $i <= 4; $i++) {
    $ion = -1 * $i;
    
    for ($n = 3; $n <= 4; $n++) {
        
# For each set of values of $i and $n, create filenames for the output
# files. In perl, variable names within double-quoted strings are
# replaced by the value of the variable.  So, for example, if $i=3 and
# $n=4, the string "nlr-$i-$n.out" gets evaluated as "nlr-3-4.out".
# This is an easy way to index the output filenames by the grid
# variables.  The $params string just lists the values of the grid
# parameters so that they can be written to the logfile and the
# screen.  

        $params = "log U = $ion, log n(H) = $n";
        $outfile = "nlr-$i-$n.out";
        $ovrfile = "nlr-$i-$n.ovr";
        $rltfile = "nlr-$i-$n.rlt";
        
        print "Cloudy run begun with $params.\n";
        
# Now, start the cloudy process running, and give it the filehandle 
# CLOUDY so that we can send commands to it from perl.  The pipe (|)
# symbol before the pathname to cloudy tells perl that any 
# "print CLOUDY" arguments will be piped into the cloudy process.
# Change the pathname to the cloudy executable to point to wherever
# cloudy lives in your directory structure.  The "> $outfile" tells
# perl that the standard output of cloudy will be redirected to the
# file whose name is given by the variable $outfile, defined above.

        open (CLOUDY, "|/home/sanchez/sda1/cloudy/c07.02.01/source/cloudy.exe > $outfile") ||
            die "Can't open output file: $!";
        
# Enter the commands to cloudy.  Each command is given as a
# double-quoted string argument to the "print CLOUDY" command.
# Each string needs to be terminated with the newline "\n" character,
# just like in C/C++.  Within a double-quoted string, we can include 
# variable names, which get evaluated as their current values.  So if
# $n=4, then the string "hden $n" gets sent to cloudy as the command 
# "hden 4".  A double-quote character can be included within a string
# by preceding with a backslash: \".  The final command to cloudy is 
# a newline character by itself, which tells cloudy that the last 
# command has been entered, and the calculation begins.  When it's
# finished, close the CLOUDY process.  Don't forget to end each
# command line with a semicolon; if you leave it out, perl will probably 
# generate a syntax error. 

        print CLOUDY "title nlr model with $params \n";
        print CLOUDY "table agn \n";
        print CLOUDY "constant density\n";
        print CLOUDY "hden $n \n";
        print CLOUDY "iterate to convergence\n";
        print CLOUDY "ionization parameter $ion \n";
        print CLOUDY "punch overview last file = \"$ovrfile\"\n";
        print CLOUDY "punch results last file = \"$rltfile\"\n";
        print CLOUDY "print last\n";
        print CLOUDY "\n";
        
        close (CLOUDY);
        
# Now, print warning messages to the logfile and to the screen.
# Unix shell commands can be accessed within perl by enclosing in
# backquotes, so to get the current date and time as a string one
# can use `date`.   The final line of the $outfile file contains the
# numbers of warnings and cautions: we will copy that line to the
# variable $warnings and print it to the logfile and to the screen.
# When the entire grid is done, close the logfile and end.

        $date = `date`;
        chomp $date;
        $warnings = `tail -1 $outfile`;
        print LOGFILE "Cloudy run with $params: \n";
        print LOGFILE "Output file $outfile created on $date. \n";
        print LOGFILE $warnings;
        print LOGFILE "\n";
        print $warnings; 
        print "Cloudy run completed with $params. \n";
        print "\n";
    }
}

close (LOGFILE);

# If you're trying this for the first time, take a look at the .out
# files that have been generated, to see how the "print CLOUDY"
# commands above have been translated into actual cloudy input lines.
# Also, look at the cloudy.log file to see if there were any warnings.

# by A. Barth, abarth@cfa.harvard.edu
