#!/usr/bin/perl
use strict;
use PDF::Reuse;

prFile("out.pdf");

for(@ARGV) {
prDoc($_);
}

prEnd(); 
