#!/usr/bin/perl -w

############################################
#                                          #
# Script to compile the Itsas PFC template #
#                                          #
############################################

use strict;

#
# Check that main file exists
#
my $f = 'note'; # Name of main file (you can change it, of course)
die "File $f.tex does not exist!\n" unless (-f "$f.tex");

#
# Get the page size from the .tex
#
my $psize = 'a4';                                          # A4 by default
$psize    = 'b5' if (`grep documentclass $f.tex` =~ /b5/); # B5, if specified

#
# Define PDF viewer
#
my $seepdf = 'evince';

#
# Other variables
#
my $tmp    = 'tmp';
my $outpdf = $f.'.pdf';

#
# Create the to-do list
#
my @list;

push(@list,"cp -f $f.tex $tmp.tex");                                         # temporal copy of file to compile
push(@list,"latex $tmp");                                                    # run latex
#push(@list,"bibtex $tmp");                                                   # run bibtex
#push(@list,"latex $tmp");                                                    # second compilation, for cross-references
push(@list,"latex $tmp");                                                    # run latex again
push(@list,"dvips -t $psize $tmp.dvi -o $tmp.ps > /dev/null");               # convert DVI to PS
push(@list,"ps2pdf $tmp.ps");                                                # convert PS to PDF
push(@list,"mv -f $tmp.pdf $outpdf");                                        # save a copy PDF in $outpdf
push(@list,"$seepdf $outpdf &");                                             # open PDF and go on
push(@list,"rm -f $tmp.*");                                                  # delete rubish

#
# To time the steps
#
my $hasi = time(); # +%s date of beginning time
my $oldt = $hasi;
my $str  = '';     # string with all steps

#
# Run commands in list
#
foreach (@list)
{
  my $do = $_;

  my $ok  = 'OK';          # whether step run smoothly
  my $sys = system "$do";  # execute command, and catch exit status
  my $t   = time();
  my $dt  = $t-$oldt;      # seconds this step takes
  my $dt0 = $t-$hasi;      # seconds taken up to this point
  $oldt   = $t;            # ending time for this step
  if ($sys) { $ok = 'KO' } # not-OK if exit status was not 0

  #
  # Save info for later printing
  #
  $str .= sprintf "%-50s   [ %2s ] %6is %6is\n",$do,$ok,$dt,$dt0;

  last if $sys;            # exit if last step failed
};

#
# Print abstract
#
printf "\nABSTRACT [ %1s ]:\n\n%-50s %8s %7s %7s\n\n%1s\n",$f,'Command','Status','Step','Total',$str;
