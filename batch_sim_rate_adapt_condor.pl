#!/usr/bin/perl

use strict;

#####
## mother script
my $mother_script = "batch_sim_rate_adapt";

#####
## files
my @files = ('rx_syms_plane1.dat', 'rx_syms_plane2.dat', 'rx_syms_plane3.dat', 'rx_syms_plane4.dat');


my $file_i = 0;
foreach my $file (@files) {
    $file_i ++;
    print $file_i.". ".$file."\n";


    ## sh
    system("sed 's/$mother_script/$mother_script$file_i/g;' $mother_script.mother.sh > $mother_script$file_i.sh");


    ## condor
    system("sed 's/$mother_script/$mother_script$file_i/g;s/condor/condor$file_i/g;' $mother_script.mother.condor > $mother_script$file_i.condor");

    ## matlab
    system("sed 's/$mother_script/$mother_script$file_i/g;s/XXXXX/$file/g;s/batch_sim/batch_sim$file_i/g;' $mother_script.mother.m > $mother_script$file_i.m");


    ## submit jobs
    system("condor_submit $mother_script$file_i.condor");
}

print `condor_q`;

