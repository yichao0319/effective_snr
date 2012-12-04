#!/usr/bin/perl

use strict;

#####
## mother script
my $mother_script = "batch_analyze_csi";

#####
## files
# my @files = ('RXDATA.static.d0.pdat', 'RXDATA.static.d6.pdat', 'RXDATA.mobile.2.pdat', 'RXDATA.mobile.b1.pdat', 'RXDATA.mobile.b2.pdat', 'RXDATA.s96.pdat', 'RXDATA.s96-2.pdat', 'RXDATA.static.highSNR.pdat', 'RXDATA.static.highSNR.2.pdat', 'RXDATA.static.highSNR.d1.pdat', 'RXDATA.static.highSNR.d1.2.pdat', 'RXDATA.static.highSNR.d2.pdat', 'RXDATA.static.highSNR.d6.pdat', 'RXDATA.static.highSNR.d6.2.pdat', 'RXDATA.mobile.highSNR.pdat', 'RXDATA.mobile.highSNR.1.pdat', 'RXDATA.mobile.highSNR.b1.pdat', 'RXDATA.mobile.highSNR.b2.pdat');
# my @files = ('RXDATA.5MHz.pdat', 'RXDATA.10MHz.pdat');
my @files = ('RXDATA.5MHz.2.pdat');


my $file_i = 0;
foreach my $file (@files) {
    $file_i ++;
    print $file_i.". ".$file."\n";


    ## sh
    system("sed 's/batch_analyze_csi/batch_analyze_csi$file_i/g;' $mother_script.mother.sh > $mother_script$file_i.sh");


    ## condor
    system("sed 's/batch_analyze_csi/batch_analyze_csi$file_i/g;s/condor/condor$file_i/g;' $mother_script.mother.condor > $mother_script$file_i.condor");

    ## matlab
    system("sed 's/XXXXX/$file/g;s/batch_out/batch_out$file_i/g;s/batch_err_ratio/batch_err_ratio$file_i/g;' $mother_script.mother.m > $mother_script$file_i.m");


    ## submit jobs
    system("condor_submit $mother_script$file_i.condor");
}

print `condor_q`;

