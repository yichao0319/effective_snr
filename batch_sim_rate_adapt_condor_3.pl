#!/usr/bin/perl

use strict;

#####
## mother script
my $mother_script = "batch_sim_rate_adapt_3";

#####
## files
my @files = 1 .. 4;
my @method_SNRs = ("preambleSNR", "allSNR", "allSNRoracle");
# my @method_RAs = ("Thresholding");
my @method_RAs = ("Thresholding", "Probability_1", "Probability_2");
my @method_PREDs = ("currPkt", "prevPkt", "EWMA", "HW");
my @thresholds1 = (map 0.01*$_, 80..100);
my @thresholds2 = (map 0.01*$_, 95..100);
my @thresholds3 = (map 0.01*$_, 95..100);
my %thresholds = ($method_RAs[0] => \@thresholds1, $method_RAs[1] => \@thresholds2, $method_RAs[2] => \@thresholds3);


my $file_i = 0;
foreach my $file_ind (@files) {

    ## Oracle
    $file_i ++;
    print $file_i.". ".$file_ind."\n";


    ## sh
    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;' $mother_script.mother.sh > tmp_$mother_script$file_i.sh");


    ## condor
    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;s/condor/condor$file_i/g;' $mother_script.mother.condor > tmp_$mother_script$file_i.condor");


    ## matlab
    my $method_SNR = $method_SNRs[0];
    my $method_RA = "Oracle";
    my $method_PRED = $method_PREDs[0];
    my $threshold = 0;
    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;s/file_ind/$file_ind/g;s/method_SNR/$method_SNR/g;s/method_RA/$method_RA/g;s/method_PRED/$method_PRED/g;s/threshold/$threshold/g;' $mother_script.mother.m > tmp_$mother_script$file_i.m");


    ## submit jobs
    system("condor_submit tmp_$mother_script$file_i.condor");


    foreach my $method_SNR (@method_SNRs) {
        foreach my $method_RA (@method_RAs) {
            foreach my $method_PRED (@method_PREDs) {
                foreach my $threshold (@{$thresholds{$method_RA}}) {
                    $file_i ++;
                    print $file_i.". ".$file_ind."\n";


                    ## sh
                    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;' $mother_script.mother.sh > tmp_$mother_script$file_i.sh");


                    ## condor
                    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;s/condor/condor$file_i/g;' $mother_script.mother.condor > tmp_$mother_script$file_i.condor");


                    ## matlab
                    system("sed 's/$mother_script/tmp_$mother_script$file_i/g;s/file_ind/$file_ind/g;s/method_SNR/$method_SNR/g;s/method_RA/$method_RA/g;s/method_PRED/$method_PRED/g;s/threshold/$threshold/g;' $mother_script.mother.m > tmp_$mother_script$file_i.m");


                    ## submit jobs
                    system("condor_submit tmp_$mother_script$file_i.condor");
                }
            }
        }
    }  
}

print `condor_q`;

