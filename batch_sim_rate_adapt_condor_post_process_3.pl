#!/usr/bin/perl

use strict;


#####
## files
my $dir = "./OUTPUT_sim_3";

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
    
    open FH, "<$dir/condor$file_i.output";
    while(<FH>) {
        if($_ =~ /^\s+(\d+\.*\d*e*\+*\d*)/) {
            # print $_;
            # print "".($1+0)."\n";
            my $tput = $1+0;
            print "$file_ind, XXX, Oracle, XXX, ".$tput."\n";
            last;
        }
    }
    close FH;
    # print "condor$file_i.output\n";
    # return;

    foreach my $method_SNR (@method_SNRs) {
        foreach my $method_RA (@method_RAs) {
            foreach my $method_PRED (@method_PREDs) {
                print "$file_ind, $method_SNR, $method_RA, $method_PRED, ";
                foreach my $threshold (@{$thresholds{$method_RA}}) {
                
                    $file_i ++;

                    open FH, "<$dir/condor$file_i.output" or die $!;
                    while(<FH>) {
                        if($_ =~ /^\s+(\d+\.*\d*e*\+*\d*)/) {
                            my $tput = $1+0;
                            # print "$file_ind, $method_SNR, $method_RA, $method_PRED, ".$1."\n";
                            print $tput.", ";
                        }
                    }
                    close FH;

                }
                print "\n";
            }
        }
    }  
}

