#!/bin/perl 

use strict;

my $DEBUG0 = 0;
my $DEBUG1 = 1;
my $DEBUG2 = 1;


#####
## input parameters
if($#ARGV != -1) {
    die "wrong number of input\n";
}


#####
## global variables
my $raw_dir = "./rawofdm.mobile.traces";
my $output_dir = "./PARSEDDATA";

my @rx_symbol_files = (
    # "RXDATA.static.d0", "RXDATA.static.d6", 
    # "RXDATA.mobile.1", "RXDATA.mobile.2", 
    # "RXDATA.mobile.b1", "RXDATA.mobile.b2"
    # "rxdata_tmp",
    # "tmp_rxdata_s16",
    # "tmp_rxdata_s32"
    # "tmp_rxdata_s128"
    # "tmp_rxdata_s200",
    # "tmp_rxdata_s400"
    # "tmp_rxdata_s240"
    # "RXDATA.static.s16",
    # "RXDATA.static.s32",
    # "RXDATA.static.s64",
    # "RXDATA.s96",
    # "RXDATA.s96-2"
    # "RXDATA.static.highSNR"
    # "RXDATA.static.highSNR.2",
    # "RXDATA.static.highSNR.d1",
    # "RXDATA.static.highSNR.d1.2",
    # "RXDATA.static.highSNR.d2",
    # "RXDATA.static.highSNR.d6",
    # "RXDATA.static.highSNR.d6.2",
    # "RXDATA.mobile.highSNR",
    # "RXDATA.mobile.highSNR.1",
    # "RXDATA.mobile.highSNR.b1",
    # "RXDATA.mobile.highSNR.b2"
    # "RXDATA.5MHz",
    # "RXDATA.10MHz"
    "RXDATA.5MHz.2"
    );
my @snr_files = (
    # "SNRDATA.static.d0", "SNRDATA.static.d6", 
    # "SNRDATA.mobile.1", "SNRDATA.mobile.2", 
    # "SNRDATA.mobile.b1", "SNRDATA.mobile.b2"
    # "snrdata_tmp",
    # "tmp_snrdata_s16",
    # "tmp_snrdata_s32"
    # "tmp_snrdata_s200",
    # "tmp_snrdata_s400"
    # "tmp_snrdata_s240"
    # "SNRDATA.static.s16",
    # "SNRDATA.static.s32",
    # "SNRDATA.static.s64"
    # "SNRDATA.s96",
    # "SNRDATA.s96-2"
    );

#####
## main

## rx symbols
foreach my $rx_symbol_file (@rx_symbol_files) {
    my $cmd = "parse.symbols $raw_dir/$rx_symbol_file > $output_dir/$rx_symbol_file.pdat";
    `$cmd`;
}

## rx symbols
foreach my $snr_file (@snr_files) {
    my $cmd = "parse.snr $raw_dir/$snr_file > $output_dir/$snr_file.pdat";
    `$cmd`;
}

#####
## functions
