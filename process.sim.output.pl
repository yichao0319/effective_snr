###############################################
# Yi-Chao Chen @ UT Austin CS
#
# e.g.
#   perl process.sim.output.pl tmp_rx_syms_plane.dat.out
#

#!/bin/perl



use strict;

#####
## constant
my $DEBUG0 = 0;
my $DEBUG1 = 0;
my $NUM_MODULATIONS = 4;
my @TABLE_TPUT = (1, 2, 4, 6);

#####
## variables
if(@ARGV != 1) {
    die "wrong number of input: ".scalar(@ARGV)."\n";
}
my $filename = $ARGV[0];
print $filename."\n" if($DEBUG1);

my $num_schemes = 10;
## info{scheme1 | scheme2 | ...}{pkt_i}{throughput | modulation | BPSK_BER}
my %pkt_info = ();

my $num_pkts = 0;
## $frame_err_rate{scheme_i}{modulation_i}{frame_err_rate}
my %frame_err_rate = ();
## scheme_info{scheme1 | scheme2 | ...}{throughputs | modulations | num_under_selection | num_over_selection}
my %scheme_info = ();
## $frame_prr2ber{scheme_i}{modulation_i}{good | bad}{BERs}
my %frame_prr2ber = ();




#####
## main

my $input_dir = './OUTPUT_sim';
open FH, "<$input_dir/$filename" or dir $!;
<FH>;
<FH>;
while(<FH>) {
    print $_ if($DEBUG1);


    #####
    ## read file
    my @tmp = split(/[:,;]/, $_);
    print join(", ", @tmp) if($DEBUG0);

    $num_pkts = $tmp[0];
    my @actual_BERs = @tmp[1 .. 4];
    print ">".join(",", @actual_BERs)."\n" if($DEBUG0);
    my @scheme_selections = @tmp[5, 6, 15, 16, 25, 26, 35, 36, 45, 46];
    print ">".join(",", @scheme_selections)."\n" if($DEBUG0);
    my @effsnr_formula_BERs = @tmp[7 .. 10];
    my @effsnr_formula_PRRs = @tmp[11 .. 14];
    my @effsnr_threshold_BERs = @tmp[17 .. 20];
    my @effsnr_threshold_PRRs = @tmp[21 .. 24];
    my @entiresnr_formula_BERs = @tmp[27 .. 30];
    my @entiresnr_formula_PRRs = @tmp[31 .. 34];
    my @entiresnr_threshold_BERs = @tmp[37 .. 40];
    my @entiresnr_threshold_PRRs = @tmp[41 .. 44];
    print ">".join(",", @entiresnr_threshold_BERs)."\n" if($DEBUG0);


    #####
    ## 1) frame error rate
    foreach my $mod_i (0 .. $NUM_MODULATIONS-1) {
        ##   i) actual frame error rate
        $frame_err_rate{'1oracle'}{$mod_i} ++ if($actual_BERs[$mod_i] > 0);
        ##   ii) estimated frame error rate
        $frame_err_rate{'2effsnr_formula'}{$mod_i} += (1 - $effsnr_formula_PRRs[$mod_i]);
        $frame_err_rate{'3effsnr_threshold'}{$mod_i} += (1 - $effsnr_threshold_PRRs[$mod_i]);
        $frame_err_rate{'4entiresnr_formula'}{$mod_i} += (1 - $entiresnr_formula_PRRs[$mod_i]);
        $frame_err_rate{'5entiresnr_threshold'}{$mod_i} += (1 - $entiresnr_threshold_PRRs[$mod_i]);
    }


    #####
    ## 2) throughput
    ##   i) best throughput
    my $selected_mod = $scheme_selections[0] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'01oracle_curr'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##   ii) throughput of other rate adaptation schemes
    ##       oracle with prediction
    $selected_mod = $scheme_selections[1] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'02oracle_pred'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       eff SNR using formula with current pkt
    $selected_mod = $scheme_selections[2] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'03effsnr_formula_curr'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       eff SNR using formula with prediction
    $selected_mod = $scheme_selections[3] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'04effsnr_formula_pred'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       eff SNR using threshold with current pkt
    $selected_mod = $scheme_selections[4] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'05effsnr_threshold_curr'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       eff SNR using threshold with prediction
    $selected_mod = $scheme_selections[5] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'06effsnr_threshold_pred'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       entire SNR using formula with current pkt
    $selected_mod = $scheme_selections[6] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'07entiresnr_formula_curr'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       entire SNR using formula with prediction
    $selected_mod = $scheme_selections[7] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'08entiresnr_formula_pred'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       entire SNR using threshold with current pkt
    $selected_mod = $scheme_selections[8] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'09entiresnr_threshold_curr'}{throughput} += $TABLE_TPUT[$selected_mod];
    }
    ##       entire SNR using threshold with prediction
    $selected_mod = $scheme_selections[9] - 1;
    if($actual_BERs[$selected_mod] == 0) {
        $scheme_info{'10entiresnr_threshold_pred'}{throughput} += $TABLE_TPUT[$selected_mod];
    }


    #####
    ## 3) find out frame PRR to BER mapping
    foreach my $mod_i (0 .. $NUM_MODULATIONS-1) {
        if($actual_BERs[$mod_i] == 0) {
            ## correct pkt
            push(@{$frame_prr2ber{'2effsnr_formula'}{$mod_i}{good}}, $effsnr_formula_BERs[$mod_i]);
            push(@{$frame_prr2ber{'3effsnr_threshold'}{$mod_i}{good}}, $effsnr_threshold_BERs[$mod_i]);
            push(@{$frame_prr2ber{'4entiresnr_formula'}{$mod_i}{good}}, $entiresnr_formula_BERs[$mod_i]);
            push(@{$frame_prr2ber{'5entiresnr_threshold'}{$mod_i}{good}}, $entiresnr_threshold_BERs[$mod_i]);
        }
        else {
            ## corrupted
            push(@{$frame_prr2ber{'2effsnr_formula'}{$mod_i}{bad}}, $effsnr_formula_BERs[$mod_i]);
            push(@{$frame_prr2ber{'3effsnr_threshold'}{$mod_i}{bad}}, $effsnr_threshold_BERs[$mod_i]);
            push(@{$frame_prr2ber{'4entiresnr_formula'}{$mod_i}{bad}}, $entiresnr_formula_BERs[$mod_i]);
            push(@{$frame_prr2ber{'5entiresnr_threshold'}{$mod_i}{bad}}, $entiresnr_threshold_BERs[$mod_i]);
        }
        ##   ii) estimated frame error rate
        
    }

}
close FH;


#####
## 1) frame error rate
foreach my $scheme (sort {$a <=> $b} keys(%frame_err_rate)) {
    foreach my $mod_i (sort {$a <=> $b} (keys %{$frame_err_rate{$scheme}})) {
        $frame_err_rate{$scheme}{$mod_i} /= $num_pkts;
    }
}
print_frame_err_rate(\%frame_err_rate) if($DEBUG1);
print_frame_err_rate2(\%frame_err_rate);



#####
## 2) throughput
print_throughput(\%scheme_info) if($DEBUG1);
print_throughput2(\%scheme_info);



#####
## 3) find out frame PRR to BER mapping
print_frame_prr2ber(\%frame_prr2ber) if($DEBUG1);

1;


#####
## functions

sub print_frame_err_rate {
    my ($ref_frame_err_rate) = @_;

    print "\nframe error rate: \n";
    foreach my $scheme (sort {$a <=> $b} keys(%{$ref_frame_err_rate})) {
        print "  $scheme: \n";
        foreach my $mod_i (sort {$a <=> $b} (keys %{$ref_frame_err_rate->{$scheme}})) {
            print "    modulation $mod_i = ".$ref_frame_err_rate->{$scheme}{$mod_i}."\n";
        }
    }
}

sub print_frame_err_rate2 {
    my ($ref_frame_err_rate) = @_;

    print "\nframe error rate: \n";
    foreach my $scheme (sort {$a <=> $b} keys(%{$ref_frame_err_rate})) {
        print "$scheme, ";
        foreach my $mod_i (sort {$a <=> $b} (keys %{$ref_frame_err_rate->{$scheme}})) {
            print $ref_frame_err_rate->{$scheme}{$mod_i}.", ";
        }
        print "\n";
    }
}


sub print_throughput {
    my ($ref_scheme_info) = @_;

    print "\nthroughput: \n";
    foreach my $scheme (sort {$a <=> $b} keys(%{$ref_scheme_info})) {
        print "  $scheme = ".$ref_scheme_info->{$scheme}{throughput}." \n";
    }
}

sub print_throughput2 {
    my ($ref_scheme_info) = @_;

    print "\nthroughput: \n";
    my $first = -1;
    foreach my $scheme (sort {$a <=> $b} keys(%{$ref_scheme_info})) {
        if($first == -1) {
            $first = $ref_scheme_info->{$scheme}{throughput};
        }
        print "".($ref_scheme_info->{$scheme}{throughput} / $first).", ";
    }
    print "\n";
}


sub print_frame_prr2ber {
    my ($ref_frame_ppr2ber) = @_;

    print "\nframe PPR to BER mapping: \n";
    foreach my $scheme (sort {$a <=> $b} keys(%{$ref_frame_ppr2ber})) {
        print "  $scheme: \n";
        foreach my $mod_i (sort {$a <=> $b} (keys %{$ref_frame_ppr2ber->{$scheme}})) {
            my $good_avg = "n/a";
            my $good_stdev = "n/a";
            my $bad_avg = "n/a";
            my $bad_stdev = "n/a";

            if(defined(@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{good}})) {
                $good_avg = average(\@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{good}});
                $good_stdev = stdev(\@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{good}});
            }
            
            if(defined(@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{bad}})) {
                $bad_avg = average(\@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{bad}});
                $bad_stdev = stdev(\@{$ref_frame_ppr2ber->{$scheme}{$mod_i}{bad}});
            }

            print "    modulation $mod_i: good ($good_avg, $good_stdev); bad ($bad_avg, $bad_stdev)\n";
        }
    }
}




sub average{
    my($data) = @_;
    if (not @$data) {
        return 0;
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}

sub stdev{
    my($data) = @_;
    if(@$data == 1){
        return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach(@$data) {
        $sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / (@$data-1)) ** 0.5;
    return $std;
}
