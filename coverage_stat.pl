#! /usr/bin/perl

use strict;
use warnings;

use feature qw(say);

my ($normal_IDT_stat, $normal_RMG_stat, $normal_CMR_stat, $tumor_IDT_stat, $tumor_RMG_stat, $tumor_CMR_stat, $followup_IDT_stat, $followup_RMG_stat, $followup_CMR_stat, $out_file) = @ARGV;

my @metrics_names = qw(
    MEAN_TARGET_COVERAGE     
    MEDIAN_TARGET_COVERAGE   
    ON_TARGET_BASES  
    ZERO_CVG_TARGETS_PCT     
    PCT_USABLE_BASES_ON_TARGET   
    PCT_TARGET_BASES_100X    
    PCT_TARGET_BASES_30X     
    PCT_TARGET_BASES_10X     
    PCT_TARGET_BASES_1X  
    PCT_EXC_BASEQ    
    PCT_EXC_DUPE     
    PCT_EXC_MAPQ     
    PCT_EXC_OFF_TARGET   
    PCT_EXC_OVERLAP  
);

my %stats = (
    normal   => {
        IDT => $normal_IDT_stat,
        RMG => $normal_RMG_stat,
        CMR => $normal_CMR_stat,
    },
    tumor    => {
        IDT => $tumor_IDT_stat,
        RMG => $tumor_RMG_stat,
        CMR => $tumor_CMR_stat,
    },
    followup => {
        IDT => $followup_IDT_stat,
        RMG => $followup_RMG_stat,
        CMR => $followup_CMR_stat,
    },
);

open (my $out_fh, ">", "$out_file")
    or die "couldn't open $out_file to write";

my $header = join "\t", 'Sample', 'ROI_set', @metrics_names;
say $out_fh $header;

for my $type (sort keys %stats) {
    for my $roi (sort keys %{$stats{$type}}) {
        open (my $stat_fh, $stats{$type}->{$roi}) or die "failed to open $type $roi coverage_stat";
        my (@metrics_keys, @metrics_values);
    
        while (<$stat_fh>) {
            chomp;
            next if /^(#|\s+)/;
            @metrics_keys = split /\t/, $_ if /^BAIT_SET/;
            @metrics_values = split /\t/, $_ if /^(AML|xgen)/;
            last if @metrics_keys and @metrics_values;
        }
        close $stat_fh;

        my %metrics;
        @metrics{@metrics_keys} = @metrics_values;
    
        my @values;
        for my $metrics_name (@metrics_names) {
            push @values, $metrics{$metrics_name};
        }
    
        my $values = join "\t", $type, $roi, @values;
        say $out_fh $values;
    }
}

close $out_fh;

