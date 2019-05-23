#! /usr/bin/perl

use strict;
use warnings;

use feature qw(say);

die("Wrong number of arguments. Provide normal alignment metrics, tumor alignment metrics, followup alignment metrics and output file") unless @ARGV == 4;
my ($normal_stat, $tumor_stat, $followup_stat, $out_file) = @ARGV;

my @metrics_names = qw(
    TOTAL_READS
    PF_READS_ALIGNED
    PCT_PF_READS_ALIGNED
    PF_ALIGNED_BASES
    PF_MISMATCH_RATE
    READS_ALIGNED_IN_PAIRS
    PCT_READS_ALIGNED_IN_PAIRS
);

my %stats = (
    normal   => $normal_stat,
    tumor    => $tumor_stat,
    followup => $followup_stat,
);

open (my $out_fh, ">", "$out_file")
    or die "couldn't open $out_file to write";

my $header = join "\t", 'Sample', @metrics_names;
say $out_fh $header;

for my $type (sort keys %stats) {
    open (my $stat_fh, $stats{$type}) or die "failed to open $type alignment_stat";
    my (@metrics_keys, @metrics_values);
    
    while (<$stat_fh>) {
        chomp;
        next if /^(#|\s+)/;
        @metrics_keys = split /\t/, $_ if /^CATEGORY/;
        @metrics_values = split /\t/, $_ if /^PAIR/;
        last if @metrics_keys and @metrics_values;
    }
    close $stat_fh;

    my %metrics;
    @metrics{@metrics_keys} = @metrics_values;
    
    my @values;
    for my $metrics_name (@metrics_names) {
        push @values, $metrics{$metrics_name};
    }
    
    my $values = join "\t", $type, @values;
    say $out_fh $values;
}

close $out_fh;

