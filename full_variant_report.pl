#!/usr/bin/perl

use strict;
use warnings;

use feature qw(say);

die("Wrong number of arguments. Provide variant tsv, followup snv bam_readcount, followup indel bam_readcount, pindel region vcf, output tsv") unless @ARGV == 5;
my ($tsv, $followup_snv_bamrc, $followup_indel_bamrc, $pindel_region_vcf, $output_tsv) = @ARGV;

my $snv_bamrc_info   = get_bamrc_info($followup_snv_bamrc);
my $indel_bamrc_info = get_bamrc_info($followup_indel_bamrc);

open(my $tsv_fh, $tsv) or die("couldn't open $tsv for read");
open(my $out_fh, ">", "$output_tsv") or die("couldn't open $output_tsv for write");

chomp (my $header = `/usr/bin/head -1 $tsv`);
$header =~ s/\r//g;
my @headers = split /\t/, $header;
push @headers, qw(FOLLOWUP.AD FOLLOWUP.DP FOLLOWUP.AF);
say $out_fh join "\t", @headers;

while (<$tsv_fh>) {
    chomp;
    next if /^CHROM/;
    $_ =~ s/\r//g;
    my %info;
    @info{@headers} = split /\t/, $_;;
    my @alts = split /,/, $info{ALT};
    if (length($info{REF}) == 1 and length($alts[0]) == 1) { #SNV
        my $id = $info{CHROM}.'__'.$info{POS};
        my $followup_bamrc_info = $snv_bamrc_info->{$id};
        my $followup_AD = $followup_bamrc_info->{allele_ct}->{$info{REF}}.','.$followup_bamrc_info->{allele_ct}->{$alts[0]};
        my $followup_DP = $followup_bamrc_info->{depth};
        my $followup_AF = sprintf("%.5f", $followup_bamrc_info->{allele_ct}->{$alts[0]}/$followup_DP);
        say $out_fh $_."\t$followup_AD\t$followup_DP\t$followup_AF";
    }
    else { #INDEL
        my ($pos, $bamrc_str);
        if (length($info{REF}) > length($alts[0])) { #DEL
            $pos = $info{POS} + 1;
            $bamrc_str = '-'.substr($info{REF}, length($alts[0]));
        }
        else { #INS
            $pos = $info{POS};
            $bamrc_str = '+'.substr($alts[0], length($info{REF}));
        }
        my $id = $info{CHROM}.'__'.$pos;
        my $followup_bamrc_info = $indel_bamrc_info->{$id};
        unless ($followup_bamrc_info) {
            warn "$id does not have indel bam_readcount info available";
            say $out_fh $_."\tNA\tNA\tNA";
            next;
        }
        my $alt_ct = exists $followup_bamrc_info->{allele_ct}->{$bamrc_str} ? $followup_bamrc_info->{allele_ct}->{$bamrc_str} : 0;
        my $followup_AD = $followup_bamrc_info->{allele_ct}->{$followup_bamrc_info->{ref}}.','.$alt_ct;
        my $followup_DP = $followup_bamrc_info->{depth};
        my $followup_AF = sprintf("%.5f", $alt_ct/$followup_DP);
        say $out_fh $_."\t$followup_AD\t$followup_DP\t$followup_AF";
    }
}
close $tsv_fh;

open(my $vcf_fh, "/bin/gunzip -c $pindel_region_vcf |") or die("couldn't open $pindel_region_vcf to read");

while (<$vcf_fh>) {
    chomp;
    next if /^#/;
    my ($chr, $pos, undef, $ref, $alt, undef, undef, undef, $format, $normal_info, $tumor_info) = split /\t/, $_;
    my @format_keys = split /:/, $format;
    my @normal_info = split /:/, $normal_info;
    my @tumor_info  = split /:/, $tumor_info;
        
    my (%normal_hash, %tumor_hash);
    @normal_hash{@format_keys} = @normal_info;
    @tumor_hash{@format_keys}  = @tumor_info;

    my (undef, $normal_DP, $normal_AF) = get_pindel_region_info($normal_hash{AD});
    my ($tumor_alt_ct, $tumor_DP,  $tumor_AF) = get_pindel_region_info($tumor_hash{AD});

    next unless $tumor_alt_ct >= 5 and $tumor_AF >= 0.01;

    my %tsv_info = (
        CHROM       => $chr,
        POS         => $pos,
        REF         => $ref,
        ALT         => $alt,
        set         => 'pindelRegion',
        SYMBOL      => 'FLT3',
        'NORMAL.GT' => $normal_hash{GT},
        'NORMAL.AD' => $normal_hash{AD},
        'NORMAL.DP' => $normal_DP,
        'NORMAL.AF' => $normal_AF,
        'TUMOR.GT'  => $tumor_hash{GT},
        'TUMOR.AD'  => $tumor_hash{AD},
        'TUMOR.DP'  => $tumor_DP,
        'TUMOR.AF'  => $tumor_AF,
    );

    my @contents;
    for my $tsv_header (@headers) {
        if (exists $tsv_info{$tsv_header}) {
            push @contents, $tsv_info{$tsv_header};
        }
        else {
            push @contents, '.';
        }
    }
    say $out_fh join "\t", @contents;
}

close $vcf_fh;
close $out_fh;


sub get_bamrc_info {
    my $bamrc = shift;
    my $bamrc_info;
    
    open(my $bamrc_fh, $bamrc) or die "failed to open $bamrc for read";
    while (<$bamrc_fh>) {
        chomp;
        my ($chr, $pos, $ref, $depth, $rest) = $_ =~ /^(.+)\t(\d+)\t([AaCcGgTtRrYyKkMmSsWwBbDdHhVvNn])\t(\d+)\t(.+)$/x; 
        my @allele_info = split /\t/, $rest;
        my $allele_ct;
        for my $allele_str (@allele_info) {
            my @info = split /:/, $allele_str;
            $allele_ct->{$info[0]} = $info[1];
        }
        $bamrc_info->{$chr.'__'.$pos} = {
            ref => $ref,
            depth => $depth,
            allele_ct => $allele_ct,
        };
    }
    close $bamrc_fh;
    return $bamrc_info;
}

sub get_pindel_region_info {
    my $AD_str = shift;
    my ($ref_ct, $alt_ct) = split /,/, $AD_str;
    my $DP = $ref_ct + $alt_ct;
    my $AF = sprintf("%.5f", $alt_ct/$DP);
    return ($alt_ct, $DP, $AF);
}
