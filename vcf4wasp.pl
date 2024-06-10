#!/usr/bin/perl -w
use strict;

my $snv_file = shift or die "Usage: perl $0 <variant_file>\n";
die "Input file is not readable" unless -r $snv_file;

open( F, $snv_file ) or die 'Cannot open infile';
while ( <F> ) {
    if ( m/^\#/ ) {
        print;
        next;
    }
    chomp;
    my ( $chr, $pos, $rs, $ref, $alt, $varq, $filter, $info, $format ) = split /\t/;
    next unless $ref =~ m/^[GATC]$/ and $alt =~ m/^[GATC]$/; #Only SNVs
    my ($ac) = $info =~ m/AC\=(\d+)/;
    my ($an) = $info =~ m/AN\=(\d+)/;
    my ($dp4) = $info =~ m/DP4\=([\d\,]+)/;
    my ( $ref_pos, $ref_neg, $alt_pos, $alt_neg ) = split /\,/, $dp4;
    my $total = $ref_pos + $ref_neg + $alt_pos + $alt_neg;
    my $af = ( $alt_pos+$alt_neg) / $total; 
    warn join( "\t", $pos, $dp4, $af ), "\n";
    next if $af < 0.01 or $af > 0.99 ; # Skip likely monomorphic sites
    print join( "\t", $chr, $pos, $rs, $ref, $alt, '.', '.', '.', 'GT', '0|1' ), "\n";
}
close F;

