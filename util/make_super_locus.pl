#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Overlap_piler;


my $usage = "usage: $0 all_exons.gtf gene_name gene_id transcript_id\n\n";

my $exons_gtf_file = $ARGV[0] or die $usage;
my $gene_name = $ARGV[1] or die $usage;
my $gene_id = $ARGV[2] or die $usage;
my $trans_id = $ARGV[3] or die $usage;


my $ORIENT;
my $CHR;

my @coords;

open (my $fh, $exons_gtf_file) or die $!;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    if ($x[2] eq "exon") {
        my ($chr, $lend, $rend, $orient) = ($x[0], $x[3], $x[4], $x[6]);
        if ($CHR) {
            if ($chr ne $CHR || $orient ne $ORIENT) {
                die "Error, inconsistent chr or orient info";
            }
        }
        else {
            $CHR = $chr;
            $ORIENT = $orient;

        }
        push (@coords, [$lend, $rend]);
    }
}
close $fh;

my @collapsed_coords = &Overlap_piler::simple_coordsets_collapser(@coords);

foreach my $coordset (@collapsed_coords) {
    my ($lend, $rend) = @$coordset;

    print join("\t", $CHR, "SuperLocus", "exon", $lend, $rend, ".", $ORIENT, ".",
               "gene_id \"$gene_id\"; transcript_id \"$trans_id\"; gene_name \"$gene_name\";") . "\n";
}

exit(0);


    
