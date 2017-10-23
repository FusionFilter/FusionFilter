#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;


my $usage = "\n\n\tusage: $0 ctat_genome_lib_dir\n\n";

my $ctat_genome_lib_dir = $ARGV[0] or die $usage;


main: {

    my $blast_idx = "$ctat_genome_lib_dir/blast_pairs.idx";
    unless (-s $blast_idx) {
        die "Error, cannot locate file: $blast_idx";
    }

    my $gene_spans_file = "$ctat_genome_lib_dir/ref_annot.gtf.gene_spans";
    unless (-s $gene_spans_file) {
        die "Error, cannot locate file: $gene_spans_file";
    }
        
    my $idx = new TiedHash( { 'use' => $blast_idx} );

    my %gene_to_span_info = &parse_gene_spans($gene_spans_file);

    my $counter = 0;
    foreach my $gene_pair ($idx->get_keys()) {
        $counter++;
        print STDERR "\r[$counter]   " if $counter % 100 == 0;
        my ($geneA, $geneB) = split(/--/, $gene_pair);
        if (&overlap($geneA, $geneB, \%gene_to_span_info)) {
            print STDERR "-must prune: $geneA--$geneB\n";
        }
    }

        
    print STDERR "-done updating blast index\n\n";
    
    exit(0);
    
}


####
sub parse_gene_spans {
    my ($gene_spans_file) = @_;


    my %gene_spans;
    
    open(my $fh, $gene_spans_file) or die "Error, cannot open file: $gene_spans_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($chr, $lend, $rend, $gene_id) = ($x[1], $x[2], $x[3], $x[5]);

        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        
        $gene_spans{$gene_id} = [$chr, $lend, $rend];
    }
    close $fh;

    return(%gene_spans);

}


####
sub overlap {
    my ($geneA, $geneB, $gene_to_span_info_href) = @_;

    my $geneA_info_aref = $gene_to_span_info_href->{$geneA};
    my ($chr_A, $lend_A, $rend_A) = @$geneA_info_aref;

    my $geneB_info_aref = $gene_to_span_info_href->{$geneB};
    my ($chr_B, $lend_B, $rend_B) = @$geneB_info_aref;

    if ($chr_A eq $chr_B
        &&
        $lend_A < $rend_B && $rend_A > $lend_B) {
        return(1);
    }
    else {
        return(0);
    }
}

