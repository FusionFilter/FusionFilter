#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;


my $usage = "\n\n\tusage: $0 blastn.gene_pairs.gz output_index_filename\n\n";

my $blast_pairs = $ARGV[0] or die $usage;
my $output_index_filename = $ARGV[1] or die $usage;

main: {

    my $idx = new TiedHash( { create => $output_index_filename } );
    
    open (my $fh, "gunzip -c $blast_pairs | ") or die "Error, cannot read $blast_pairs  ";
    while (<$fh>) {
        chomp;
        my ($geneA, $geneB, $per_id, $Evalue) = split(/\s+/);
        if ($geneA eq $geneB) { next; }
        
        my $token = "{$geneA|$geneB|pID:$per_id|E:$Evalue}";
        
        $idx->store_key_value("$geneA--$geneB", $token);
        $idx->store_key_value("$geneB--$geneA", $token);
                
    }
    
    print STDERR "-done creating index file: $output_index_filename\n\n";

    exit(0);
    
}

