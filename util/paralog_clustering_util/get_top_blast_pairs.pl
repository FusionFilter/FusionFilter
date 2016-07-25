#!/usr/bin/env perl

use strict;
use warnings;



main: {
    

    # m blastn.outfmt6.grouped.geneSym | sort -k4,4g -k3,3gr > blastn.outfmt6.grouped.geneSym.sorted


    my %data;
    
    open (my $fh, "blastn.outfmt6.grouped.geneSym.sorted") or die $!;
    while (<$fh>) {
        my $line = $_;
        chomp;
        my @x = split(/\t/);
        my $geneA = $x[0];
        my $geneB = $x[1];
        
        my $token = join("$;", sort ($geneA, $geneB) );

        unless ($data{$token}) {
            $data{$token} = 1;
            print $line;
        }

    }

    exit(0);
}



        
