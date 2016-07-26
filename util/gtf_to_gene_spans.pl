#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 annots.gtf\n\n";

my $annots_gtf = $ARGV[0] or die $usage;

main: {

    my %data;

    open (my $fh, $annots_gtf) or die "Error, cannot open file $annots_gtf";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        unless (/\w/) { next; }
        my @x = split(/\t/);
        my $chr = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        
        my $info = $x[8];
        $info =~ /gene_id \"([^\"]+)\"/ or die "Error, cannot extract gene_id from $info";
        
        my $gene_id = $1;

        my $gene_name = "";
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_name = $1;
        }
        
        

        push (@{$data{$gene_id}->{coords}}, $lend, $rend);
        $data{$gene_id}->{chr} = $chr;
        $data{$gene_id}->{orient} = $orient;
        $data{$gene_id}->{gene_name} = $gene_name;
        
    }
    close $fh;

    foreach my $gene (keys %data) {
        my $chr = $data{$gene}->{chr};
        my $orient = $data{$gene}->{orient};

        my @coords = @{$data{$gene}->{coords}};

        @coords = sort {$a<=>$b} @coords;

        my $lend = shift @coords;
        my $rend = pop @coords;
        my $gene_name = $data{$gene}->{gene_name} || ".";

        print join("\t", $gene, $chr, $lend, $rend, $orient, $gene_name) . "\n";
    }

    exit(0);
}

