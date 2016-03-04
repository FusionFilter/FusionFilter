#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use Data::Dumper;
use List::Util qw(max);

my $MAX_PROMISCUITY = 3;  # perhaps a poor choice of words, but still a best fit IMHO.
my $MIN_PCT_DOM_PROM = 20;

my $usage = <<__EOUSAGE__;


###################################################################################################
#
# Required:
#
#  --fusion_preds <string>        preliminary fusion predictions
#                                 Required formatting is:  
#                                 geneA--geneB (tab) junction_read_count (tab) spanning_read_count (tab) ... rest
#
# Optional: 
##
#
#  --max_promiscuity <int>               maximum number of partners allowed for a given fusion. Default: $MAX_PROMISCUITY
#
#  --min_pct_dom_promiscuity <int>       for promiscuous fusions, those with less than this support of the dominant scoring pair 
#                                        are filtered prior to applying the max_promiscuity filter.
#                                        (default: $MIN_PCT_DOM_PROM)
#
####################################################################################################

Two output files are generated:

    \${fusion_preds}.post_promisc_filter  : contains those fusion predictions excluding the promiscuous entries.

    \${fusion_preds}.post_promisc_filter.info   : contains all input fusions and promiscuous entries are annotated accordingly.



__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;
my $DEBUG;

&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,

              'max_promiscuity=i' => \$MAX_PROMISCUITY,

              'min_pct_dom_promiscuity=i' => \$MIN_PCT_DOM_PROM,
                   
              'debug' => \$DEBUG,
    );

if (@ARGV) {
    die "Error, dont recognize arguments: @ARGV";
}


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file) {
    die $usage;
}



=input_format:

0       ETV6--NTRK3
1       84
2       18
3       ONLY_REF_SPLICE
4       ETV6^ENSG00000139083.6
5       chr12:12022903:+
6       NTRK3^ENSG00000140538.12
7       chr15:88483984:-
8       comma-delim list of junction reads
9       comma-delim list of spanning frags

=cut


main: {

    my $final_preds_file = "$fusion_preds_file.post_promisc_filter";
    open (my $final_ofh, ">$final_preds_file") or die "Error, cannot write to $final_preds_file";
    
    my $filter_info_file = "$fusion_preds_file.post_promisc_filter.info";
    open (my $filter_ofh, ">$filter_info_file") or die "Error, cannot write to $filter_info_file";
    
    my @fusions;
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    unless ($header =~ /^\#/) {
        die "Error, file $fusion_preds_file doesn't begin with a header line";
    }
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my $line = $_;
        my @x = split(/\t/);
        my $fusion_name = $x[0];
        my $J = $x[1];
        my $S = $x[2];


        my ($geneA, $geneB) = split(/--/, $fusion_name);

        #my $score = sqrt($J**2 + $S**2);
        my $score = $J*4 + $S;
        
        
        my $fusion = { fusion_name => $fusion_name,
                       
                       J => $J,
                       S => $S,
                       
                       geneA => $geneA,
                       geneB => $geneB,
                       
                       score => $score, 
                       
                       sum_JS => $J + $S,
                       
                       line => $line,
        };
    
        push (@fusions, $fusion); 
    }
    close $fh;
    
    if ($DEBUG) {
        print STDERR "Incoming fusions: " . Dumper(\@fusions);
    }
    
    
    # generate outputs
    
    print $filter_ofh $header;
    print $final_ofh $header;
    
    @fusions = reverse sort {$a->{score} <=> $b->{score} 
                             ||
                                 $b->{fusion_name} cmp $a->{fusion_name}  # more stable sorting

    } @fusions;
    
    @fusions = &remove_promiscuous_fusions(\@fusions, $filter_ofh, $MAX_PROMISCUITY, $MIN_PCT_DOM_PROM);
    

    foreach my $fusion (@fusions) {
        my ($geneA, $geneB) = split(/--/, $fusion->{fusion_name});

        my $line = $fusion->{line};
        
        print $final_ofh "$line\n";
        print $filter_ofh "$line\n";
    }
    
    close $filter_ofh;
    close $final_ofh;
    
    

    exit(0);
}


####
sub remove_promiscuous_fusions {
    my ($fusions_aref, $filter_ofh, $max_promiscuity, $min_pct_prom_dom) = @_;
            
    my @filtered_fusions = &filter_promiscuous_low_pct_prom_dom($fusions_aref, $filter_ofh, $max_promiscuity, $min_pct_prom_dom);

    if ($DEBUG) {
        print STDERR "After filtering proms low pct dom: " . Dumper(\@filtered_fusions);
    }

    
    @filtered_fusions = &filter_remaining_promiscuous_fusions(\@filtered_fusions, $filter_ofh, $max_promiscuity);
    
    
    if ($DEBUG) {
        print STDERR "After removing remaining promiscuous fusions: " . Dumper(\@filtered_fusions);
    }
    
    return(@filtered_fusions);
}

####
sub filter_promiscuous_low_pct_prom_dom {
    my ($fusions_aref, $filter_ofh, $max_promiscuity, $min_pct_prom_dom) = @_;
    
    my %max_sum_support_fusion_partner = &get_max_sum_support_fusion_partner($fusions_aref);
    
    my %partners = &count_fusion_partners($fusions_aref);
    
    my @ret_fusions;
    
    foreach my $fusion (@$fusions_aref) {
        
        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};
        
        my $sum_JS = $fusion->{sum_JS};
        
        my $num_geneA_partners = scalar(keys %{$partners{$geneA}});
        my $num_geneB_partners = scalar(keys %{$partners{$geneB}});
        
        
        if (&is_promiscuous($num_geneA_partners, $num_geneB_partners, $max_promiscuity)) {
            
            my $max_partner_support = max($max_sum_support_fusion_partner{$geneA}, $max_sum_support_fusion_partner{$geneB});

            my $pct_prom_dom = sprintf("%.1f", $sum_JS / $max_partner_support * 100);
            if ($pct_prom_dom < $min_pct_prom_dom) {
                
                print $filter_ofh "#" . $fusion->{line} . "\tFILTERED DUE TO reached max promiscuity ($max_promiscuity), num_partners($geneA)=$num_geneA_partners and num_partners($geneB)=$num_geneB_partners AND having only $pct_prom_dom support ($sum_JS) of max partner support ($max_partner_support)\n";
                next;
            }
        }
        
        # ok, keeping it.
        push (@ret_fusions, $fusion);
        
    }
    
    return(@ret_fusions);
    
}

####
sub filter_remaining_promiscuous_fusions {
    my ($fusions_aref, $filter_ofh, $max_promiscuity) = @_;
    
    ## any remaining promiscuous ones are to be tossed.
    
    my %partners = &count_fusion_partners($fusions_aref);

    my @ret_fusions;
    
    foreach my $fusion (@$fusions_aref) {

        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

        my $num_geneA_partners = scalar(keys %{$partners{$geneA}});
        my $num_geneB_partners = scalar(keys %{$partners{$geneB}});
        
        if (&is_promiscuous($num_geneA_partners, $num_geneB_partners, $max_promiscuity)) {
            
            print $filter_ofh "#" . $fusion->{line} . "\tFILTERED DUE TO reached max promiscuity ($max_promiscuity), num_partners($geneA)=$num_geneA_partners and num_partners($geneB)=$num_geneB_partners\n";
        }
        else {
            push (@ret_fusions, $fusion);
        }
    }
    

    return(@ret_fusions);
}

####
sub count_fusion_partners {
    my ($fusions_aref) = @_;

    my %partners;
    
    foreach my $fusion (@$fusions_aref) {

        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

        $partners{$geneA}->{$geneB}++;
        $partners{$geneB}->{$geneA}++;

    }
    
    return(%partners);
}

####
sub get_max_sum_support_fusion_partner {
    my ($fusions_aref) = @_;

    my %partner_to_max_support;

    foreach my $fusion (@$fusions_aref) {
        my ($geneA, $geneB) = ($fusion->{geneA}, $fusion->{geneB});
        
        my $sum_support = $fusion->{sum_JS} or confess "Error, no sum_JS for " . Dumper($fusion);
        
        if ( (! exists $partner_to_max_support{$geneA}) || $partner_to_max_support{$geneA} < $sum_support) {
            $partner_to_max_support{$geneA} = $sum_support;
        }

        if ( (! exists $partner_to_max_support{$geneB}) || $partner_to_max_support{$geneB} < $sum_support) {
            $partner_to_max_support{$geneB} = $sum_support;
        }

    }

    return(%partner_to_max_support);
}


####
sub is_promiscuous {
    my ($num_geneA_partners, $num_geneB_partners, $max_promiscuity) = @_;
    

    if ($num_geneA_partners > $max_promiscuity || $num_geneB_partners > $max_promiscuity) {
        return(1);
    }
    else {
        return(0);
    }
}
