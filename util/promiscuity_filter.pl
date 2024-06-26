#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use Data::Dumper;
use List::Util qw(max);
use lib ("$FindBin::Bin/../lib");
use Gene_overlap_check;
use DelimParser;

my $MAX_PROMISCUITY = 10;
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
#  --genome_lib_dir <string>      CTAT genome lib directory
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
my $genome_lib_dir = "";

my $EXCLUDE_LOCI_OVERLAP_CHECK = 0;

&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,

              'max_promiscuity=i' => \$MAX_PROMISCUITY,

              'min_pct_dom_promiscuity=i' => \$MIN_PCT_DOM_PROM,
                   
              'genome_lib_dir=s' => \$genome_lib_dir,

              'debug|d' => \$DEBUG,
              
              'exclude_loci_overlap_check' => \$EXCLUDE_LOCI_OVERLAP_CHECK,

    );

if (@ARGV) {
    die "Error, dont recognize arguments: @ARGV";
}


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file && $genome_lib_dir) {
    die $usage;
}


main: {

    my $final_preds_file = "$fusion_preds_file.post_promisc_filter";
    open (my $final_ofh, ">$final_preds_file") or die "Error, cannot write to $final_preds_file";
    
    my $filter_info_file = "$fusion_preds_file.post_promisc_filter.info";
    open (my $filter_ofh, ">$filter_info_file") or die "Error, cannot write to $filter_info_file";
    
    my @fusions;
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";


    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_headers = $delim_parser->get_column_headers(); 


    my $final_ofh_writer = new DelimParser::Writer($final_ofh, "\t", \@column_headers);
    my $filter_ofh_writer = new DelimParser::Writer($filter_ofh, "\t", [@column_headers, "FilterReason"]);
    
        
    while (my $row = $delim_parser->get_row()) {
        my $fusion_name = $row->{'#FusionName'};

        my $J = 0;
        my $S = 0;

        if (defined $row->{JunctionReadCount}) { 
        
            $J = ( (!defined($row->{est_J})) || $row->{est_J} eq "NA") ? $row->{JunctionReadCount} : $row->{est_J};
            $S = ( (!defined($row->{est_S})) || $row->{est_S} eq "NA") ? $row->{SpanningFragCount} : $row->{est_S};
        }
        
        my $num_LR = $row->{num_LR} || $row->{num_reads} || 0;
        
        
        if ($J !~ /\d/) { $J = 0; }
        if ($S !~ /\d/) { $S = 0; }
        
        if ($num_LR eq "NA") {
            $num_LR = 0;
        }
        
        my ($geneA, $geneB) = split(/--/, $fusion_name);

        my $score = $num_LR + $J*4 + $S;

	#print STDERR "$fusion_name\tscore: $score\n";
        
        my $fusion = { fusion_name => $fusion_name,
                       
                       J => $J,
                       S => $S,
                       
                       geneA => $geneA,
                       geneB => $geneB,
                       
                       score => $score, 
                       
                       sum_JS => $num_LR + $J + $S,
                       
                       row => $row,
        };
    
        push (@fusions, $fusion); 
    }
    close $fh;
    
    if ($DEBUG) {
        print STDERR "Incoming fusions: " . Dumper(\@fusions);
    }
    
    
    # generate outputs
    
    @fusions = reverse sort {$a->{score} <=> $b->{score} 
                             ||
                                 $b->{fusion_name} cmp $a->{fusion_name}  # more stable sorting

    } @fusions;

    
    @fusions = &remove_promiscuous_fusions(\@fusions, $filter_ofh_writer, $MAX_PROMISCUITY, $MIN_PCT_DOM_PROM);
    

    foreach my $fusion (@fusions) {
        my ($geneA, $geneB) = split(/--/, $fusion->{fusion_name});

        my $row = $fusion->{row};

	$final_ofh_writer->write_row($row);
	
    }
    
    close $filter_ofh;
    close $final_ofh;
    
    

    exit(0);
}


####
sub remove_promiscuous_fusions {
    my ($fusions_aref, $filter_ofh_writer, $max_promiscuity, $min_pct_prom_dom) = @_;
            
    my @filtered_fusions = &filter_promiscuous_low_pct_prom_dom($fusions_aref, $filter_ofh_writer, $max_promiscuity, $min_pct_prom_dom);

    if ($DEBUG) {
        print STDERR "After filtering proms low pct dom: " . Dumper(\@filtered_fusions);
    }


    @filtered_fusions = &filter_remaining_promiscuous_fusions(\@filtered_fusions, $filter_ofh_writer, $max_promiscuity);
    
    
    if ($DEBUG) {
        print STDERR "After removing remaining promiscuous fusions: " . Dumper(\@filtered_fusions);
    }
    
    return(@filtered_fusions);
}

####
sub filter_promiscuous_low_pct_prom_dom {
    my ($fusions_aref, $filter_ofh_writer, $max_promiscuity, $min_pct_prom_dom) = @_;
    
    my %max_sum_support_fusion_partner = &get_max_sum_support_fusion_partner($fusions_aref);
    
    my %partners = &count_fusion_partners($fusions_aref);
    
    my @ret_fusions;
    
    foreach my $fusion (@$fusions_aref) {
        
        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};
        
        my $sum_JS = $fusion->{sum_JS};
	#my $score = $fusion->{score};
	        
        my $num_geneA_partners = scalar(keys %{$partners{$geneA}});
        my $num_geneB_partners = scalar(keys %{$partners{$geneB}});
        
        print STDERR "Fusion: $geneA--$geneB, partnersA: $num_geneA_partners, partnersB: $num_geneB_partners\n" if $DEBUG;
        
        if (&is_promiscuous($num_geneA_partners, $num_geneB_partners, $max_promiscuity)) {
            
            my $max_partner_support = max($max_sum_support_fusion_partner{$geneA}, $max_sum_support_fusion_partner{$geneB});

            my $pct_prom_dom = sprintf("%.1f", $sum_JS / $max_partner_support * 100);
            if ($pct_prom_dom < $min_pct_prom_dom) {

		my $row = $fusion->{row};
		$row->{FilterReason} = "FILTERED DUE TO reached max promiscuity ($max_promiscuity), num_partners($geneA)=$num_geneA_partners and num_partners($geneB)=$num_geneB_partners AND having only $pct_prom_dom support ($sum_JS) of max partner support ($max_partner_support)\n";

		$filter_ofh_writer->write_row($row);

		next; # important, don't keep it below.
            }
        }
        
        # ok, keeping it.
        push (@ret_fusions, $fusion);
        
    }
    
    return(@ret_fusions);
    
}

####
sub filter_remaining_promiscuous_fusions {
    my ($fusions_aref, $filter_ofh_writer, $max_promiscuity) = @_;
    
    print STDERR "-filter_remaining_promiscuous_fusions\n" if $DEBUG;
    
    ## any remaining promiscuous ones are to be tossed.
    
    my %partners = &count_fusion_partners($fusions_aref);

    my @ret_fusions;
    
    foreach my $fusion (@$fusions_aref) {

        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

        my @geneA_partners = keys %{$partners{$geneA}};
        my @geneB_partners = keys %{$partners{$geneB}};
        
        my $num_geneA_partners = scalar(@geneA_partners);
        my $num_geneB_partners = scalar(@geneB_partners);

	#my $score = $fusion->{score};
		
        print STDERR "Fusion: $geneA--$geneB, num_A_partners: $num_geneA_partners, num_B_partners: $num_geneB_partners\n" if $DEBUG;
        
        if ( (! $EXCLUDE_LOCI_OVERLAP_CHECK) # for testing
             && 
             &is_promiscuous($num_geneA_partners, $num_geneB_partners, $max_promiscuity)) {

            print STDERR "-looks potentially promiscuous, ressessing based on loci counting\n" if $DEBUG;

            ## reassess by counting chromosomal loci
            if ($num_geneA_partners > $max_promiscuity) {
                $num_geneA_partners = &reassess_loci_count(@geneA_partners);
            }
            if ($num_geneB_partners > $max_promiscuity) {
                $num_geneB_partners = &reassess_loci_count(@geneB_partners);
            }
            print STDERR "\treassessed: Fusion: $geneA--$geneB, num_A_partners: $num_geneA_partners, num_B_partners: $num_geneB_partners\n" if $DEBUG;


        }
    
        if (&is_promiscuous($num_geneA_partners, $num_geneB_partners, $max_promiscuity)) {

	    my $row = $fusion->{row};

	    $row->{FilterReason} = "FILTERED DUE TO reached max promiscuity ($max_promiscuity), num_partners($geneA)=$num_geneA_partners and num_partners($geneB)=$num_geneB_partners\n";

	    $filter_ofh_writer->write_row($row);
	    
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
        
        my $sum_support = $fusion->{sum_JS};
        unless ($sum_support > 0) {
            print STDERR "Warning! , no sum_JS for " . Dumper($fusion);
        }
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


####
sub reassess_loci_count {
    my (@genes) = @_;

    my %loci;
    my $gene_overlap_checker = new Gene_overlap_check("$genome_lib_dir/ref_annot.gtf.gene_spans");
    
    foreach my $gene (@genes) {
        my $gene_span_info_struct = $gene_overlap_checker->get_gene_span_info($gene);

        my $chr = $gene_span_info_struct->{chr};
        my $midpt = ($gene_span_info_struct->{lend} + $gene_span_info_struct->{rend})/2;

        my $locus_interval = int($midpt/1e5); # use 100k intervals

        my $locus_token = "$chr:$locus_interval";
        $loci{$locus_token}++;
    }

    my $num_loci = scalar(keys %loci);

    return($num_loci);
    
}
        
