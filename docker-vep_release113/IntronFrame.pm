=head1 LICENSE

Copyright (c) 2025 David H. Spencer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=head1 CONTACT

 David H. Spencer

=cut

=head1 NAME

 IntronFrame

=head1 SYNOPSIS

 mv IntronFrame.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin IntronFrame

=head1 DESCRIPTION

 A VEP plugin that gets the coding frame of introns relative
 to a variant position and the transcript strand.

=cut

package IntronFrame;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub feature_types {
    return ['Transcript'];
}

sub variant_feature_types {
    return ['BaseVariationFeature'];
}

sub get_header_info {
    return {
        IntronFrame => "Phase of upstream and downstream splice junctions relative to variant position(s)",
    };
}

sub run {
    my ($self, $tva) = @_;

    my $tr = $tva->transcript;
    return {} unless $tr;

    my $vf;
    my $is_sv = 0;
    if ($tva->can('variation_feature')) {
        $vf = $tva->variation_feature;
    } elsif ($tva->can('structural_variation_feature')) {
        $vf = $tva->structural_variation_feature;
        $is_sv = 1;
    } else {
        return {};
    }

    return {} unless $vf;

    if ($is_sv) {
        my $type = '';
        if ($vf->can('class_SO_term')) {
            $type = $vf->class_SO_term;
        }
        
        if ($type =~ /deletion/i || $type =~ /duplication/i || $type =~ /copy_number/i) {

            # If variant start < transcript start and variant end > transcript end
            if ($vf->start < $tr->start && $vf->end > $tr->end) {
                return {};
            }

            my ($start_phase, $start_feature) = $self->_get_feature_info($tr, $vf->start);
            my ($end_phase, $end_feature) = $self->_get_feature_info($tr, $vf->end);

            # If start and end are in the same feature (and feature is found) dont return phase info
            if (defined($start_feature) && defined($end_feature) && $start_feature eq $end_feature) {
                return {};
            }

            # Otherwise, report both if they are different. This means the breakpoints are in an exon.
            my @results;
            push @results, $start_phase if $start_phase ne '';
            push @results, $end_phase if $end_phase ne '';
            
            my %seen;
            my @unique_results = grep { !$seen{$_}++ } @results;

            return { IntronFrame => join('&', @unique_results) } if @unique_results;
            return {};
        }
    }

    # Fallback for small variants and other SVs (BND, INS)
    my @positions = ($vf->start);

    my @results;
    foreach my $pos (@positions) {
        my $res = $self->get_phase($tr, $pos);
        push @results, $res if $res ne '';
    }
    
    return {} unless @results;

    return { IntronFrame => join('&', @results) };
}

sub get_phase {
    my ($self, $tr, $pos) = @_;
    my ($phase, undef) = $self->_get_feature_info($tr, $pos);
    return $phase;
}

sub _get_feature_info {
    my ($self, $tr, $pos) = @_;

    # Get exons in 5' to 3' order (transcription order)
    my @exons = @{$tr->get_all_Exons};
    my $strand = $tr->strand;

    # Check if pos is in an intron
    for my $i (0 .. $#exons - 1) {
        my $e1 = $exons[$i];
        my $e2 = $exons[$i+1];
        
        my $in_intron = 0;
        if ($strand == 1) {
            $in_intron = 1 if ($pos > $e1->end && $pos < $e2->start);
        } else {
            $in_intron = 1 if ($pos < $e1->start && $pos > $e2->end);
        }

        if ($in_intron) {
             # Upstream exon is e1 (5'). Downstream is e2 (3').
             # We want e1 downstream phase, e2 upstream phase.
             my $up   = defined($e1->end_phase) ? $e1->end_phase : -1;
             my $down = defined($e2->phase) ? $e2->phase : -1;
             return ("$up", "intron$i");
        }
    }

    return ("", undef);
}

1;
