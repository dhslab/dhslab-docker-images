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

 TranscriptCoordinates

=head1 SYNOPSIS

 mv TranscriptCoordinates.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TranscriptCoordinates

=head1 DESCRIPTION

 A VEP plugin that adds the genomic start and end coordinates of the
 transcript to the output for each transcript consequence.

=cut

package TranscriptCoordinates;

# Using the exact structure from the working TSSDistance.pm model
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

# Constructor, calls the constructor of the parent class
sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

# Defines the new fields to be added to the VEP output header
sub get_header_info {
    return {
        'T_START' => 'Start coordinate of the transcript',
        'T_END'   => 'End coordinate of the transcript',
    };
}

# Defines that this plugin works on Transcript features
sub feature_types {
    return ['Transcript'];
}

# Defines that this plugin also considers the variant feature
# This ensures the 'run' method receives the TranscriptVariationAllele object
sub variant_feature_types {
    return ['BaseVariationFeature'];
}

# Main plugin method
sub run {
    my ($self, $tva) = @_;

    # $tva is a Bio::EnsEMBL::Variation::TranscriptVariationAllele object
    # Get the transcript object from the $tva
    my $transcript = $tva->transcript;

    # Return the start and end coordinates if a transcript is present
    if ($transcript) {
        return {
            'T_START' => $transcript->start,
            'T_END'   => $transcript->end,
        };
    }

    # Otherwise return an empty hash
    return {};
}

1;