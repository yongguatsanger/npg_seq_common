# Author:        John O'Brien and Marina Gourtovaia
# Created:       26 January 2010
#

package npg_common::sequence::reference::base_count;

use Moose;
use MooseX::ClassAttribute;
use Carp;
use IO::File;
use Fatal qw( open close );
use MooseX::Storage;
use Readonly;

with Storage( 'traits' => ['OnlyWhenBuilt'],
              'format' => 'JSON',
              'io'     => 'File' );

our $VERSION = '0';

## no critic (Documentation::RequirePodAtEnd ProhibitParensWithBuiltins)

=head1 NAME

npg_common::sequence::reference::base_count

=head1 VERSION

=head1 SYNOPSIS

 my $bc = npg_common::data::reference::base_count->new(reference_path => q[my.ref]);
 $bc->run;
 
=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=cut


Readonly::Scalar my $SEQ_NAME_LINE_TAG  => 62;

Readonly::Array  my @IUPAC_BASE_CODES     => qw(A C G T M R W S Y K V H D B N);
Readonly::Scalar my $GAP_CHAR => q[_];
Readonly::Scalar my $GAP_CHAR_1 => q[.];
Readonly::Scalar my $N_EQUIVALENT => q[X];
Readonly::Scalar my $COUNTS_HASH  => q[counts];
Readonly::Scalar my $NA           => 1;


=head2 output_base_codes

Class attribute. Returns a ref to a list of output base codes

=cut
class_has '_output_base_codes' => (isa             => 'ArrayRef',
                                   is              => 'ro',
                                   required        => 0,
                                   default         => sub {return \@IUPAC_BASE_CODES;},
                                   reader          => 'output_base_codes',
);

=head2 translation_table

Class attribute. Returns a ref to a hash with a translation table from ASCII characters
to base output codes.

=cut
class_has '_translation_table' => (isa             => 'HashRef',
                                   is              => 'ro',
                                   required        => 0,
                                   reader          => 'translation_table',
                                   lazy_build      => 1,
);
sub _build__translation_table {
    my $hash = {};
    my @codes = ();

    foreach my $code (@IUPAC_BASE_CODES) {
        $hash->{ord($code)}    = $code;
        $hash->{ord(lc $code)} = $code;
    }
    $hash->{ord($N_EQUIVALENT)}    = q[N];
    $hash->{ord(lc $N_EQUIVALENT)} = q[N];

    $hash->{ord($GAP_CHAR)}   = q[];
    $hash->{ord($GAP_CHAR_1)} = q[];

    return $hash;
}

=head2 reference_path

Reference sequence file path

=cut
has 'reference_path' => ( isa           => 'Str',
                          is            => 'ro',
                          required      => 1,
                          documentation => 'Path to a reference fasta file',
);

=head2 _summary

A private variable with a public reader. A reference to a hash with information
about a sequence as a whole.

=cut
has '_summary'     => ( isa           => 'HashRef',
                        is            => 'rw',
                        required      => 0,
                        reader        => 'summary',
                        default       => sub { my $h={'ref_length' => 0, 'counts' => {}}; return $h;},
);

=head2 run

Performs the base count. The results can be accessed through the bin and summary methods.

=cut
sub run {
    my ($self) = @_;
    if ( !-r $self->reference_path) {
        croak $self->reference_path . q[does not exist or is not readable];
    }
    my $fh = IO::File->new( $self->reference_path(), q[<] );
    $self->_parse_reference($fh);
    $fh->close();
    return;
}

sub _parse_reference {
    my ($self, $fh) = @_;

    while (my $line = <$fh>) {
        $line =~ s/\s//gmsx;
        if (!$line) { next; }
        # Convert a line to an array of ascii values for all characters in a string
        my @ascii_values = unpack('C*', $line);

        # Trap header lines.
        if ($ascii_values[0] != $SEQ_NAME_LINE_TAG) {
            foreach my $ascii (@ascii_values) {
                if (!exists $self->translation_table->{$ascii}) {
                    croak "$line contains unexpected character $ascii";
		}
                my $char_out = $self->translation_table->{$ascii};
                if ($char_out) { #if this char is not the gap char
                    $self->_summary->{$COUNTS_HASH}->{$char_out}++;
	        }
                $self->_summary->{ref_length}++;
            }
        }
    }
    return;
}


__PACKAGE__->meta->make_immutable;
no Moose;

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item MooseX::ClassAttribute

=item MooseX::Storage

=item Carp

=item IO::File

=item Fatal qw( open close )

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Marina Gourtovaia E<lt>mg8@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2013 GRL, by Marina Gourtovaia

This file is part of NPG.

NPG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
