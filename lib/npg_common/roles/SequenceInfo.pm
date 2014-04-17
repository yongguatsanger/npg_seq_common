#############
# Created By: ajb
# Created On: 2011-06-14

package npg_common::roles::SequenceInfo;
use Moose::Role;
use Carp;

use Readonly;

our $VERSION = '0';

Readonly::Scalar our $FAILED_TO_FIND_INDEX  => -1;
Readonly::Scalar our $SANGER_QUALITY_OFFSET => 33;

## no critic (Documentation::RequirePodAtEnd)
=head1 NAME

npg_common::roles::SequenceInfo

=head1 VERSION

=head1 SYNOPSIS

  package MyPackage;
  use Moose;
  with qw{npg_common::roles::SequenceInfo};

  no Moose;
  1;

=head1 DESCRIPTION

A helper role for methods which can provide useful information on a sequence string

=head1 SUBROUTINES/METHODS

=head2 indices_of_base

passing in a dna string and a base, returns an array of the indices of those bases within the string

=cut

sub indices_of_base {
  my ( $self, $dna, $base, $case_sensitive ) = @_;

  if ( ! $case_sensitive ) {
    $dna = lc $dna;
    $base = lc $base;
  }

  my @base_positions;
  my $offset = 0;
  while ( ( my $idx_base = index $dna, $base, $offset ) != $FAILED_TO_FIND_INDEX ) {
    push @base_positions, $idx_base;
    $offset = $idx_base + 1;
  }

  return @base_positions;
}

=head2 convert_to_quality_score

Given an ascii character, converts it to the quality score it represents.
You can just provide the character, in which case it is assumed to use sanger offset of 33, or a hashref with both the quality_character and an offset (i.e. should you wish to use illumina offset of 64)

  my $iQualityScore = $class->convert_to_quality_score( '+' );

  my $iQualityScore = $class->convert_to_quality_score( {
    quality_character => 'B',
    offset => 64,
  } );

=cut

sub convert_to_quality_score {
  my ( $self, $args ) = @_;

  if ( ! defined $args ) {
    croak q{no details provided};
  }

  my $char;
  my $offset;
  if ( ! ref $args ) {
    $char = $args;
  }

  if ( ref $args ) {
    $char = $args->{quality_character};
    $offset = $args->{offset};
  }

  $offset ||= $SANGER_QUALITY_OFFSET;
  return ( ( ord $char ) - $offset );
}

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item English -no_match_vars

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Andy Brown

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2011 GRL, by Andy Brown (ajb@sanger.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
