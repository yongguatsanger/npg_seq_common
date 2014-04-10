#############
# Created By: gq1
# Created On: 2010-03-17

package npg_common::roles::run::lane::tag_info;
use Moose::Role;

use strict;
use warnings;
use Carp qw(carp cluck croak confess);
use Perl6::Slurp;
use English qw{-no_match_vars};

with qw{npg_tracking::illumina::run::short_info npg_tracking::illumina::run::folder};
with qw{npg_common::roles::log};
with qw{npg_tracking::glossary::tag};

use Readonly;
our $VERSION = '0';
## no critic (RequirePodAtEnd)

Readonly::Array our @LANE_RANGE_ARRAY => ( 1..8 );
Readonly::Scalar our $LSF_INDEX_MULTIPLIER => 1000;

=head1 NAME

npg_common::roles::run::lane::tag_info

=head1 VERSION


=head1 SYNOPSIS

  package MyPackage;
  use Moose;
  with qw{npg_common::roles::run::lane::tag_info};

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 get_tag_index_list

read the lane tag file and return a list of tag index plus 0

=cut

has q{lane}       => ( isa => q{Int}, is => q{rw}, required => 0 );

sub get_tag_index_list {
  my ( $self, $lane ) = @_;

  $lane ||= $self->lane();

  if ( ! $lane ){
    croak 'no run lane number given';
  }

  my $tag_list_hash = $self->tag_list( $lane );

  my $tag_index_list = [ 0 ];

  foreach my $index ( sort keys %{ $tag_list_hash } ) {
    push @{ $tag_index_list }, $index;
  }

  return $tag_index_list;
}

=head2 tag_list

return a hashref, tag_index as key and tag as value

=cut

sub tag_list {
  my ( $self, $lane ) = @_;

  $lane ||= $self->lane();


  if(!$lane){
    croak 'no run lane number given';
  }

  my $lane_tag_file = $self->runfolder_path().qq{/lane_$lane.tag};
  my %tag_list;
  if(-e $lane_tag_file){
     my @lane_tag_lines = slurp $lane_tag_file;

     %tag_list = map {reverse split /\s+/mxs} @lane_tag_lines;

  }else{
    $self->log("The tag info file not available for lane $lane");
  }
  return \%tag_list;
}

=head2 is_multiplexed_lane

for a multiplexed run, returns 1 if this lane is multiplexed, 0 otherwise

=cut

sub is_multiplexed_lane {
  my ($self, $lane) = @_;

  $lane ||= $self->lane();

  if(!$lane){
    croak 'no run lane number given';
  }

  my $tag_index_list = $self->get_tag_index_list($lane);
  my $num_tag_index = scalar @{$tag_index_list};
  if( $num_tag_index > 1){
    return 1;
  }

  return 0;
}

=head2 all_lane_tag_pairs

store of lane_tag_pairs to be used

this is an array of hashrefs as follows

[ {
  position => x,
  tag_index => y,
},
...]

to populate this manually

  $oClass->add_lane_tag_pair( {
    position => $iX,
    tag_index => $iY,
  } );

or populate by utilising add_lane_tag_pairs/add_all_lane_tag_pairs

=cut

has q{lane_tag_pairs} => (
  traits     => [ q{Array} ],
  isa         => q{ArrayRef[HashRef[Int]]},
  is          => q{ro},
  init_arg    => undef,
  default    => sub { [] },
  handles    => {
    all_lane_tag_pairs    => q{elements},
    add_lane_tag_pair     => q{push},
    count_lane_tag_pairs  => q{count},
    has_no_lane_tag_pairs => q{is_empty},
    clear_lane_tag_pairs  => q{clear},
  },
);

=head2 add_lane_tag_pairs

works out and adds the lane_tag pairs to the arrayref lane_tag_pairs based on the lane tag file if it is a multiplexed lane
will croak if it can't find a lane position, either provided or via $self->lane()

  eval {
    $oClass->add_lane_tag_pairs();
  } or do {
    ...error handling in case you have not provided $self->lane()...
  };

  $oClass->add_lane_tag_pairs( $iLanePosition );

=cut

sub add_lane_tag_pairs {
  my ( $self, $lane ) = @_;

  $lane ||= $self->lane();

  if ( ! $lane ) {
    croak q{No lane provided. Please supply a lane position, or use add_all_lane_tag_pairs};
  }

  if ( $self->is_multiplexed_lane( $lane ) ) {
    my $tag_index_list = $self->get_tag_index_list( $lane );
    foreach my $tag_index ( @{ $tag_index_list } ) {

      # test that the tag index is within the allowed parameters
      # otherwise croak, as we must be able to deal with the index
      eval {
        $self->tag_index( $tag_index );
        1;
      } or do {
        croak $EVAL_ERROR;
      };

      $self->add_lane_tag_pair( {
        position => $lane,
        tag_index => $tag_index,
      } );
    }
  }

  return 1;
}

=head2 add_all_lane_tag_pairs

works out and adds the lane_tag pairs to the arrayref lane_tag_pairs based on the lane tag file if it is a multiplexed lane
for all the lanes.
If you have an arrayref $self->lanes(), it will use the contents of that, else loops through [ 1..8 ]

  $oClass->add_all_lane_tag_pairs();

For both add_lane_tag_pairs and add_all_lane_tag_pairs, they will croak if the find index values which are outside the allowed
parameters. For this, see role npg_tracking::glossary::tag

=cut

sub add_all_lane_tag_pairs {
  my ( $self ) = @_;

  my $lanes;
  if ( $self->can( q{lanes} ) && scalar @{ $self->lanes() } ) {
    $lanes = $self->lanes();
  } else {
    $lanes = \@LANE_RANGE_ARRAY;
  }

  foreach my $lane ( @{ $lanes } ) {
    $self->add_lane_tag_pairs( $lane );
  }

  return 1;
}

=head2 lsf_job_array_from_lane_tag_pairs

Generates an lsf_job_array from the lane_tag_pairs found. This does run add_all_lane_tag_pairs if no lane_tag_pairs already provided

The array will look as follows

  [1000-1006, 2000,2010,2100, 3000-3004, 3006, 5000, 6000, 7000, 7062-7083]

Which would represent
  Lanes 4 and 8 not being multiplexed,
  L1 with tags 1->6,
  L2 with tags 10 and 100,
  L3 with tags 1->4 and 6
  L5 and 6 multiplexed but with no tags apart from 0 (as highly unlikely as this is going to be!)
  L7 with tags 62->83

Should there actually be a tag with index 999

=cut

sub lsf_job_array_from_lane_tag_pairs {
  my ( $self ) = @_;

  if ( $self->has_no_lane_tag_pairs() ) {
    $self->add_all_lane_tag_pairs();
  }

  my @all_lane_tag_pairs = $self->all_lane_tag_pairs();

  my @lsf_indices;
  foreach my $pair ( @all_lane_tag_pairs ) {
    my $index = ( $pair->{position} * $LSF_INDEX_MULTIPLIER ) + $pair->{tag_index};
    push @lsf_indices, $index;
  }

  @lsf_indices = sort { $a <=> $b } @lsf_indices;

  return $self->create_array_string( @lsf_indices );
}

=head2 create_array_string

takes an array of integers, and then converts them to an LSF job array string for appending to a job_name

string takes the format

  [1,4-7,10...]

  my $sArrayString = $oClass->create_array_string( 1,4,5,6,7,10... );

=cut

sub create_array_string {
  my ( $self, @lsf_indices ) = @_;

  my ( $start_run, $end_run, $ret );

  $ret = q{};
  foreach my $entry ( @lsf_indices ) {

    # have we already started looping through
    if ( defined $end_run ) {

      # if the number is consecutive, increment end of the run
      if ( $entry == $end_run + 1 ) {
        $end_run = $entry;

      # otherwise, finish up that run, which may just be a single number
      } else {
        if ( $start_run != $end_run ) {
          $ret .= q{-} . $end_run;
        }
        $ret .= q{,} . $entry;
        $start_run = $end_run = $entry;
      }

    # we haven't looped through at least once, so set up
    } else {
      $ret .= $entry;
      $start_run = $end_run = $entry;

    }

  }

  if ( $start_run != $end_run ) {
    $ret .= q{-} . $end_run ;
  }
  $ret = q{[} . $ret . q{]};
  return $ret;
}

=head2 position_decode_string
=head2 tag_index_decode_string

these return the code string which can be put in a command line to decode the LSB_JOBINDEX back to position and tag_index
respectively, assuming that you used this code to generate them, including backticks

  my $sPostitionDecodeString = $oClass->position_decode_string();
  my $sTagIndexDecodeString  = $oClass->tag_index_decode_string();

=cut

sub position_decode_string {
  my ( $self ) = @_;

  return q{`echo $} . q{LSB_JOBINDEX/} . $LSF_INDEX_MULTIPLIER . q{ | bc`};

}

sub tag_index_decode_string {
  my ( $self ) = @_;

  return q{`echo $} . q{LSB_JOBINDEX%} . $LSF_INDEX_MULTIPLIER . q{ | bc`};

}

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item Readonly

=item Perl6::Slurp

=item npg_common::roles::log

=item npg_tracking::illumina::run::short_info

=item npg_tracking::illumina::run::folder

=item npg_tracking::glossary::tag

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Guoying Qi

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
