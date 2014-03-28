#########
# Author:        Marina Gourtovaia
# Maintainer:    $Author$
# Created:       05 May 2009
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::diagram::visio_histo_google;

use strict;
use warnings;
use Readonly;
use Carp;
use English qw(-no_match_vars);
use POSIX;

## no critic(ProhibitParensWithBuiltins Capitalization ProhibitMixedCaseSubs ProhibitManyArgs)

our $VERSION = '0';

Readonly::Scalar our $DELIMETER      => q[&];
Readonly::Scalar our $GOOGLE_URL     => q[http://chart.apis.google.com/chart?];
Readonly::Scalar our $GOOGLE_MAX_ENCODED_VALUE => 4095;
Readonly::Scalar our $DEFAULT_Y_MAX  => 100;
Readonly::Scalar our $BASE_ENCODE    => 64;
Readonly::Scalar our $YMAX_TOLERANCE => 0.01;
Readonly::Array  our @GOOGLE_ENCODE  => ( q{A}..q{Z}, q{a}..q{z}, 0..9, q{-}, q{.} );

sub new {
    my $class = shift;
    my $self = shift || {};
    bless $self, $class;

    $self->{encode}                 ||= 0;
    $self->{y_max}                  ||= $DEFAULT_Y_MAX;
    $self->{cht_chart_type}         ||= ['cht' , 'bvg'];
    $self->{chco_chart_colour}      ||= ['chco', '4D89F9'];
    $self->{chtt_chart_title}       ||= ['chtt', 'batch+X,run+Y,+lane+Z'];
    $self->{chs_chart_size}         ||= ['chs' , '300x250'];
    $self->{chd_chart_data}         ||= ['chd' , 't:20,30,40,23'];
    $self->{chds_chart_min_max}     ||= ['chds', '0,40'];
    $self->{chxt_chart_axises}      ||= ['chxt', 'x,y'];
    $self->{chxr_chart_axis_labels} ||= ['chxr', '0,1,4|1,0,40,10'];
    $self->{chbh_bar_props}         ||= ['chbh', '5,1,1']; #absolute units, 5 px bar, 1px distance between bars

    return $self;
}

sub chdl_chart_legend {
  my ( $self, $chart_legend ) = @_;

  if ( $chart_legend ) {
    if ( ! ref $chart_legend ) {
      $self->{chart_legend} = [ 'chdl', $chart_legend ];
    } elsif ( ref $chart_legend eq q{ARRAY} ) {
      $self->{chart_legend} = [ 'chdl', (join q{|}, @{ $chart_legend } ) ];
    } else {
      $self->{chart_legend} = [ 'chdl', (join q{|}, sort keys %{ $chart_legend } ) ];
    }
  }

  return $self->{chart_legend};
}

sub encode {
  my ( $self, $encode ) = @_;
  if ( $encode ) {
    $self->{encode} = $encode;
  }
  return $self->{encode};
}

sub y_max {
  my ( $self, $y_max ) = @_;
  if ( $y_max ) {
    $self->{y_max} = $y_max;
  }
  return $self->{y_max};
}

sub get_diagram_string {

    my ( $self, $no_data ) = @_;
    my @params = ();
    foreach my $key (sort keys %{$self}) {
        next if $key eq q{encode};
        next if $key eq q{y_max};
        next if ( $no_data && $key =~ /chtt|chs|chd|chds|chxt|chxr|chbh/xms );

        if ($key ne q[data] && $key ne q[xmax] && $key ne q[xmin]) {
            my $pair = $self->{$key};
            if ($pair && @{$pair} > 0) {
                push @params, $pair->[0] . q[=] . $pair->[1];
            }
        }
    }
    if ( $no_data ) {
      push @params, 'chs=70x250';
    }
    return $GOOGLE_URL . join("$DELIMETER", @params);
}

sub _encode_value {
  my ( $self, $value ) = @_;

  my $value_is_greater_than_0 = $value > 0 ? 1 : 0;
  my $y_max = $self->y_max();

  if ( $value >= $y_max ) {
    if ($value - $y_max > $YMAX_TOLERANCE) {
      croak qq[Value to encode $value is much bigger than max y $y_max (tolerance $YMAX_TOLERANCE)];
    }
    $value = $GOOGLE_MAX_ENCODED_VALUE;
  } else {
    $value = floor ( $value * $GOOGLE_MAX_ENCODED_VALUE / $y_max );
  }

  # we want to see a small bar if there is something.
  # If the y_max is too big, the above calculation could some values disappear
  if ( ! $value && $value_is_greater_than_0 ) {
    $value = 1;
  }

  my $char1 = @GOOGLE_ENCODE[ ( floor $value/$BASE_ENCODE ) ];
  my $char2 = @GOOGLE_ENCODE[ ( $value % $BASE_ENCODE ) ];

  return $char1 . $char2;
}

sub set_data { ## no critic (Subroutines::ProhibitExcessComplexity)

    my ( $self, $data ) = @_;
    my $raw_data;
    if ( ref $data eq q{ARRAY} ) {
      @{ $raw_data } = @{ $data };
    } else {
      %{ $raw_data } = %{ $data };
    }
    $self->{data} = $raw_data;

    my $prepared_data;
    my $data_set_join_char    = $self->encode() ? q{,}
                              :                   q{|}
                              ;
    my $data_point_join_char  = $self->encode() ? q{}
                              :                   q{,}
                              ;
    if ( ref $data eq q{ARRAY} ) {
      if ( ! ref $data->[0] ) {
        if ( $self->encode() ) {
          @{ $data } = map { $self->_encode_value( $_ ) } @{ $data };
        }
        $prepared_data = join( $data_point_join_char, @{ $data } );
      } else {
        my @sets;
        foreach my $data_set ( @{ $data } ) {
          if ( scalar @{ $data_set } ) {
            if ( $self->encode() ) {
              @{ $data_set } = map { $self->_encode_value( $_ ) } @{ $data_set };
            }
            push @sets, ( join( $data_point_join_char, @{ $data_set } ) );
          } else {
            if ( $self->encode() ) {
              push @sets, q{AA};
            } else {
              push @sets, 0;
            }
          }
        }
        $prepared_data = join( $data_set_join_char, @sets );
      }
    } else {
      my @sets;
      foreach my $key ( sort keys %{ $data } ) {
        if ( scalar @{ $data->{$key} } ) {
          if ( $self->encode() ) {
            @{ $data->{$key} } = map { $self->_encode_value( $_ ) } @{ $data->{$key} };
          }
          push @sets, ( join( $data_point_join_char, @{ $data->{$key} } ) );
        } else {
          if ( $self->encode() ) {
            push @sets, q{AA};
          } else {
            push @sets, 0;
          }
        }
      }
      $prepared_data = join( $data_set_join_char, @sets );
    }

    my $encoding_char = $self->encode() ? q{e}
                      :                   q{t}
                      ;

    $self->{chd_chart_data} = ['chd' , $encoding_char . q{:} . $prepared_data];

    # we'll reset this just in case we have encoded,
    # but other code doesn't expect $data to be changed
    if ( ref $data eq q{ARRAY} ) {
      @{ $data } = @{ $raw_data };
    } else {
      %{ $data } = %{ $raw_data };
    }

    return 1;
}


sub set_axisY_min_max {

    my ($self, $min, $max) = @_;
    $self->{chds_chart_min_max} = ['chds', $min . q[,] . $max];
    return 1;
}


sub set_chart_title {

    my ($self, $title) = @_;
    $self->{chtt_chart_title} = ['chtt', join(q[+], split(q[ ],  $title))];
    return 1;
}


sub set_chart_labels {

    my ($self, $xmin, $xmax, $xdelta, $ymin, $ymax, $ydelta) = @_;

    $self->{xmax} = $xmax;
    $self->{xmin} = $xmin;
    my @x = (0, $xmin, $xmax, $xdelta);
    my @y = (1, $ymin, $ymax, $ydelta);
    my $labels = join(q[,], @x) . q[|] . join(q[,], @y);
    $self->{chxr_chart_axis_labels} = ['chxr', $labels];
    return 1;
}


sub set_chart_size {

    my ($self, $width, $height) = @_;
    $self->{chs_chart_size} = ['chs' , $width . 'x' . $height];
    return 1;
}


sub set_bar_size {

    my ($self, $bar_width, $bar_distance) = @_;
    if (!defined $bar_distance || $bar_distance < 0) {
        $bar_distance = 0;
    }
    if (!$bar_width) {
        carp q[Bar width has to be set for a google chart];
    }
    if ($bar_width <= 0) {
        carp qq[Bar width for a google chart should be positive; value given $bar_width];
    }
    $self->{chbh_bar_props}         = ['chbh', $bar_width . ',1,' . $bar_distance];
    return 1;
}


sub set_vert_band_markers {

    my ($self, $bands) = @_;
    if (!exists $self->{data} || !$self->{data}) {
        croak q[Please set data array first];
    }
    if (!exists $self->{xmax} || !exists $self->{xmin}) {
        croak q[Please set axis labels first];
    }

    my $num_bands = @{$bands};
    if ($num_bands % 2 != 0) {
        croak q[The length of the arg array should be an even number];
    }

    ## no critic (ProhibitBooleanGrep)
    if ((grep {$_<$self->{xmin} || $_>$self->{xmax}} @{$bands})> 0) {
        croak q[Band values should be within the X axis min and max values];
    }

    my $count = 0;
    my $labels = q[];
    my $range = $self->{xmax} - $self->{xmin};
    my $colour = q[80C65A];
    while ($count < $num_bands-1) {
        if ($count == 0) {$colour = q[FF0000];}
        my $lower = ($bands->[$count]-$self->{xmin})/$range;
        $count++;
        my $upper = ($bands->[$count]-$self->{xmin})/$range;
        if ($count != 1) {$labels = $labels . q[|];}
        $labels = $labels . q[R,] . $colour . q[,0,] . $lower . q[,] . $upper;
    }
    $self->{chm_markers} = ['chm', $labels];
    return 1;
}




1;
__END__

=head1 NAME

npg_common::diagram::visio_histo_google - an API for producing bar diagrams with Google Chart API

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  use npg_common::diagram::visio_histo_google;

=head1 DESCRIPTION

npg_common::diagram::visio_histo_google

=head1 USAGE

  my $dia = npg_common::diagram::visio_histo_google->new();
  my $array = [1,2,2,3,5,5,7];
  $dia->set_title(q[My Chart]);
  $dia->set_data($array);
  $dia->set_axisY_min_max(0,7);
  my $url = $dia->get_diagram_string();

=head1 OPTIONS

=head1 SUBROUTINES/METHODS

=head2 new returns a npg_common::diagram::visio_histo_google object
  $dia = npg_common::diagram::visio_histo_google->new();

=head2 get_diagram_string returns a URL for generating a Google chart

=head2 set_chart_title sets the title of the chart

=head2 set_data  sets data

=head2 set_axisY_min_max sets axis Y min and max values

=head2 set_chart_labels - setd label for axises

=head2 set_chart_size - sets chart width and height

=head2 set_bar_size - sets bar width and a distance between bars

=head2 set_vert_band_markers

=head2 encode

Accessor to instruct the url generation to encode the dataset to provide url compression related to the base64 encoding in the google api.

Can optionally be provided on construction. Defaults to 0 (off).

=head2 y_max

Accessor to set y_max to be used, primarily for encoding.

Can optionally be provided on construction. Defaults to 100.

=head2 chdl_chart_legend

store the chart legend in the object. Will convert on storage to an arrayref ready for inclusion in the google uri.

Returns this arrayref. Due to the conversion, you cannot add this on object construction.

Can take a string, which it will treat as the string for the uri, an arrayref, which it will use the | char to join together to make the string for the uri, or a hashref, in which case it will sort the keys, and use these joined by the | for the string

To avoid possible problems, only use the hashref if your data has come from this, and you want it sorted by key

  my $aGoogleChartLegend = $oDia->chdl_chart_legend( $sLegend );
  my $aGoogleChartLegend = $oDia->chdl_chart_legend( $aLegend );
  my $aGoogleChartLegend = $oDia->chdl_chart_legend( $hData );

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item Carp

=item English -no_match_vars

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

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

=cut

