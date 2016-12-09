#############
# Created By: ajb
# Created On: 2009-10-08

package npg_common::fastqcheck;

use Moose;
use Moose::Util::TypeConstraints;
use Carp;
use Perl6::Slurp;
use Math::Round qw(round);
use English qw{-no_match_vars};
use Readonly;
use File::Basename;

use npg_qc::Schema;

our $VERSION = '0';
## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::fastqcheck

=head1 VERSION
i
=head1 SYNOPSIS

  use npg_common::fastqcheck;

  my $oFastqCheck = npg_common::fastqcheck->new({
    fastqcheck_path => $sFqCheck_filename,
  });

=head1 DESCRIPTION

This object provides some methods on top of a given fastqcheck file.

=head1 SUBROUTINES/METHODS

=cut

Readonly::Scalar our $EXPECTED_BINS_TOTAL         => 1000;

Readonly::Scalar our $POSITION_OF_QVAL_LOWEST     => 6;
Readonly::Scalar our $NUM_OTHER_FIELDS            => 6;
Readonly::Scalar our $TOP_ROW_NUM_FIELDS          => 9;
Readonly::Scalar our $MIN_NUM_WORDS               => 5;
Readonly::Scalar our $NUM_HEADERS                 => 4;
Readonly::Array  our @BASES                       => qw/  A    C    G    T    N /;

Readonly::Scalar our $THREE   => 3;
Readonly::Scalar our $FOUR    => 4;
Readonly::Scalar our $FIVE    => 5;
Readonly::Scalar our $SIX     => 6;
Readonly::Scalar our $SEVEN   => 7;
Readonly::Scalar our $EIGHT   => 8;

Readonly::Scalar our $RESULT_CLASS_NAME   => q[Fastqcheck];


subtype 'FastqcheckFile' => as 'Str'
                         => where { /[.]fastqcheck$/imsx }
                         => message {"$_ does not have .fastqcheck extension"};

subtype 'NonNegativeInteger'  => as 'Int'
                              => where { $_ >= 0; }
                              => message { "$_ must be a non-negative integer" };

subtype 'NonNegativeNum'  => as 'Num'
                          => where { $_ >= 0; }
                          => message { "$_ must be a non-negative number" };

=head2 fastqcheck_path

A path to a .fastqcheck file

=cut
has 'fastqcheck_path' => ( isa           => 'FastqcheckFile',
                           is            => 'ro',
                           required      => 0,
                          );

=head2 schema

DBIx schema object for the NPG QC database

=cut
has 'schema' =>    ( isa        => 'Object',
                     is         => 'ro',
                     required   => 0,
                     lazy_build => 1,
                   );
sub _build_schema {
    return npg_qc::Schema->connect();
}

=head2 db_lookup

A boolean flag indicating whether the files should be looked up in the db, defaults to false

=cut
has 'db_lookup' => ( isa        => 'Bool',
                     is         => 'ro',
                     required   => 0,
                     default    => 0,
                   );

=head2 file_content

Fastqcheck file as a string

=cut
has 'file_content' => ( isa        => 'Maybe[Str]',
                        is         => 'ro',
                        required   => 0,
                        lazy_build => 1,
                      );
sub _build_file_content {
    my $self = shift;

    my $text;
    if (!$self->db_lookup) {
        if ($self->fastqcheck_path) {
            $text = slurp($self->fastqcheck_path);
        }
    } else {
        my ($fname, $dir, $ext) = fileparse($self->fastqcheck_path);
        my $row = $self->schema->resultset($RESULT_CLASS_NAME)
            ->search({file_name => $fname,})->first;
        if ($row) {
            $text = $row->file_content;
	} else {
            croak $self->fastqcheck_path . ' is not in the long-term storage';
	}
    }
    return $text;
}


=head2 _lines_list

A list of lines found in the fastqcheck file. A private attribute.

=cut
has '_lines_list'     => ( isa           => 'ArrayRef',
                           is            => 'ro',
                           required      => 0,
                           lazy_build    => 1,
                         );
sub _build__lines_list {
    my $self = shift;
    my @all_lines = ();
    my $text = $self->file_content;
    if (defined $text) {
        @all_lines = split /\n/xms, $text;
    }
    return \@all_lines;
}


=head2 _is_empty

  Returns the number of total pf bases. A private attribute.
  my $zero_flag = $oFastqCheck->_is_empty();

=cut
has '_is_empty'  =>  (isa           => 'Bool',
                                is            => 'ro',
                                required      => 0,
                                lazy_build    => 1,
                               );
sub _build__is_empty {
    my $self = shift;

    my $result = 0;
    my @words =  split /\s+/xms, $self->_lines_list->[0];
    if (@words != $TOP_ROW_NUM_FIELDS) {
        if(@words >= $MIN_NUM_WORDS && $words[0] == 0 && $words[2] == 0)  {
            $result = 1;
        } else {
            croak $self->fastqcheck_path . q[ does not have expected structure];
	}
    }
    return $result;
}

=head2 total_pf_bases

  Returns the number of total pf bases. A private attribute with a public accessor.
  my $pf_bases = $oFastqCheck->total_pf_bases;

=cut
has '_total_pf_bases'   =>   (isa           => 'NonNegativeInteger',
                              is            => 'ro',
                              required      => 0,
                              lazy_build    => 1,
                              reader        => 'total_pf_bases',
			     );
sub _build__total_pf_bases {
    my $self = shift;
    my @words =  split /\s+/xms, $self->_lines_list->[0];
    return $words[2];
}


=head2 num_reads

  Returns the number of sequences found. A private attribute with a public accessor.
  my $iSequences = $oFastqCheck->num_reads;

=cut
has '_num_reads'    => (isa => 'NonNegativeInteger',
                            is => 'ro',
                            required   => 0,
                            reader     => 'num_reads',
                            lazy_build => 1,);

sub _build__num_reads {
    my $self = shift;

    if ($self->_is_empty) {
        return 0;
    }
    my @words =  split /\s+/xms, $self->_lines_list->[0];
    return $words[0];
}


=head2 average_read_length

  Returns the average cycle count. A private attribute with a public accessor.
  my $count = $oFastqCheck->average_read_length;

=cut
has '_average_read_length'    =>  (isa => 'NonNegativeNum',
                                  is => 'ro',
                                  required   => 0,
                                  reader     => 'average_read_length',
                                  lazy_build => 1,);

sub _build__average_read_length {
    my $self = shift;

    if ($self->_is_empty) {
        return 0;
    }
    my @words =  split /\s+/xms, $self->_lines_list->[0];
    if (@words < $SIX) {
        croak $self->fastqcheck_path . q[ does not have expected structure];
    }
    return $words[$FIVE];
}


=head2 max_read_length

  Returns the max cycle count. A private attribute with a public accessor.
  my $count = $oFastqCheck->max_read_length;

=cut
has '_max_read_length'    => (isa => 'NonNegativeInteger',
                              is => 'ro',
                              required   => 0,
                              reader     => 'max_read_length',
                              lazy_build => 1,
                             );
sub _build__max_read_length {
    my $self = shift;

    if ($self->_is_empty) {
        return 0;
    }
    my @words =  split /\s+/xms, $self->_lines_list->[0];
    if (@words < $EIGHT) {
        croak $self->fastqcheck_path . q[ does not have expected structure];
    }
    return $words[$SEVEN];
}


has '_non_negative_number'  => (isa        => 'NonNegativeNum',
                                is         => 'ro',
                                required   => 0,
                                writer     => '_check_number',
                                default    => 0,
			       );

has '_non_negative_int'  =>    (isa        => 'NonNegativeInteger',
                                is         => 'ro',
                                required   => 0,
                                writer     => '_check_int',
                                default    => 0,
			       );


=head2 _max_threshold

  Returns the max available threshold in the file. A private attribute.
  my $threshold = $oFastqCheck->_max_threshold;

=cut
has '_max_threshold'     =>    (isa        => 'Int',
                                is         => 'ro',
                                required   => 0,
                                lazy_build => 1,
			       );
sub _build__max_threshold {
    my $self = shift;
    if ($self->_is_empty) {
        return 0;
    }
    my @temp = split /\s+/xms, $self->_lines_list->[$NUM_HEADERS - 1];
    return (@temp - $NUM_OTHER_FIELDS - 2);
}


=head2 BUILD

  The last method called before a new object is returned by the constructor.
  Checks that the fastqcheck file is available and has expected structure.

=cut
sub BUILD {
    my $self = shift;

    my $path = $self->fastqcheck_path;
    if (!$path && !defined $self->file_content) {
      croak q[Either fastqcheck_path or file_content should be supplied];
    }

    my @lines = @{$self->_lines_list};
    if ( !@lines ) { croak qq[No data in file $path]; }

    if ( !$self->_is_empty ) {
        if (@lines < $NUM_HEADERS) {
             croak qq[$path does not have expected structure];
	}
        $self->total_pf_bases; # from first line

        if ($lines[$THREE] !~ /^Total/smx) {
            croak qq[$path does not have expected structure];
        }

        if (scalar @lines < $self->max_read_length + $NUM_HEADERS) {
            croak qq[Wrong number of lines in $path];
        }
    }
}


=head2 read_length

  Returns the numbers of cycles (the same as the max_read_length method). A convinience method.

=cut
sub read_length {
    my $self = shift;
    return $self->max_read_length;
}


=head2 qx_yield

  Returns a reference to an array with numbers of bases at and above qualities
  that are supplied in the argument array.

  Getting the number of bases at and above Q20 and Q30 for the whole sequence:
    my $qualities = [20, 30];
    my $qvals_list = $fqObject->qx_yield($qualities);
    my $num_bases_at_and_above_Q20 = $qvals_list->[0];
    my $num_bases_at_and_above_Q20 = $qvals_list->[1];

  Getting the number of bases at and above Q20 and Q30 for a particular cycle:
    my $qualities = [20, 30];
    my $cycle = 5;
    my $qvals_list = $fqObject->qx_yield($qualities, $cycle);
    my $num_bases_at_and_above_Q20 = $qvals_list->[0];
    my $num_bases_at_and_above_Q20 = $qvals_list->[1];

=cut
sub qx_yield {
    my ($self, $thresholds, $cycle_number) = @_;

    if (ref($thresholds) ne q[ARRAY]) {
        croak q[Wrong type of argument in ] .  __PACKAGE__ . q[::get_qvals];
    }

    my @results = ();

    if ($self->_is_empty) {
        foreach my $threshold (@{$thresholds}) {
            push @results, 0;
	}
    } else {
        my $zero_position = $POSITION_OF_QVAL_LOWEST;
        my @temp;
        if (!defined $cycle_number) {
            @temp = split /\s+/xms, $self->_lines_list->[$NUM_HEADERS - 1];
	} else {
            @temp = @{$self->_split_cycle_line($cycle_number)};
            $zero_position++;
	}
        my $end = @temp - 2;
        # if we are really unlucky and everything has been rounded one way
        ##no critic (ProhibitParensWithBuiltins ProhibitMagicNumbers)
        my $margin = int(($end - $zero_position) * 0.5 + 1);
        ##use critic

        my $total = 0;
        foreach my $position ($zero_position .. $end) {
            $self->_check_int($temp[$position]);
            $total += $temp[$position];
        }

        if (abs($EXPECTED_BINS_TOTAL - $total) > $margin) {
            croak qq[Bins total $total either too small or large, allowed difference from $EXPECTED_BINS_TOTAL is $margin];
	}

        foreach my $threshold (@{$thresholds}) {
            my $result = 0;
            if ($self->_check_int($threshold) && $threshold <= $self->_max_threshold) {
                my $start = $threshold + $zero_position;
                for my $i ($start .. $end) {
	            $result += $temp[$i];
                }
                $result = ($self->total_pf_bases * $result)/$total;
            }
            push @results, round($result);
        }
    }

    return \@results;
}


=head2 base_percentages

  Returns a reference to a hash where the keys are the bases (A, C, G, T, N) and
  the value for each key is the percent of the relevant base.

  Getting the total values:
    my $bases = $fqObject->base_percentages();
    foreach my $key (keys %{$bases}) {
        print $bases->{$key};
    }

  Getting the values for a particular cycle:
    my $cycle = 5;
    my $bases = $fqObject->base_percentages($cycle);

=cut
sub base_percentages {
    my ($self, $cycle_number) = @_;

    my $bh = {};
    if ($self->_is_empty) {
        foreach my $base (@BASES) {
            $bh->{$base} = 0;
	}
    } else {

        my $offset = 1;
        my @line;
        if (defined $cycle_number) {
            @line = @{$self->_split_cycle_line($cycle_number)};
            $offset++;
        } else {
            @line = split /\s+/xms, $self->_lines_list->[$NUM_HEADERS - 1];
        }

        foreach my $base (@BASES) {
            $bh->{$base} = $self->_check_number($line[$offset]);
            $offset++;
        }
    }
    return $bh;
}


=head2 gc_percentage

  Returns min and max gc percent

  Getting the total values:
    my ($min, $max) = $fqObject->gc_percentage();

  Getting the values for a particular cycle:
    my $cycle = 5;
    my ($min, $max) = $fqObject->qc_percentage($cycle);

=cut
sub gc_percentage {
    my($self, $cycle_number) = @_;

    my $bh = $self->base_percentages($cycle_number);
    my $gc = 0;
    if (exists $bh->{G}) {
        $gc += $bh->{G};
    }
    if (exists $bh->{C}) {
        $gc += $bh->{C};
    }

    my $max_gc = $gc;
    if (exists $bh->{N}) {
        $max_gc += $bh->{N};
    }

    return $gc, $max_gc;
}


sub _split_cycle_line {
    my ($self, $num) = @_;

    if (!defined $num || $num <= 0 || $num > $self->max_read_length) {
        croak qq[Invalid cycle number $num];
    }
    $self->_check_int($num);

    my @temp = split /\s+/xms, $self->_lines_list->[$NUM_HEADERS - 1 + $num ];

    if ($temp[1] != $num) {
        croak "Failed to locate line for cycle $num";
    }
    return \@temp;
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item Moose::Util::TypeConstraints

=item Math::Round

=item Perl6::Slurp

=item File::Basename

=item Carp

=item English -no_match_vars

=item Readonly

=item npg_qc::Schema

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Marina Gourtovaia E<lt>mg8@sanger.ac.ukE<gt> and Andy Brown E<lt>ajb@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) GRL 2010 by Andy Brown (ajb@sanger.ac.uk) and Marina Gourtovaia (mg8@sanger.ac.uk)

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
