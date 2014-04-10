#############
# Created By: Marina Gourtovaia
# Created On: 6 September 2010

package npg_common::roles::run::lane::file_names;

use Moose::Role;
use Carp;
use npg_tracking::glossary::tag;

use Readonly;
our $VERSION = '0';

sub file_ext {
    my ($self, $file_extension) = @_;
    return $file_extension ? q[.].$file_extension : q[];
}

sub generate_filename {
    my ($self, $file_extension, $end) = @_;
    return [$self->create_filename($file_extension, $end)];
}

sub create_filename {
    my ($self, $file_extension, $end) = @_;

    my $tag_label = $self->tag_label;
    my $type_string = $self->can(q[sequence_type]) && $self->sequence_type ? q[_] . $self->sequence_type : q[];
    my $extension =  $self->file_ext($file_extension);

    my $end_string = q[];
    if ($end) {
        if ($end eq q[1] || $end eq q[2] || $end eq q[t]) {
	    $end_string = q[_] . $end;
        } else {
            croak qq[Unrecognised end string $end];
	}
    }
    return $self->id_run . q[_] . $self->position . $end_string . $tag_label . $type_string . $extension;
}

sub nonsplit2split {
    my ($self, $fname, $non_human, $extension) = @_;

    if (!$fname) {
        croak q[File name argument is not defined in nonsplit2split];
    }

    if (!defined $non_human) {
        return $fname;
    }
    $extension = $self->file_ext($extension);
    if ($extension) {
        $fname = substr $fname, 0, -(length $extension);
    }

    my $tag_delim = $npg_tracking::glossary::tag::TAG_DELIM;
    my ($base, $tag_index) = $fname =~ /(.+)\Q$tag_delim\E(\d+)$/smx;
    if (!$base) {
        $base = $fname;
        $tag_index = undef;
    }

    my $human_split_string = $non_human ? q[nonhuman] : q[human];
    $base .= q[_] . $human_split_string;
    if (defined $tag_index) {
        $base .= $tag_delim . $tag_index;
    }
    return $base . $extension;
}

no Moose::Role;

1;

__END__

=head1 NAME

npg_common::roles::run::lane::file_names

=head1 VERSION


=head1 SYNOPSIS

=head1 DESCRIPTION

A role for generating file names

=head1 SUBROUTINES/METHODS

=head2 file_ext - A helper function to generate a file extension.

=head2  generate_filename - Given an end value and the file extension generates a file name for the source
following convention id-run_position_end.extension. The default default end is 1. Possible file ends are
1, 2 and t. Returns a reference to an array of (one) filename.

  $obj->generate_filename(q[fastq]);
  my $file_end = 2;
  $obj->generate_filename(q[fastq], $file_end);

=head2  create_filename - Given an end value and the file extension generates a file name for the source
following convention id-run_position_end.extension. The default default end is 1. Possible file ends are
1, 2 and t. Returns the filename

  $obj->generate_filename(q[fastq]);
  my $file_end = 2;
  $obj->generate_filename(q[fastq], $file_end);

=head2 nonsplit2split - Given a file name and its extension, generates a name for one of files coming out of a human/non-human split

  my $non_human;
  print $obj->nonsplit2split('222_1_1.fastq', $non_human, 'fastq')   #prints 222_1_1.fastq

  $non_human = 1;
  print $obj->nonsplit2split('222_1_1.fastq', $non_human, 'fastq')   #prints 222_1_1_nonhuman.fastq
  print $obj->nonsplit2split('222_1_1#3.fastq', $non_human, 'fastq') #prints 222_1_1_nonhuman#3.fastq
  print $obj->nonsplit2split('222_1_1#3.fastq', $non_human) #prints 222_1_1_nonhuman#3

  $non_human = 0;
  print $obj->nonsplit2split('222_1_1.fastq', $non_human, 'fastq')   #prints 222_1_1_human.fastq
  print $obj->nonsplit2split('222_1_1#3.fastq', $non_human, 'fastq') #prints 222_1_1_human#3.fastq
  print $obj->nonsplit2split('222_1_1#3.fastq', $non_human) #prints 222_1_1_human#3

=head2

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item npg_tracking::glossary::tag

=item Carp

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Marina Gourtovaia

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
