#############
# Created By: Marina Gourtovaia
# Created On: 23 April 2010

package npg_common::run::file_finder;

use Moose;
use Carp;
use English qw{-no_match_vars};
use File::Spec::Functions qw(catfile);
use File::Basename;
use Readonly;

use npg_qc::Schema;

with    qw/ npg_tracking::glossary::lane
            npg_tracking::glossary::tag
            npg_common::roles::run::lane::file_names
            npg_tracking::illumina::run::short_info
            npg_tracking::illumina::run::folder
          /;

our $VERSION = '0';

Readonly::Scalar our $FILE_EXTENSION      => q[fastq];
Readonly::Scalar our $RESULT_CLASS_NAME   => q[Fastqcheck];


has 'db_lookup' =>  (
                             isa            => 'Bool',
                             is             => 'ro',
                             required       => 0,
                             writer         => '_set_db_lookup',
                             default        => 1,
			    );

has 'file_extension'   =>   (
                             isa            => 'Maybe[Str]',
                             is             => 'rw',
                             required       => 0,
                             default        => $FILE_EXTENSION,
		            );

has 'with_t_file'      =>   (isa            => 'Bool',
                             is             => 'ro',
                             required       => 0,
                             default        => 0,
		            );

has 'lane_archive_lookup' =>   (isa            => 'Bool',
                                is             => 'ro',
                                required       => 0,
                                default        => 1,
		               );

has 'qc_schema' =>    ( isa        => 'npg_qc::Schema',
                        is         => 'ro',
                        required   => 0,
                        lazy_build => 1,
                      );
sub _build_qc_schema {
    my $self = shift;
    my $schema = npg_qc::Schema->connect();
    return $schema;
}

has 'globbed'         =>   ( isa            => 'HashRef',
                             is             => 'ro',
                             required       => 0,
                             lazy_build     => 1,
		           );
sub _build_globbed {
    my $self = shift;

    my $hfiles = {};
    if ( $self->file_extension eq q[fastqcheck] && $self->db_lookup) {
        my $pattern = join q[_], $self->id_run, $self->position;
        $pattern .= q[%];
        if ($self->file_extension) {
            $pattern .= $self->file_extension;
	}
        my @rows = $self->qc_schema->resultset($RESULT_CLASS_NAME)->search(
	              {file_name => {'like' => $pattern,},},
	              {columns => 'file_name',},
                   )->all;
        foreach my $row (@rows) {
            my $fname = $row->file_name;
            $hfiles->{$fname} = $fname;
	}
    }

    if ((scalar keys %{$hfiles}) == 0) {
        my $path = (defined $self->tag_index && $self->lane_archive_lookup) ?
    File::Spec->catfile($self->archive_path, $self->lane_archive) : $self->archive_path;
        my $glob = catfile($path, q[*]);
        if ($self->file_extension) {
            $glob .= $self->file_extension;
	}
        my @files = glob "$glob";
        foreach my $file (@files) {
            my ($fname, $dir, $ext) = fileparse($file);
            $hfiles->{$fname} = $file;
        }
        $self->_set_db_lookup(0);
    }
    return $hfiles;
}


sub BUILD {
    my $self = shift;
    $self->_test_options_compatibility();
    return;
}


sub _test_options_compatibility {
    my $self = shift;
    if (defined $self->tag_index && $self->with_t_file) {
        croak 'tag_index and with_t_file attributes cannot be both set';
    }
}


sub files {
    my $self = shift;

    my $fnames = {};
    my $f = $self->create_filename($self->file_extension);
    if (exists $self->globbed->{$f}) {
        $fnames->{forward} =  $self->globbed->{$f};
    }

    if (!exists $fnames->{forward}) {
        my $forward  =  $self->create_filename($self->file_extension, 1);
        if (exists $self->globbed->{$forward}) {
            $fnames->{forward} =  $self->globbed->{$forward};
        }
        my $reverse =  $self->create_filename($self->file_extension, 2);
        if (exists $self->globbed->{$reverse}) {
	    $fnames->{reverse} =  $self->globbed->{$reverse};
        }
    }

    if ($self->with_t_file) {
        my $tag = $self->create_filename($self->file_extension, q[t]);
        if (exists $self->globbed->{$tag}) {
             $fnames->{tags} = $self->globbed->{$tag};
        }
    }
    return $fnames;
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 NAME

npg_common::run::file_finder

=head1 VERSION

=head1 SYNOPSIS

  my $finder = npg_common::run::file_finder->new(id_run => 2222, position => 1);

=head1 DESCRIPTION

Locates the files for a lane.

=head1 SUBROUTINES/METHODS

=head2 file_extension - an attribute, defaults to fastq, can be empty

=head2 with_t_file - a boolean attribute, defaults to false, determines whether a file for tags will be looked up.

=head2 qc_schema - DBIx schema object for the NPG QC database

=head2 db_lookup - a boolean attribute defining whether a lookup in teh qc db should be performed.
Is reset by the files method to show whether the file names do come from the db lookup. The
default initial value is true.

=head2 lane_archive_lookup - a boolean attribute indicating whether the files for tags (plexes) are
expected to be in the lane archive under the archive folder; defaults to true;

=head2 globbed - a lazily buils hash ref containing all actually available file names for a lane

=head2 files - Returns a reference to a hash containing paths to\names of a forward and reverse(if any) input files and
also a file for a tag if requested. The possible keys are forward, reverse, and tags.

  my $hash = $finder->files;
  my $forward = $hash->{forward};

=head2 BUILD - some sanity checking before returning a reference to a newly created object to the caller; for internal use only

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item Carp

=item English -no_match_vars

=item Readonly

=item File::Spec::Functions -catfile

=item File::Basename

=item npg_qc::Schema

=item npg_common::roles::run::lane::file_names

=item npg_tracking::illumina::run::short_info

=item npg_tracking::illumina::run::folder

=item npg_tracking::glossary::lane

=item npg_tracking::glossary::tag

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Marina Gourtovaia

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
