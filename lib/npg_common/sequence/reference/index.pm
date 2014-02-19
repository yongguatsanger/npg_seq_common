
#########
# Author:        Marina Gourtovaia
# Maintainer:    $Author: gq1 $
# Created:       12 July 2011
# Last Modified: $Date: 2011-06-30 16:05:58 +0100 (Thu, 30 Jun 2011) $
# Id:            $Id: find.pm 13501 2011-06-30 15:05:58Z gq1 $
# $HeadURL: svn+ssh://svn.internal.sanger.ac.uk/repos/svn/new-pipeline-dev/useful_modules/branches/prerelease-37.0/lib/npg_common/sequence/reference/roles/find.pm $
#

package npg_common::sequence::reference::index;

use strict;
use warnings;
use Carp;
use File::Path qw(make_path);
use Cwd qw(cwd);
use File::Spec::Functions;
use Moose;

with 'npg_tracking::data::reference::list';

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision: 13501 $ =~ /(\d+)/smx; $r; };
## no critic (Documentation::RequirePodAtEnd ProhibitBacktickOperators)

=head1 NAME

npg_common::sequence::reference::index

=head1 VERSION

$Revision: 13501 $

=head1 SYNOPSIS

 my $i = npg_common::sequence::reference::index->new(all_species => 0,
                                                     aligner => q[smalt],
                                                     aligner_options => q[index -k 13 -s 4],
                                                     target_ref_repository => q[mypath],
                                                     create_tree => 1,
                                                    );
 $i->generate;

=head1 DESCRIPTION

A class for generating aligner-specific binary indices for a collection of references.

=head1 SUBROUTINES/METHODS

=cut


=head2 defaults

Only indices for default strains will be generated

=cut
has 'defaults' =>   (isa       =>'Bool',
                     is        => 'ro',
                     required  => 0,
                     default   => 0,
		    );

=head2 create_tree

If this option is set, a ref repository directory tree will be created

=cut
has 'create_tree' =>   (isa       =>'Bool',
                        is        => 'ro',
                        required  => 0,
                        default   => 0,
		       );

=head2 aligner

aligner name

=cut
has 'aligner' =>      (isa       => 'Str',
                       is        => 'ro',
                       required  => 0,
                       default   => q[smalt],
		      );

=head2 aligner_options

aligner options

=cut
has 'aligner_options' =>      (isa       => 'Str',
                               is        => 'ro',
                               required  => 0,
                               default   => q[index],
       	                      );

=head2 target_ref_repository

path to the reference repository where the index files have to be created

=cut
has 'target_ref_repository' =>      (isa       => 'Str',
                                     is        => 'ro',
                                     required  => 1,
		                    );

sub _split_key {
  my ($self, $key) = @_;
  return split /:/smx, $key;
}

sub _apath {
  my ($self, $species, $strain, $aligner) = @_;
  my $path = catfile($species, $strain, q[all]);
  if ($aligner) {
    $path = catfile($path, $aligner);
  }
  return $path;
}

=head2 generate

Generates aligner-specific binary indices for a collection of references

=cut
sub generate {
  my $self = shift;

  if($self->create_tree) {
    $self->create_ref_repository();
  }

  foreach my $pair (sort keys %{$self->repository_contents}) {
    if ($self->defaults && !$self->repository_contents->{$pair}->{default}) { next; }
    my ($species, $strain) = $self->_split_key($pair);
    my $tools_dir = catfile($self->ref_repository, $self->_apath($species, $strain));
    my $ref_root = $self->ref_file_prefix($tools_dir);
    my $target = catfile($self->target_ref_repository, $self->_apath($species, $strain, $self->aligner), $ref_root);
    my $source = catfile($tools_dir, q[fasta], $ref_root);
    my $command = join q[ ], $self->aligner, $self->aligner_options, $target, $source;
    carp $command;
    `$command`;
  }
  return;
}

=head2 create_ref_repository

Creates a directory tree for the target reference repository

=cut
sub create_ref_repository {
  my $self = shift;

  my $path = $self->target_ref_repository;
  if (!-d $path) {
    croak qq[$path directory does not exist];
  }
  if (!-w $path) {
    croak qq[Cannot write to $path];
  }

  foreach my $pair (sort keys %{$self->repository_contents}) {
    my $strain_is_default = $self->repository_contents->{$pair}->{default} ? 1 : 0;
    if ($self->defaults && !$strain_is_default) { next; }
    my ($species, $strain) = $self->_split_key($pair);
    my $apath = catfile($path, $self->_apath($species, $strain, $self->aligner));
    make_path($apath);

    my $fasta_dir = $self->_apath($species, $strain, q[fasta]);
    my $target_fpath = catfile($path, $fasta_dir);
    my $fpath =  catfile($self->ref_repository, $fasta_dir);
    `ln -s $fpath $target_fpath`;

    if ($strain_is_default) {
      my $spath = catfile($path, $species, $strain);
      my $default = catfile($path, $species, q[default]);
      `ln -s $spath $default`;
    }
  }
  return;
}

no Moose;
1;

__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item Moose

=item Carp

=item File::Spec::Functions

=item Readonly

=item File::Path

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Author: Marina Gourtovaia E<lt>mg8@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2011 GRL, by Marina Gourtovaia

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

