#########
# Author:        Jennifer Liddle (js10)
# Maintainer:    $Author$
# Created:       2012 02 07
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::run::Interop;

use strict;
use warnings;
use Moose;
use Carp;
use Try::Tiny;
use File::Spec;

use npg_common::irods::Loader;

with qw{MooseX::Getopt npg_common::roles::log};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $DEFAULT_ROOT_DIR   => q{/seq/};
Readonly::Array our @FILE_LIST => ('RunInfo.xml', 'runParameters.xml', 'InterOp');

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::run::Interop

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 SUBROUTINES/METHODS

=head2 irods_root

Root directory for irods

=cut
has 'irods_root' => (isa => 'Str',
                     is  => 'rw',
                     default => $DEFAULT_ROOT_DIR,
                    );
=head2 runfolder_path

Directory to look for files

=cut
has 'runfolder_path' => (isa => 'Str',
                         is  => 'rw',
                         required   => 1,
                        );

=head2 id_run

The run ID

=cut
has 'id_run' => (isa => 'Int',
                 is  => 'rw',
                 required   => 1,
                );

=head2 process

=cut

sub process {
  my $self = shift;

  my $id_run = $self->id_run();

  $self->log("Loading Interop files for run $id_run");

  foreach my $file (@FILE_LIST) {

    $self->log("---- Loading file $file");

    my $meta  = {
                  id_run => $id_run,
                };

    my $loader = npg_common::irods::Loader->new({
       file        => File::Spec->rel2abs($file,$self->runfolder_path()),
       collection  => $self->irods_root() . $id_run,
       meta_data   => $meta,
    });

    try {
      $loader->run();
    } catch {
      if ( /No such file/smi ) {;
        $self->log("Warning: $_");
      } else {
        croak $_;
      }
    }
  }

  return 1;
}

no Moose;

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item MooseX::Getopt

=item Carp

=item Readonly

=item Try::Tiny

=item File::Spec

=item npg_common::irods::Loader

=item npg_common::roles::log

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Guoying Qi

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
