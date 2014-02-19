#############
# $Id$
# Created By: ajb
# Mast Maintained By: $Author$
# Created On: 2009-10-01
# Last Changed On: $Date$
# $HeadURL$

package npg_common::roles::run::fs_resource;
use Moose::Role;

use strict;
use warnings;
use File::Spec::Functions qw(splitdir catfile catdir);
use Carp qw(carp cluck croak confess);
use Sys::Filesystem::MountPoint qw(:all);

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

requires q{runfolder_path};

has fs_resource => (
  is => 'ro',
  isa => 'Str',
  lazy_build => 1,
  documentation => 'filesystem resource string, can be worked out, so only set if you need to override this',
);

sub _build_fs_resource {
  my ($self) = @_;

  if ($ENV{TEST_FS_RESOURCE}) { return $ENV{TEST_FS_RESOURCE}; }

  my $r = join '_', grep {$_} splitdir path_to_mount_point($self->runfolder_path());
  return join q(),$r=~/\w+/xsmg;
}

1;
__END__

=head1 NAME

npg_common::roles::run::fs_resource

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  package MyPackage;
  use Moose;
  with qw{npg_tracking::illumina::run::short_info
          npg_tracking::illumina::run::folder
          npg_common::roles::run::fs_resource};

  my $oPackage = MyPackage->new({
    runfolder_path => q{/string/to/a/run_folder},
  });

  my $oPackage = MyPackage->new({
    run_folder     => q{run_folder},
  });

=head1 DESCRIPTION

This package needs to have something provide a runfolder_path, either declared in your class,
or because you have previously declared with npg_tracking::illumina::run::folder, which is the preferred
option.

=head1 SUBROUTINES/METHODS

=head2 fs_resource - returns the fs_resource for the given runfolder_path

  my $sFSResource = $oPackage->fs_resource();

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item Readonly

=item Cwd

=item File::Spec::Functions

=item Sys::Filesystem::MountPoint

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2009 Andy Brown (ajb@sanger.ac.uk)

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
