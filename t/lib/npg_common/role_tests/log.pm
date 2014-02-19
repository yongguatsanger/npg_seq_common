#############
# $Id$
# Created By: ajb
# Mast Maintained By: $Author$
# Created On: 2009-10-01
# Last Changed On: $Date$
# $HeadURL$

package npg_common::role_tests::log;
use Moose;
use Carp;
use English qw{-no_match_vars};
use Readonly;

with qw{npg_common::roles::log};

Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mxs; $r; };



no Moose;
__PACKAGE__->meta->make_immutable;
1;
__END__

=head1 NAME

npg_common::role_tests::log

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item Carp

=item English -no_match_vars

=item Readonly

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
