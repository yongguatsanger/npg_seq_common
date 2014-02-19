#!/usr/bin/env perl
#########
# Author:        Jennifer Liddle (js10)
# Maintainer:    $Author$
# Created:       01 March 2012
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

#########################
# This script copies the Interop files into iRODS
##########################

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::irods::run::Interop;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/smx; $r; };

npg_common::irods::run::Interop->new_with_options()->process();

exit 0;

__END__

=head1 NAME

irods_interop_loader.pl

=head1 VERSION

$LastChangedRevision$

=head1 USAGE

irods_interop_loader.pl --id_run 1234

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

=head1 EXIT STATUS

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item FindBin

=item npg_common::irods::run::Interop

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Jennifer Liddle (js10@sanger.ac.uk)

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2012 GRL, by Jennifer Liddle

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

