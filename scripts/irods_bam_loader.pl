#!/usr/bin/env perl
#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       18 January 2010
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::irods::run::Bam;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/smx; $r; };

npg_common::irods::run::Bam->new_with_options()->process();

exit 0;

__END__

=head1 NAME

irods_bam_loader.pl

=head1 VERSION

$LastChangedRevision$

=head1 USAGE

irods_bam_loader.pl --id_run 5169

or

irods_bam_loader.pl --archive_path archive_path

or

irods_bam_loader.pl --id_run 5169 --positions 1 [--exclude_bam_index]

to exclude bam index files in the loading, use --exclude_bam_index flag

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

=item npg_common::irods::run::Bam

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

