#!/usr/bin/env perl
#########
# Author:        gq1
# Maintainer:    $Author: gq1 $
# Created:       2012-05-31
# Last Modified: $Date: 2010-10-27 15:18:11 +0100 (Wed, 27 Oct 2010) $
# Id:            $Id: irods_bam_loader.pl 11500 2010-10-27 14:18:11Z gq1 $
# $HeadURL: svn+ssh://svn.internal.sanger.ac.uk/repos/svn/new-pipeline-dev/useful_modules/branches/prerelease-45.0/scripts/irods_bam_loader.pl $
#

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::irods::BamDeletion;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision: 11500 $ =~ /(\d+)/smx; $r; };

npg_common::irods::BamDeletion->new_with_options()->process();

exit 0;

__END__

=head1 NAME

irods_bam_deletion.pl --bam_file irods_bam_file.bam --complete_deletion

=head1 VERSION

$LastChangedRevision: 11500 $

=head1 USAGE

irods_bam_deletion

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

=item npg_common::irods::BamDeletion

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2012 GRL, by Guoying Qi

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

