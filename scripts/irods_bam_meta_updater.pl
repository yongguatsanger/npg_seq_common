#!/usr/bin/env perl
#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-07-12
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::irods::BamMetaUpdater;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/smx; $r; };

npg_common::irods::BamMetaUpdater->new_with_options()->process();

exit 0;

__END__

=head1 NAME

irods_bam_meta_updater.pl

=head1 VERSION

$LastChangedRevision$

=head1 USAGE

irods_bam_meta_updater.pl                 
	--verbose
	--id_run                 the id_run filter                
	--min_id_run             the minimum id_run filter
	--max_id_run             the maximum id_run filter
	--lane                   the lane number filter
	--tag_index              the plex tag index number filter

	--output                 The file name to store imeta sub commands, which will be printed out if not given

	--not_dry_run            If true, imeta sub commands will be run and irods meta data will be updated
	
	--not_check_warehouse     do not check warehouse for lims data change
	
   --check_bam               check bam header for meta data change
   
   --check_md5               check md5 value in irods meta data match the value reported by ils command

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

=item npg_common::irods::BamMetaUpdater

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2011 GRL, by Guoying Qi

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

