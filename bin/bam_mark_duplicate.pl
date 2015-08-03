#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );
use npg_common::sequence::BAM_MarkDuplicate;

our $VERSION = '0';

npg_common::sequence::BAM_MarkDuplicate->new_with_options()->process();

exit 0;

__END__

=head1 NAME

bam_mark_duplicate.pl

=head1 USAGE

 bam_mark_duplicate.pl --input_bam 4783_5.bam --output_bam 4783_6_mk.bam --metrics_json 4783_6_bam_flagstats.json

 --id_run 4783 --position 5 [ --tag_index 2 --human_split human]
 --sort_input
 --change_bam_header
 --replace_file 
 --bamcheck_flags

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

This script marks duplicates in a bam file using biobambam
sorting the input bam if requied.
It also generates bam flagstats json file.

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

=item lib

=item FindBin

=item npg_common::sequence::BAM::MarkDuplicate

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015 GRL

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

