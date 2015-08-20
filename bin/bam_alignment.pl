#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::sequence::BAM_Alignment;

our $VERSION = '0';

npg_common::sequence::BAM_Alignment->new_with_options()->generate();

exit 0;

__END__

=head1 NAME

bam_alignment.pl

=head1 VERSION

=head1 USAGE

bam_alignment.pl --input  input.bam --is_paired_read --output_prefix output --reference --human_reference --phix_reference --ref_dict picard_reference_dict [--no_alignment]

or

bam_alignment.pl --input  input.bam --output_prefix output --id_run 6551 --position 1 --tag_index

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

This script takes a bam file, possible with phix alignment
doing alignment against the given reference
and split phix and non-consented human reads if required.

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

=item npg_common::sequence::BAM_Alignment

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

