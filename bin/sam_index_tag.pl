#!/usr/bin/env perl
#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       09 Aug 2010
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );

use npg_common::sequence::SAM_Index_Tag;

use Readonly;
our $VERSION = '0';

npg_common::sequence::SAM_Index_Tag->new_with_options()->process();

exit 0;

__END__

=head1 NAME

sam_index_tag.pl

=head1 VERSION


=head1 USAGE
  
  samtools view -h 5008_1#2.bam | scripts/sam_index_tag.pl --index_fastq_file 5008_1_t.fastq
  
or scripts/sam_index_tag.pl --index_fastq_file 5008_1_t.fastq --sam 5008_1#2.sam

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

This script adds indexing tag sequence to the sam from fastq file
The order of the short reads should match the order in sam file
SAM file name can be given via command or sam contents passed from stdin

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

=item npg_common::sequence::SAM_Index_Tag

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

