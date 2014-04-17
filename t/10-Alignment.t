#########
# Author:        gq1
# Created:       2009-04-17
#

use strict;
use warnings;
use Test::More tests => 4;
use Test::Exception;
use File::Temp qw/tempdir/;

use_ok('npg_common::Alignment');

{
  my $bwa = npg_common::Alignment->new(bwa_cmd => q{t/bin/bwa});
               
  isa_ok($bwa, 'npg_common::Alignment', 'object test');

  my $reference_bwa = 't/data/references/Human/NCBI36/all/bwa/Homo_sapiens.NCBI36.48.dna.all.fa';
  my $query_fastq = 't/data/fastq_split/2605_1_1.fastq';
  my $sam = join q[/], tempdir( CLEANUP => 1 ), q[2605_1_1.sam];

  lives_ok {$bwa->bwa_align_se({fastq=>$query_fastq, bam_out=>$sam, ref_root=>$reference_bwa})}
    "no croak for bwa alignment single";

  lives_ok {$bwa->bwa_align_pe({fastq1=>$query_fastq, fastq2=>$query_fastq, sam_out=>$sam, fork_align=>1, ref_root=>$reference_bwa})}
    "no croak for bwa alignment paired ";
}

1;
