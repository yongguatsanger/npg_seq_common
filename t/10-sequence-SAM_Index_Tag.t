#########
# Author:        gq1
# Created:       2010-08-03
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 4;
use Test::Exception;

use_ok('npg_common::sequence::SAM_Index_Tag');

{
  my $fastq_t = 't/data/sequence/SecondCall/4394_1_t.fastq';
  my $sam_file = 't/data/sequence/SecondCall_paired/4394_1.sam';
  my $sam = npg_common::sequence::SAM_Index_Tag->new({
                 sam                => $sam_file,
                 index_fastq_file   => $fastq_t,
               });
               
  isa_ok($sam, 'npg_common::sequence::SAM_Index_Tag', 'object test');
  
  open my $fastq_t_fh, q{<}, $fastq_t or croak "Can not open file $fastq_t";
  my $sam_line = q{IL37_4394:1:1:1115:7227	77	*	0	0	*	*	0	0	TCATAGNATTTCTTGTAATTCTCCTTGTCTTCCACCAGCTCAGTGAAGAGCTCA	III@EI$DCAA?A?AAFCC=BDEF9C@CI:IFII?D@7ADF+CAC?,=B@D:DA};
  is($sam->_add_index_tag($sam_line, $fastq_t_fh), $sam_line.'	RT:Z:GATCCGAT	QT:Z:,,A#A(?#', 'correct sam line returned with indexint tag sequence');
  
  lives_ok{$sam->process();} q{processing ok};

}

1;
