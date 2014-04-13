#########
# Author:        gq1
# Created:       2009-04-17
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 7;
use npg_common::Hit_Sequence_SAM;

use_ok('npg_common::Hit_Sequence_SAM');
  
{
  my $hit_sequence_sam = npg_common::Hit_Sequence_SAM->new({
                                        filename =>'t/data/fastq_split/2518_2_5_merge_part.sam',
                                        hit_ids  => {},
                                        nothit_ids =>{},
                                        });

  isa_ok($hit_sequence_sam, 'npg_common::Hit_Sequence_SAM', 'object test');

  eval{
    $hit_sequence_sam->parsing_file();
    1;
  };
  is($EVAL_ERROR, q{}, 'no croak when parsing sam file');

  my $hit_ids = $hit_sequence_sam->hit_ids();
  isa_ok($hit_ids, 'HASH', 'hit_ids');
  is ($hit_sequence_sam->num_sequences_hit(), 12, 'correct number hit short reads');

  my $non_hit_ids = $hit_sequence_sam->non_hit_ids();
  isa_ok($non_hit_ids, 'HASH', 'non_hit_ids');
  is($hit_sequence_sam->num_sequences_nonhit(), 14, 'correct number non hit short reads');
}

1;
__END__
