use strict;
use warnings;
use Test::More tests => 48;
use Test::Exception;
use Test::Deep;
use Perl6::Slurp;
use npg_qc::Schema;
use Moose::Meta::Class;
use npg_testing::db;

use_ok('npg_common::fastqcheck');

my $db_creator = Moose::Meta::Class->create_anon_class(
          roles => [qw/npg_testing::db/])->new_object({});
my $schema;
lives_ok { $schema = $db_creator->create_test_db('npg_qc::Schema') } 'test db created';

{
  my $fq;
  lives_ok { $fq = npg_common::fastqcheck->new(fastqcheck_path => 
                             q{t/data/2549_1_2.fastqcheck}) } q{fq_object created ok};
  isa_ok($fq, q{npg_common::fastqcheck}, q{$fq});
}

{
  throws_ok { 
    npg_common::fastqcheck->new(
                      fastqcheck_path => q[file.fastqcheck],
                      schema => $schema
                               ) 
            } qr/No such file or directory/,  'error when file is not found';

  throws_ok {
    npg_common::fastqcheck->new(
                      fastqcheck_path => q[file.fastqcheck]
                               )
            } qr/No such file or directory/,  'error when file is not found';

  throws_ok {
    npg_common::fastqcheck->new(
                      fastqcheck_path => q[file.fastqcheck],
                      schema => $schema,
                      db_lookup => 0
                               )
            } qr/No such file or directory/,  'error when file is not found';

  throws_ok {
    npg_common::fastqcheck->new(
                      fastqcheck_path => q[file.fastqcheck],
                      schema => $schema,
                      db_lookup => 1
                               )
            } qr/is not in the long-term storage/,  'error when file is not found';
}


{
  throws_ok { npg_common::fastqcheck->new(); } qr/Either fastqcheck_path or file_content should be supplied/,  'error path to fastqcheck is not supplied';
  throws_ok { npg_common::fastqcheck->new(fastqcheck_path => q[file.fastq]); } qr/does not pass the type constraint/,  'error on file with wrong extension';
}



{
  throws_ok { npg_common::fastqcheck->new(fastqcheck_path => q[t/data/wrong.fastqcheck]); } qr/expected structure/,  'error on file with wrong structure';
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/2549_1_2.fastqcheck});

  is($fq->num_reads, 9_855_354, q{number of sequences is correct});
  is($fq->total_pf_bases(), 532_189_116, q{total_pf_bases is correct});
  is($fq->average_read_length(), '54.00', q{average cycle count is correct});
  is($fq->max_read_length(), 54, q{max cycle count is correct});
}

{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/2549_1_2.fastqcheck});
  my $expected = {A => 26.9, C => 22.9, G => 23.3, T=> 26.5,  N => 0.3,};
  cmp_deeply ($fq->base_percentages, $expected, 'total base percentage hash is OK');
}


{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/2549_1_2.fastqcheck});

  throws_ok {$fq->qx_yield(q[sdsdsd])} qr/Wrong type of argument/, 
                     'error on the argument to qx_yield not an array reference';

  my @thresholds = qw(-2 45);
  throws_ok {$fq->qx_yield(\@thresholds)} qr/-2 must be a non-negative integer/, 'error when threshold is negative';

  @thresholds = qw(20 24.5);
  throws_ok {$fq->qx_yield(\@thresholds)} qr/24.5 must be a non-negative integer/, 'error when threshold is a float';

  @thresholds = qw(20 dada);
  throws_ok {$fq->qx_yield(\@thresholds)} qr/dada must be a non-negative integer/, 'error when threshold is a string';
}


{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/2549_1_2.fastqcheck});

  my @thresholds = qw(20 30 35);
  is (join(q[;], @{$fq->qx_yield(\@thresholds)}), q[469992217;235266530;84371445], 'Q20, Q30 and Q35 values computed');

  @thresholds = qw(42);
  is (join(q[;], @{$fq->qx_yield(\@thresholds)}), q[0], 'Q42 value computed');

  @thresholds = qw(38);
   is (join(q[;], @{$fq->qx_yield(\@thresholds)}), q[540843], 'Q38 value computed');

  @thresholds = qw(43);
  is (join(q[;], @{$fq->qx_yield(\@thresholds)}), q[0], 'non-existing Q43 value computed');

}



{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/zero.fastqcheck});
  is($fq->num_reads, 0, 'zero sequences for an empty fastq');
  is($fq->average_read_length, 0, 'zero av cycle count for an empty fastq');
  is($fq->max_read_length, 0, 'zero max cycle count for an empty fastq');
  is($fq->total_pf_bases, 0, 'zero total pf bases count for an empty fastq');
  
  my $expected = {A => 0, C => 0, G => 0, T=> 0,  N => 0,};
  cmp_deeply ($fq->base_percentages, $expected, 'total base percentage hash for an empty fastq');

  my @thresholds = qw(42);
  is (join(q[;], @{$fq->qx_yield(\@thresholds)}), q[0], 'zero returned for an empty fastq, one threshold as an arg');
  
  push @thresholds, 20;
  is (join(q[ ], @{$fq->qx_yield(\@thresholds)}), q[0 0], 'zeros returned for an empty fastq, two thresholds as an arg');
}


{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/test1.fastqcheck});
  my @thresholds = qw(42);
  is($fq->_max_threshold, 42, '_max_threshold correct');
  is (join(q[;], @{$fq->qx_yield(\@thresholds)}), 1617596, 'last available quality OK');
}

{
  my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/test1.fastqcheck});
  throws_ok {$fq->base_percentages(-2)}  qr/Invalid cycle number -2/, 
                            'error when base % for a negative cycle are requested';
  throws_ok {$fq->base_percentages(0)}  qr/Invalid cycle number 0/, 
                            'error when base % for cycle 0 are requested';
  throws_ok {$fq->base_percentages(78)}  qr/Invalid cycle number 78/, 
                            'error when base % for cycle 0 are requested';
  my $expected = {A => 25.2, C => 24.8, G => 23.9, T=> 25.7,  N => 0.3,};
  cmp_deeply ($fq->base_percentages(4), $expected, 'base percentages for cycle 4');
}

{
   my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/test1.fastqcheck});
   is($fq->qx_yield([42], 5)->[0], 1086100, 'last available qval for cycle 5');
}

{
   my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/test1.fastqcheck});
   is(join(q[ ], $fq->gc_percentage()), '46.2 46.5', 'total qc content');
   is(join(q[ ], $fq->gc_percentage(6)), '46.4 46.7', 'qc content for lane 6');
}

{
   my $fq = npg_common::fastqcheck->new(fastqcheck_path => q{t/data/test2.fastqcheck});
   throws_ok {$fq->qx_yield([25])} qr/does not pass the type constraint/, 'error when one of the bin values is negative';
}

{  
  my $fname = q[2549_1_2.fastqcheck];
  my $text = slurp q[t/data/] . $fname;
  lives_ok { $schema->resultset('Fastqcheck')->update_or_create(
             {file_name => $fname, file_content => $text,
              id_run => 2549, position => 1, section => 'reverse'})
           } 'test db updated';

  my $fq;
  throws_ok {$fq = npg_common::fastqcheck->new(
        fastqcheck_path => $fname, schema => $schema)
            } qr/No such file or directory/,  'error when file is not found';
  lives_ok  { $fq = npg_common::fastqcheck->new(
        fastqcheck_path => $fname, schema => $schema, db_lookup => 1)
            } 'object created ok';

  is($fq->num_reads, 9_855_354, q{number of sequences is correct});
  is($fq->total_pf_bases(), 532_189_116, q{total_pf_bases is correct});
  is($fq->average_read_length(), '54.00', q{average cycle count is correct});
  is($fq->max_read_length(), 54, q{max cycle count is correct});
}

1;