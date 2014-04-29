use strict;
use warnings;
use Test::More tests => 69;
use English qw(-no_match_vars);
use Test::Exception;
use File::Temp qw/tempdir/;
use Cwd qw/cwd/;
use File::Compare;

my $source1 = q[t/data/s1.fastq];
my $source2 = q[t/data/s2.fastq];
my $dir = tempdir(CLEANUP => 1);

my $empty = $dir . q[/empty_fastq.fastq];
`touch $empty`;

my @subs = qw(read_count generate_equally_spaced_reads split_reads to_fasta);
use_ok( 'npg_common::extractor::fastq', @subs);

foreach my $sub (@subs) {
  can_ok(__PACKAGE__, $sub);
}

{
  my $count = read_count(q[t/data/ss1.fastq]);
  is($count, '26', 'read count when fastqcheck file is present');

  $count = read_count($source2);
  is($count, '26', 'read count when fastqcheck file is not present');  
}

{
  my $l1 = q[@IL36_1937:1:1:1416:928/1];
  my $l2 = q[@IL36_1937:1:1:1416:928/2];
  ok(npg_common::extractor::fastq::_is_pair($l1, $l2), 'read names the same');
  
  $l2 = q[@IL36_1937:1:1:1755:1159/1];
  ok(!npg_common::extractor::fastq::_is_pair($l1, $l2), 'read names different');
}

{
   # 26 reads in test file
   my $dest1   = $dir . q[/sXe.fastq];
   my $dest2   = $dir . q[/sX2e.fastq];
   my $num_reads_required = 33;

   is (generate_equally_spaced_reads([$source1, $source2], [$dest1, $dest2], $num_reads_required), 
           26, 'correct num. reads returned');

   ok (compare($dest1, $source1)==0, 'output as expected for the first file');
   ok (compare($dest2, $source2)==0, 'output as expected for the second file');   
}

{
   my $dest1   = $dir . q[/sXe.fastq];
   my $dest2   = $dir . q[/sX2e.fastq];
   my $num_reads_required = 26;

   is (generate_equally_spaced_reads([$source1, $source2], [$dest1, $dest2], $num_reads_required), 
              $num_reads_required, 'correct num. reads returned');
   ok (compare($dest1, $source1) == 0, 'output as expected for the first file');
   ok (compare($dest2, $source2) == 0, 'output as expected for the second file');   
}

{
   my $dest1   = $dir . q[/sZe.fastq];
   my $dest2   = $dir . q[/sZ2e.fastq];
   my $source3 = q[t/data/unsorted.fastq];
   my $num_reads_required = 13;
   throws_ok {generate_equally_spaced_reads([$source1, $source3], [$dest1, $dest2], $num_reads_required)} 
            qr/Reads are out of order/, 'error if fastq files are unsorted';
}

{   
   my $dest1   = $dir . q[/sYe.fastq];
   my $dest2   = $dir . q[/sY2e.fastq];
   my $test = q[t/data/test.fastq];
   my $num_reads_required = 3;

   is (generate_equally_spaced_reads([$source1, $source2], [$dest1, $dest2], $num_reads_required), 
                 $num_reads_required, 'correct num. reads returned');
   ok (compare($dest1, $test)==0, 'output as expected for the first file');
   ok (compare($dest2, $test)==0, 'output as expected for the second file');
}


{
   my $s1 = q[t/data/ss1.fastq];
   my $s2 = q[t/data/ss2.fastq];
   my $dest1   = $dir . q[/tdtd1.fastq];
   my $dest2   = $dir . q[/tdtd2.fastq];
   my $test = q[t/data/test.fastq];
   my $num_reads_required = 3;
   generate_equally_spaced_reads([$s1, $s2], [$dest1, $dest2], $num_reads_required);

   ok(-e $dest1, 'output for read1 exists');
   ok(-e $dest2, 'output for read2 exists');
   ok(compare($dest1, $test)==0, 'first file extraction when the read count is taken from a fastqcheck file');
   ok(compare($dest2, $test)==0, 'second file extraction when the read count is taken from a fastqcheck file'); 
}

{
   my $s1 = q[t/data/ss1.fastq];
   my $s2 = q[t/data/ss2.fastq];
   my $s3 = q[t/data/ss3.fastq];
   my $dest1   = $dir . q[/tdta1.fastq];
   my $dest2   = $dir . q[/tdta2.fastq];
   my $dest3   = $dir . q[/tdta3.fastq];
   my $test = q[t/data/test.fastq];
   my $test3 = q[t/data/test3.fastq];
   my $num_reads_required = 3;
   generate_equally_spaced_reads([$s1, $s2, $s3], [$dest1, $dest2, $dest3], $num_reads_required);

   ok(-e $dest1, 'output for read1 exists');
   ok(-e $dest2, 'output for read2 exists');
   ok(-e $dest3, 'output for readX exists');
   ok(compare($dest1, $test)==0, 'first file extraction when the read count is taken from a fastqcheck file');
   ok(compare($dest2, $test)==0, 'second file extraction when the read count is taken from a fastqcheck file');
   ok(compare($dest3, $test3)==0, 'another file extraction'); 
}

{
   my $s1 = q[t/data/sss1.fastq];
   my $s2 = q[t/data/sss2.fastq];
   my $dest1 = $dir . q[/sLe.fastq];
   my $dest2 = $dir . q[/sL2e.fastq];
   my $num_reads_required = 3;
   my $test = q[t/data/test.fastq];

   generate_equally_spaced_reads([$s1, $s2], [$dest1, $dest2], $num_reads_required);
   
   ok(compare($dest1, $test)==0, 'first file extraction when the read count in a fastqcheck file is an illegal number');
   ok(compare($dest2, $test)==0, 'second file extraction when the read count in a fastqcheck file is an illegal number'); 
}

{
   my $s1 = q[t/data/ss1.fastq];
   my $dest1 = $dir . q[/sKe.fastq];
   my $num_reads_required = 3;
   my $test = q[t/data/test.fastq];

   generate_equally_spaced_reads([$s1], [$dest1], $num_reads_required);
   ok(compare($dest1, $test)==0, 'single file extraction when the read count is taken from a fastqcheck file');
}

{
   is (generate_equally_spaced_reads([$empty], ['t/data/empty1'], 5), 0, 'zero length for an empty file');
}

{
   my $s1 = $dir . q[/trunkated.fastq];
   `head -n 8 t/data/ss1.fastq > $s1`;
   `cp t/data/ss1.fastqcheck $dir/trunkated.fastqcheck`; 
   my $dest1 = $dir . q[/trunkated_dest.fastq];
   my $num_reads_required = 26;
   throws_ok {generate_equally_spaced_reads([$s1], [$dest1], $num_reads_required)} qr/is shorter than reported 26 reads/, 'error when the source file is trunkated';
}

{
   my $s1 = $dir . q[/trunkated.fastq];
   `head -n 7 t/data/ss1.fastq > $s1`;
   `cp t/data/ss1.fastqcheck $dir/trunkated.fastqcheck`; 
   my $dest1 = $dir . q[/trunkated_dest.fastq];
   my $num_reads_required = 13;
   throws_ok {generate_equally_spaced_reads([$s1], [$dest1], $num_reads_required)} qr/is shorter than reported 26 reads/, 'error when the source file is trunkated';
}

{
   throws_ok {split_reads()} qr/Input file name should be given/, 'error when no args are given';
   throws_ok {split_reads($empty)} qr/Read lengths for the target files should be given as an array reference/, 'error when no target read length is given';
   throws_ok {split_reads($empty, [0])} qr/Target read length should be positive/, 'error when zero target read length is given';
   throws_ok {split_reads($empty, [55, 0])} qr/Target read length should be positive/,'error when zero target read length is given';
   throws_ok {split_reads($empty, [55, -55])} qr/Target read length should be positive/,'error when negative target read length is given';  
}

{
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($empty, [37, 37])} 'splitting an empty file leaves';
   ok(-e q[37_1.fastq] && -e q[37_2.fastq], 'output for splitting an empty file exists');
   chdir $current;

   my $out1 = $dir . q[/dodo1];
   my $out2 = $dir . q[/dodo2];
   lives_ok {split_reads($empty, [37, 37], [$out1, $out2])} 'splitting an empty file leaves';
   ok(-e $out1 && -e $out2, 'output for splitting an empty file exists');
}

{
   my $out1 = $dir . q[/dodo3];
   my $out2 = $dir . q[/dodo4];
   lives_ok {split_reads($empty, [], [$out1, $out2])} 'splitting an empty file in halves leaves';
   ok(-e $out1 && -e $out2, 'output for splitting an empty file exists');
}

{
   my $current = cwd;
   chdir $dir;
   throws_ok {split_reads($current . q[/t/data/s1.fastq], [37])} qr/is too short/, 'error when the read is too short';
   chdir $current;
}

{  
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($current . q[/t/data/ss1.fastq], [3])} 'splitting reads lives';
   ok(-e q[3_1.fastq], 'output exists');
   ok(!-e q[3_2.fastq], 'second output does not exist');
   chdir $current;
   ok(compare(q[t/data/ss1_split3.fastq], $dir . q[/3_1.fastq])==0, 'correct output for a single file'); 
}

{  
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($current . q[/t/data/1008_s_1.fastq], [37,37])} 'splitting reads in two with explicit lengths lives';
   ok(-e q[37_1.fastq] && -e q[37_2.fastq], 'output exists');
   chdir $current;
   ok(compare(q[t/data/1008_1_1.fastq], $dir . q[/37_1.fastq])==0, 'correct output for a forward file');
   ok(compare(q[t/data/1008_1_2.fastq], $dir . q[/37_2.fastq])==0, 'correct output for a reverse file'); 
}

{  
   my $current = cwd;
   chdir $dir;
   my $f1 = q[half1.fastq];
   my $f2 = q[half2.fastq];
   lives_ok {split_reads($current . q[/t/data/1008_s_1.fastq], [], [$f1,$f2])} 'splitting reads in halves';
   ok(-e $f1 && -e $f2, 'output exists');
   chdir $current;
   ok(compare(q[t/data/1008_1_1.fastq], $dir . q[/] . $f1)==0, 'correct output for a forward file');
   ok(compare(q[t/data/1008_1_2.fastq], $dir . q[/] . $f2)==0, 'correct output for a reverse file'); 
}

{
  lives_ok {to_fasta(q[t/data/1008_1_1.fastq])} 'fasta out to stdout lives';
  my $out = join(q[/], $dir, q[out.fasta]);
  lives_ok {to_fasta(q[t/data/1008_1_1.fastq], $out)} 'fasta out to a file lives';
  ok( -e $out, 'output exists');
  open my $fh, q[<], $out or die $!;
  my @lines = <$fh>;
  close $fh;
  is (scalar @lines, 6, 'six line fasta file in teh output');
  my $line = $lines[0];
  $line =~ s/\s+$//;
  is ($line, q[>@IL14_1008:1:1:470:276/1], 'first fasta line');
  my $expected_line = q[AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA];
  $line = $lines[3];
  $line =~ s/\s+$//;
  is ($line, $expected_line, '4th fasta line');
  $line = $lines[5];
  $line =~ s/\s+$//;
  is ($line, $expected_line, '6th fasta line');

  my $out_empty = join(q[/], $dir, q[outempty.fasta]);
  lives_ok {to_fasta($empty, $out_empty)} 'fasta on an empty file lives';
  ok(-e $out_empty, 'output file exists');
  ok(-z $out_empty, 'output file is empty');
}

1;
