#########
# Author:        kl2
# Maintainer:    $Author $
# Created:       2012-10-31
# Last Modified: $Date $
# Id:            $Id $
# $HeadURL $
#

use strict;
use warnings;
use Test::More tests => 10;
use Test::Cmd;
use Test::Exception;
use File::Which;
use File::Temp qw/tempdir/;

my $exec = 'fastq_summ';

SKIP: {

  skip "$exec executable not found", 10, if (!which($exec));

  my $test_cmd = Test::Cmd->new(prog => 'scripts/generate_cached_fastq', workdir => '', );
  #if ($test_cmd) { print "all right 1\n"; } else { print "not all right 1\n"; }

  my $id_run = 7197;
  my $position = 7;
  my $sample_size = 10000;

  my $testdatasrc = q[t/data/scripts/generate_cached_fastq/] . $id_run;
  my $testdir = tempdir( CLEANUP => 1 );
  $testdir = $testdir . q[/];
  my $path = $testdir . $id_run;
  my $file = $path . q[/] . $id_run . q[_] . $position . q[.bam];
  system(q[cp -a ] . $testdatasrc . q[ ] . $testdir) == 0 or die $!;

#print "file: $file\n";
#print "path: $path\n";
  $test_cmd->run(args => "--path $path --file $file --sample_size $sample_size");

  ok($? == 0, 'generate_cached_fastq executed successfully');

  my $orig_fqc = $path . q[/7197_7_1.fastqcheck];
  ok(! -e $orig_fqc, "$orig_fqc should not exist");

  my $fqc = $path . q[/7197_7.fastqcheck];
  ok(-e $fqc && -s $fqc, "$fqc should exist and not be empty");

  my $empty_fqc = $path . q[/lane7/7197_7.fastqcheck];
  ok(-e $empty_fqc && ! -s $empty_fqc, "$empty_fqc should exist and be empty");

  my $fq = $path . q[/.npg_cache_10000/7197_7.fastq.10000];
  ok(-e $fq && -s $fq, "$fq should exist and not be empty");

  # now test se with ic
  $id_run = 7197;
  $position = 7;
  $sample_size = 10000;

  $testdatasrc = q[t/data/scripts/generate_cached_fastq/] . $id_run . q[_se_with_ic/];
  $path = $testdir . $id_run . q[_se_with_ic/];
  $file = $path . $id_run . q[_] . $position . q[.bam];

  system(q[cp -a ] . $testdatasrc . q[ ] . $testdir) == 0 or die $!;

  $test_cmd->run(args => "--path $path --file $file --sample_size $sample_size");

  ok($? == 0, 'generate_cached_fastq executed successfully');

  $orig_fqc = $path . q[/7197_7_1.fastqcheck];
  ok(-e $orig_fqc && -s $orig_fqc, "$orig_fqc should exist and not be empty");

  $fqc = $path . q[/7197_7.fastqcheck];
  ok(! -e $fqc, "$fqc should not exist");

  $empty_fqc = $path . q[/lane7/7197_7.fastqcheck];
  ok(-e $empty_fqc && ! -s $empty_fqc, "$empty_fqc should exist and be empty");

  $fq = $path . q[/.npg_cache_10000/7197_7_1.fastq.10000];
  ok(-e $fq && -s $fq, "$fq should exist and not be empty");

  system(q[rm -fr ] . $testdir) == 0 or die $!;

}

1;

