use strict;
use warnings;
use Test::More tests => 38;
use Test::Exception;
use Test::Deep;
use lib q{t};
use t::util;

my $util = t::util->new();
my $dir = $util->temp_directory();
my $id_run = 1008;
my $tag_index;
my $position = 1;

my @subs = qw/ generate_cache retrieve_from_cache/;
use_ok( 'npg_common::extractor::fastq', @subs);
foreach my $sub (@subs) {
  can_ok(__PACKAGE__, $sub);
}

{
  is(npg_common::extractor::fastq::_cache_dir(q[t], 3), q[t/.npg_cache_3], 'cache dir name');
  is(npg_common::extractor::fastq::_cache_dir(q[t]), q[t/.npg_cache_10000], 'cache dir name');
  local $ENV{NPG_FASTQ_CACHE} = q[t];
  is(npg_common::extractor::fastq::_cache_dir(q[t], 3), q[t], 'cache dir name');
}

{
  my $u = t::util->new();
  my $d = $u->temp_directory();
  my $file = $d . q[/1008_1.fastq];
  `cp t/data/1008_1_1.fastq $file`;

  my $cache = join(q[/], $d, q[.npg_cache_5]);
  generate_cache($d, [$file], 5);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1.fastq.3]), 'cached single exists');

  $cache = join(q[/], $d, q[.npg_cache_10000]);
  generate_cache($d, [$file]);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1.fastq.3]), 'cached single exists');
}

{
  my $u = t::util->new();
  my $d = $u->temp_directory();
  my $file = $d . q[/1008_1#1.fastq];
  `cp t/data/1008_1_1.fastq $file`;
 
  my $cache = join(q[/], $d, q[.npg_cache_5]);
  generate_cache($d, [$file], 5);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1#1.fastq.3]), 'cached single tag exists');

  $cache = join(q[/], $d, q[.npg_cache_10000]);
  generate_cache($d, [$file]);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1#1.fastq.3]), 'cached single tag exists');
}

{
  `cp t/data/1008_1_1.fastq $dir`;
  `cp t/data/1008_1_2.fastq $dir`;
 
  my $cache = join(q[/], $dir, q[.npg_cache_5]);
  generate_cache($dir, ["$dir/1008_1_1.fastq", "$dir/1008_1_2.fastq"], 5);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1_1.fastq.3]), 'cached forward exists');
  ok (-e join(q[/], $cache, q[1008_1_2.fastq.3]), 'cached reverse exists');

  $cache = join(q[/], $dir, q[.npg_cache_10000]);
  generate_cache($dir, ["$dir/1008_1_1.fastq", "$dir/1008_1_2.fastq"]);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1_1.fastq.3]), 'cached forward exists');
  ok (-e join(q[/], $cache, q[1008_1_2.fastq.3]), 'cached reverse exists');
  `rm -rf $cache`;
}

{
  my $file = $dir . q[/1008_1#1.fastq];
  `cp t/data/1008_1_1.fastq $file`;
  my $filet = $dir . q[/1008_1_t.fastq];
  `cp t/data/1008_1_1.fastq $filet`;
  my $cache = join(q[/], $dir, q[.cache]);
  local $ENV{NPG_FASTQ_CACHE} = $cache;
  generate_cache(q[t/data], ["$dir/1008_1_1.fastq", "$dir/1008_1_2.fastq", $file, $filet], 2);
  ok (-d $cache, 'cache directory exists');
  ok (-e join(q[/], $cache, q[1008_1_1.fastq.2]), 'cached forward exists');
  ok (-e join(q[/], $cache, q[1008_1_2.fastq.2]), 'cached reverse exists');
  ok (-e join(q[/], $cache, q[1008_1#1.fastq.2]), 'cached plex exists');
  ok (-e join(q[/], $cache, q[1008_1_t.fastq.2]), 'cached t exists');   
}

{
  my $cache = q[.npg_cache_5];
  my $f = join(q[/], $dir, q[1008_1_1.fastq]);
  is(retrieve_from_cache($f, 5), join(q[/], $dir, $cache, q[1008_1_1.fastq.3]), 'file retrieved from cache');
  
  my $c10000 = $dir . q[/.npg_cache_10000];
  throws_ok {retrieve_from_cache($f)} qr/Cache directory $c10000 does not exist/, 'error when cache does not exist';
  mkdir $c10000;
  is(retrieve_from_cache($f), undef, 'nothing retrieved from cache');
  my $cf10000 = join q[/], $c10000, q[1008_1_1.fastq.5];
  `touch $cf10000`;
  is(retrieve_from_cache($f), join(q[/], $c10000, q[1008_1_1.fastq.5]), 'file retrieved from cache');

  my $f_phix = join(q[/], $dir, q[1008_1_1_phix.fastq]);
  `touch $f_phix`;
  sleep 1;
  is(retrieve_from_cache($f_phix), undef, 'phix file not retrieved from cache');

  sleep 1;
  `touch $f`;
  is(retrieve_from_cache($f, 5), undef, 'nothing retrieved from cache - file old');

  $cache = join(q[/], $dir, q[.cache]);
  my $cached = join(q[/], $cache, q[1008_1_2.fastq.2]);
  local $ENV{NPG_FASTQ_CACHE} = $cache;
  is(retrieve_from_cache(join(q[/], $dir, q[1008_1_2.fastq]), 2), $cached, 'file retrieved from cache');

  sleep 1;
  `touch $cf10000`;
}

{
  my $cache = q[.npg_cache_5];
  my $f = join(q[/], $dir, q[1008_1_2.fastq]);
  my $new_file = join(q[/], $dir, q[new_file]);
  is(npg_common::extractor::fastq::_link2cache($f, $new_file, 5), 3, 'actual sample size');
  ok (-l $new_file, 'symbolic link created');

  $f = join(q[/], $dir, q[1008_1_1.fastq]);
  $new_file = join(q[/], $dir, q[new_file1]);
  is(npg_common::extractor::fastq::_link2cache($f, $new_file), 5, 'actual sample size');
  ok (-l $new_file, 'symbolic link created');
}

{
  my $u = t::util->new();
  my $d = $u->temp_directory();
  my $file = $d . q[/1008_1.bam];
  `touch $file`;

  my $cache = join(q[/], $d, q[.npg_cache_5]);
  local $ENV{PATH}="$d";
  throws_ok { generate_cache($d, [$file], 5); } qr/fastq_summ failed for file:/, q[testing generate_cache with bam file $file];
}

{
  my $u = t::util->new();
  my $d = $u->temp_directory();
  my $file = $d . q[/1008_1.bam];
  `touch $file`;

  `touch $d/fastq_summ; chmod +x $d/fastq_summ;`;

  my $cache = join(q[/], $d, q[.npg_cache_5]);
  local $ENV{PATH}="$d";
  lives_ok { generate_cache($d, [$file], 5); } q[testing generate_cache with bam file $file];
}

