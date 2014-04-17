#########
# Author:        Marina Gourtovaia
# Created:       27 January 2010
#

use strict;
use warnings;
use Test::More tests => 18;
use Test::Exception;
use Test::Deep;
use File::Temp qw/ tempdir /;

use_ok('npg_common::sequence::reference::base_count');
my $dir = tempdir( CLEANUP => 1 );

my $ref0 =  join q[/], $dir, q[simple_reference.ref];
open my $fh0, '>', $ref0 or die "cannot open $ref0";
my @l0 = qw/ 
>1
jfsdjfsdlkfj
dfsdfsdfsfsf
sdfsdfdsfsdf
sdsfsfsfsfs
sfsfsfsfsfsfd
sdfsffsfsfsf
>20
aaaaaaaaaaaaa
bbbbbbbbbbbbb
>21
asdasdasdadad
/;

print $fh0 join("\n", @l0);
close $fh0;

{
  my $bc = npg_common::sequence::reference::base_count->new(reference_path => $ref0);
  isa_ok($bc, 'npg_common::sequence::reference::base_count');
  throws_ok {npg_common::sequence::reference::base_count->new()}
             qr/is required/, 'error when path to ref is not set in the constructor';
  throws_ok { npg_common::sequence::reference::base_count->new(reference_path => q[t/data/simple])->run() } 
             qr/does not exist or is not readable/, 'error when path to ref does not exist';
  lives_ok {$bc = npg_common::sequence::reference::base_count->new(
                            reference_path => $ref0)} 'constractor with a ref path only lives';
  is ($bc->reference_path, $ref0, 'ref path set correctly');
  lives_ok {$bc->summary} 'summary call after creating an object lives';
  throws_ok {$bc->run()} qr/contains unexpected character/  , 'error when a bad char present';
}

{
  my $expected = {
          '100' => 'D', '103' => 'G', '104' => 'H', '107' => 'K', '109' => 'M', '110' => 'N', 
          '114' => 'R', '115' => 'S', '116' => 'T', '118' => 'V', '119' => 'W', '120' => 'N', 
          '121' => 'Y', '46'  => '',  '65'  => 'A', '66'  => 'B', '67'  => 'C', '68'  => 'D', 
          '71'  => 'G', '72'  => 'H', '75'  => 'K', '77'  => 'M', '78'  => 'N', '82'  => 'R', 
          '83'  => 'S', '84'  => 'T', '86'  => 'V', '87'  => 'W', '88'  => 'N', '89'  => 'Y',
          '95'  => '',  '97'  => 'A', '98'  => 'B', '99'  => 'C',};

  cmp_deeply (npg_common::sequence::reference::base_count->translation_table(), $expected, 'translation table class method ok');
  is (join(q[ ], @{npg_common::sequence::reference::base_count->output_base_codes()}),
                     q[A C G T M R W S Y K V H D B N], 'output codes returned');
}

{
  my $empty = join q[/], $dir, q[empty.ref];
  `touch $empty`;
  my $bc = npg_common::sequence::reference::base_count->new(reference_path => $empty,);
  lives_ok {$bc->run()} 'runs on empty ref';
  cmp_deeply ($bc->summary, {counts => {}, ref_length => 0,}, 'summary shows zero length and no chars');
}

my $ref =  join q[/], $dir, q[no_seq_name.ref];
open my $fh, '>', $ref or die "cannot open $ref";
print $fh "abctdaa\naCtgaGgan\n";
close $fh;

{
  my $ref =  join q[/], $dir, q[no_seq_name.ref];
  open my $fh, '>', $ref or die "cannot open $ref";
  print $fh "abctdaa\naCtgaGgan\n";
  close $fh;

  my $bc = npg_common::sequence::reference::base_count->new(reference_path => $ref);
  lives_ok {$bc->run()} 'run method OK for a sequence without header';
  my $summary = {counts => {'A' => 6, 'B' => 1, 'C' => 2, 'T' => 2, 'D' => 1, 'G' => 3, 'N' => 1,}, ref_length => 16,};
  cmp_deeply ($bc->summary, $summary, 'summary is correct');
}

{
  my $ref1 =  join q[/], $dir, q[three_seq.ref];
  open my $fh1, '>', $ref1 or die "cannot open $ref1";
  print $fh1 ">chr1\n";
  print $fh1 "abctdaa\n";
  print $fh1 "aCtgaGgan\n";
  print $fh1 ">chr2\n";
  print $fh1 "xaadWS\n";
  print $fh1 ">chr3\n";
  print $fh1 "dartYX\n";
  close $fh1;

  my $bc = npg_common::sequence::reference::base_count->new(reference_path => $ref1);
  $bc->run();
  my $summary = {counts => {'A' => 9, 'B' => 1, 'C' => 2, 'T' => 3, 'D' => 3, 'G' => 3, 'N' => 3, 'W' => 1, 'S' => 1, 'R' => 1, 'Y' => 1,}, ref_length => 28,};
  cmp_deeply ($bc->summary, $summary, 'summary is correct');

  my $frozen;
  lives_ok {$frozen = $bc->freeze()} 'freezing the object to json lives';
  my $bc2 = npg_common::sequence::reference::base_count->thaw($frozen);
  cmp_deeply ($bc, $bc2, 'object from a frozen json string looks like the original object');
}

{
  my $ref2 =  join q[/], $dir, q[with_gaps.ref];
  open my $fh2, '>', $ref2 or die "cannot open $ref2";
  print $fh2 ">chr1\n";
  print $fh2 "_abc__tdaa\n";
  print $fh2 ".aCtgaG_ganga_\n";
  close $fh2;

  my $bc = npg_common::sequence::reference::base_count->new(reference_path => $ref2);
  $bc->run();
  my $summary = {counts => {'A' => 7, 'B' => 1, 'C' => 2, 'T' => 2, 'D' => 1, 'G' => 4, 'N' => 1,}, ref_length => 24,};
  cmp_deeply ($bc->summary, $summary, 'summary is correct');
}

1;
