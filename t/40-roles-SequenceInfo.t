use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Moose::Util;
use Moose::Meta::Class;
use Test::More tests => 6;
use Test::Exception;
use Test::Deep;

BEGIN {
  use_ok(q{npg_common::roles::SequenceInfo});
}

my $dna = q{ACgTGAtgCNTgctcacnnntgACGTGCtcgtgTn};

{
  my $test = Moose::Meta::Class->create_anon_class(
    roles => [qw/npg_common::roles::SequenceInfo/],
  )->new_object();

  my @results = $test->indices_of_base($dna, q{A});
  my $expected = [0,5,15,22];

  is_deeply( \@results, $expected, q{expected indices for A, case insensitive} );

  $expected = [17,18,19,34];
  @results = $test->indices_of_base($dna, q{n}, 1);
  is_deeply( \@results, $expected, q{expected indices for n, case sensitive} );

  is( $test->convert_to_quality_score('+'), 10, q{+ converts to a quality score of 10} );
  is( $test->convert_to_quality_score( {
    quality_character => 'B',
    offset => '64',
  }), 2, q{B converts to a quality score of 2 wih offset 64} );
  throws_ok {
    $test->convert_to_quality_score();
  } qr{no[ ]details[ ]provided}, q{no details throws ok};
}
1;
