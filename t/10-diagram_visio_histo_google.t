use strict;
use warnings;
use English;
use Test::More tests => 31;
use Test::Exception;
use Test::Deep;

use_ok( 'npg_common::diagram::visio_histo_google');

{
  my $dia = npg_common::diagram::visio_histo_google->new();

  isa_ok($dia, 'npg_common::diagram::visio_histo_google');

  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=300x250&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y', 'default diagram string');
}


{
  my $dia = npg_common::diagram::visio_histo_google->new();
  my $data = [2,3,4];
  $dia->set_data($data);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:2,3,4&chds=0,40&chs=300x250&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y', 'set_data');
}

{
  my $dia = npg_common::diagram::visio_histo_google->new();
  $dia->set_axisY_min_max(2,6);
  is($dia->get_diagram_string(), q[http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:20,30,40,23&chds=2,6&chs=300x250&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y], 'set_axisY_min_max');
}

{
  my $dia = npg_common::diagram::visio_histo_google->new();
  $dia->set_chart_title(q[my chart]);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=300x250&cht=bvg&chtt=my+chart&chxr=0,1,4|1,0,40,10&chxt=x,y', 'set_title');
}


{
  my $dia = npg_common::diagram::visio_histo_google->new();
  $dia->set_chart_labels(2,5,1,4,7,1);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=300x250&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,2,5,1|1,4,7,1&chxt=x,y', 'set labels');
}


{
  my $dia = npg_common::diagram::visio_histo_google->new();
  $dia->set_chart_size(900,200);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=5,1,1&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=900x200&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y', 'set size');
}

{
  my $dia = npg_common::diagram::visio_histo_google->new();
  $dia->set_chart_size(900,200);
  $dia->set_bar_size(7, 2);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=7,1,2&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=900x200&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y', 'set bar size');
  $dia->set_bar_size(8);
  is($dia->get_diagram_string(), 'http://chart.apis.google.com/chart?chbh=8,1,0&chco=4D89F9&chd=t:20,30,40,23&chds=0,40&chs=900x200&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y', 'set bar size');
use POSIX;
  is( $dia->_encode_value( '0.50' ), q{AU}, q{encode value of 0.5 returns AU with default y_max} );
  $dia->y_max( 1000 );
  is( $dia->_encode_value( '0.01' ), q{AB}, q{encode value of 0.01 returns AB with y_max of 1000} );
  $dia->y_max( q{0.575} );
  is( $dia->_encode_value( '0.50' ), q{3o}, q{encode value of 0.5 returns 3o with y_max of 0.575} );
  is( $dia->_encode_value( '0.575' ), q{..}, q{encode value of 0.575 returns AU with y_max of 0.575} );
  is( $dia->_encode_value( '0.58' ), q{..}, q{encode value of 0.58 returns AU with y_max of 0.575} );
  throws_ok { $dia->_encode_value( '0.59' ) } qr/Value to encode 0\.59 is much bigger than max y 0\.575 \(tolerance 0\.01\)/, q{encode value of 0.59 with y_max of 0.575 throws an error};

  $dia->y_max( q{0.57} );
  is( $dia->encode() , 0, q{encode is default 0} );
  is( $dia->encode( 1 ), 1, q{encode is 1} );

  my $data_set = [0.00,0.14,0.06,0.00,0.57];
  my $ref_data_set = ref $data_set;
  my $test_no_change;
  @{ $test_no_change } = @{ $data_set };
  lives_ok { $dia->set_data( $data_set ); } q{set_data ran ok - single array data};
  my $expected_chd_chart_data = ['chd','e:AAPtGvAA..'];

  is_deeply( $dia->{chd_chart_data}, $expected_chd_chart_data, q{single array chd_chart_data is produced correctly} );
  is_deeply( $data_set, $test_no_change, q{no change to original data set - data} );
  is_deeply( ref $data_set, $ref_data_set, q{no change to original data set - memory location} );

  $data_set = [[0.00,0.14,0.06,0.00,0.57],[0.10,0.17,0.32,0.01,0.00],[0.00,0.16,0.18,0.00,0.00]];
  $expected_chd_chart_data = ['chd','e:AAPtGvAA..,LOTFj6BHAA,AAR9UNAAAA'];
  lives_ok { $dia->set_data( $data_set ); } q{set_data ran ok - multiple array data};
  is_deeply( $dia->{chd_chart_data}, $expected_chd_chart_data, q{multiple array chd_chart_data is produced correctly} );

  $data_set = {
    10 => [0.00,0.14,0.06,0.00,0.57],
    20 => [0.00,0.16,0.18,0.00,0.00],
    30 => [0.10,0.17,0.32,0.01,0.00],
  };
  $expected_chd_chart_data = ['chd','e:AAPtGvAA..,AAR9UNAAAA,LOTFj6BHAA'];
  lives_ok { $dia->set_data( $data_set ); } q{set_data ran ok - multiple array data};
  is_deeply( $dia->{chd_chart_data}, $expected_chd_chart_data, q{hash multiple array chd_chart_data is produced correctly} );

  is_deeply( $dia->chdl_chart_legend( q{foo} ), [ q{chdl}, q{foo} ], q{chdl_chart_legend returns correct with string} );
  is_deeply( $dia->chdl_chart_legend( { foo => 1, bar => 1, baz=> 1 } ), [ q{chdl}, q{bar|baz|foo} ], q{chdl_chart_legend returns correct with hashref} );
  is_deeply( $dia->chdl_chart_legend( [qw{ >30 <30 <15 N }] ), [ q{chdl}, q{>30|<30|<15|N} ], q{chdl_chart_legend returns correct with arrayref} );

  $dia->{chco_chart_colour} = ['chco', '0812F7,31F613,ECF21C,DB4400'];
  $dia->{chd_chart_data} = ['chd','e:AAPtGvAA..,AAR9UNAAAA,LOTFj6BHAA,HeLLoWOrlD'];
  is( $dia->get_diagram_string(), q{http://chart.apis.google.com/chart?chdl=>30|<30|<15|N&chbh=8,1,0&chco=0812F7,31F613,ECF21C,DB4400&chd=e:AAPtGvAA..,AAR9UNAAAA,LOTFj6BHAA,HeLLoWOrlD&chds=0,40&chs=900x200&cht=bvg&chtt=batch+X,run+Y,+lane+Z&chxr=0,1,4|1,0,40,10&chxt=x,y}, q{histogram produced with a legend if data for one is present} );
  is( $dia->get_diagram_string(1), q{http://chart.apis.google.com/chart?chdl=>30|<30|<15|N&chco=0812F7,31F613,ECF21C,DB4400&cht=bvg&chs=70x250}, q{only legend produced} );

}




    