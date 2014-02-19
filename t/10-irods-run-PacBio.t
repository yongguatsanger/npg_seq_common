#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-02-08
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 15;
use Test::Exception;
use File::Temp qw/ tempfile tempdir /;

use_ok('npg_common::irods::run::PacBio');

{
  my $pacbio = npg_common::irods::run::PacBio->new( runfolder_path=>'t/data/pacbio/superfoo/02FEB11_TestGenomes2_121' );

  my $bas1 = q{t/data/pacbio/superfoo/02FEB11_TestGenomes2_121/A01_1/Analysis_Results/m110202_150158_00127_c000055432550000000115027204131180_s1_p0.bas.h5};  
  my $bas2 = q{t/data/pacbio/superfoo/02FEB11_TestGenomes2_121/B01_1/Analysis_Results/m110202_161331_00127_c000055432550000000115027204131181_s1_p0.bas.h5};
  my $bas_h5_list = [$bas1, $bas2];

  
  is_deeply ($pacbio->file_list(), $bas_h5_list, 'correct bas h5 file list');

  my $meta1 = q{t/data/pacbio/superfoo/02FEB11_TestGenomes2_121/A01_1/m110202_150158_00127_c000055432550000000115027204131180_s1_p0.metadata.xml}; 
  is ( $pacbio->_get_meta_xml_from_bas_file($bas1), $meta1,'correct metadata xml file path');
  
  my $meta = $pacbio->read_meta_xml($meta1);
  is($meta->{run}, '02FEB11_TestGenomes2', 'correct run name');
  is($meta->{sample}, '2kb lambda (LTS)', 'correct sample name');
  is($meta->{well}, 'A01', 'correct well name');
  is($meta->{instrument_name}, '00127', 'correct instrument name');
  is($meta->{collection_number}, '1', 'correct collection number');
  is($meta->{cell_index}, '0', 'correct cell index');
  is($meta->{set_number}, '0', 'correct set number');
  

  is($pacbio->_get_des_dir($bas2), '/seq/pacbio/02FEB11_TestGenomes2_121/B01_1/Analysis_Results/', 'correction irods destination collection name for bas file');
  is($pacbio->_get_des_dir($meta1), '/seq/pacbio/02FEB11_TestGenomes2_121/A01_1/', 'correction irods destination collection name for meta xml file');
  
  #$pacbio->process();
}

{
  my $pacbio = npg_common::irods::run::PacBio->new( runfolder_path    =>'t/data/pacbio/superfoo/02FEB11_TestGenomes2_121',
                                                    collection_folder => 'A01_1',
   );

  my $bas1 = q{t/data/pacbio/superfoo/02FEB11_TestGenomes2_121/A01_1/Analysis_Results/m110202_150158_00127_c000055432550000000115027204131180_s1_p0.bas.h5};  

  my $bas_h5_list = [$bas1];
  
  is_deeply ($pacbio->file_list(), $bas_h5_list, 'correct bas h5 file list');
  
  #$pacbio->process();
}

{
  my $pacbio = npg_common::irods::run::PacBio->new( runfolder_path    =>'t/data/pacbio/superfoo/24862_627',
                                                    collection_folder => 'A01_2',
   );


  my $bas_h5_list = [qw(
    t/data/pacbio/superfoo/24862_627/A01_2/Analysis_Results/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.1.bax.h5
    t/data/pacbio/superfoo/24862_627/A01_2/Analysis_Results/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.2.bax.h5
    t/data/pacbio/superfoo/24862_627/A01_2/Analysis_Results/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.3.bax.h5
    t/data/pacbio/superfoo/24862_627/A01_2/Analysis_Results/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.bas.h5
  )];
  
  is_deeply ([sort @{$pacbio->file_list()}], $bas_h5_list, 'correct bas h5 file list');
  
  is ( $pacbio->_get_meta_xml_from_bas_file(q(t/data/pacbio/superfoo/24862_627/A01_2/Analysis_Results/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.1.bax.h5)),
       q(t/data/pacbio/superfoo/24862_627/A01_2/m131209_215015_00127_c100579142550000001823092301191431_s1_p0.metadata.xml),
       'correct metadata xml file path'
  );

  #$pacbio->process();
}
1;
__END__
