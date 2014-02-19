#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2009-08-27
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use English qw{-no_match_vars};
use Test::More tests => 107;
use Test::Exception;
use Test::Deep;
use Cwd qw(abs_path);
use File::Temp qw/ tempdir  /;

local $ENV{'NPG_WEBSERVICE_CACHE_DIR'} = q[t/data/sequence/];

use_ok('npg_common::irods::Loader');
use_ok('npg_common::irods::run::Bam');

### it's important that the user have ownership permission on /seq/npg/test1 so that ownership rights are
### propagated to the data added in these tests to allow its deletion.
my $dir = tempdir( CLEANUP => 1 );
my @comp = split '/', $dir;
my $dname = pop @comp;
my $IRODS_TEST_AREA1 = "/seq/npg/test1/$dname";
my $have_irods_execs = exist_irods_executables();
my $test_area_created = $have_irods_execs ? create_irods_test_area() : 0;

{
  my $bam = npg_common::irods::run::Bam->new(   id_run => 5174,
       archive_path => q[t/data/sequence/100818_IL32_05174/Latest_Summary/archive],
  );

  is( $bam-> _bam_file_name({id_run=>5174, position=>1}), '5174_1.bam', 'lane 1 bam file name');
  is( $bam-> _bam_file_name({id_run=>5174, position=>1, tag_index=>0}), '5174_1#0.bam', 'lane 1 tag_index 0 target bam file name');
 
  is( scalar @{$bam->file_list()}, 6, "correct number of bam files");

  ok ( ! $bam->alignment_filter_in_plex(), 'alignment_filter_in_plex not true');
  
  ok( !$bam->nonhuman_not_in_filename(), 'file name with nonhuman');

  my $names_1 = $bam->get_study_library_sample_names(5174, 1);
  
  ok ( !$names_1->{is_phix_control}, 'run 5174 lane 1 is not phix control data');

  my $study_names_1 = $names_1->{study};
  is( scalar @{$study_names_1}, 1, 'only 1 study name for run 5174 lane 1' );
  is( $study_names_1->[0], 'ZF_MrSol4000', 'study name for run 5174 lane 1');  
  is($names_1->{study_title}->[0], "Zebrafish Tilling by Illumina Sequencing", 'correct study title'); 
  is($names_1->{study_accession_number}->[0], "ERP000110", 'correct study accession nubmer');  
  is($names_1->{study_id}->[0], 307, 'correct study id');

  is($names_1->{sample_id}->[0], 33318, 'correct sample id');
  is($names_1->{sample_common_name}->[0], 'Danio rerio', 'correct sample common name');
  my $sample_names_1 = $names_1->{sample};
  is( scalar @{$sample_names_1}, 2, '2 sample names for run 5174 lane 1' );
  is( $sample_names_1->[0], '83_3_1-11', 'first sample name for run 5174 lane 1');

  is($names_1->{is_yhuman_split}, '', 'separate y chromosome data NOT defined for Zebrafish project');
  is($names_1->{human_split_type}, undef, 'human split is undef for Zebrafish project');
  

  is($names_1->{sample_public_name}->[0], '83_3_1-11', 'first sample public name');
  is( scalar @{$names_1->{sample_public_name}}, 2, '2 sample public name for run 5174 lane 1' );

  is( scalar @{$names_1->{sample_accession_number}}, 2, '2 sample public accession number for run 5174 lane 1' );
  
  is($names_1->{sample_accession_number}->[0], 'ERS009168', 'first sample accession number');
  is($names_1->{sample_accession_number}->[1], 'ERS009169', 'second sample accession number');

  is( scalar @{ $names_1->{library} }, 1, 'one library name for run 5174 lane 1');
  is( $names_1->{library}->[0], 'MrSol4000Pool18', 'library name for run 5174 lane 1' );
  
  is( scalar @{ $names_1->{library_id} }, 1, 'one library id for run 5174 lane 1');
  is( $names_1->{library_id}->[0], '151127', 'library name for run 5174 lane 1' );

  my $names_4 = $bam->get_study_library_sample_names(5174, 4);
  
  ok ( $names_4->{is_phix_control}, 'run 5174 lane 4 is phix control data');

  is( scalar @{ $names_4->{study} }, 0, 'no study name for run 5174 lane 4' );
  is($names_4->{library}->[0], 'Multiplex-phix', 'library name for run 5174 lane 4');
  is($names_4->{sample}->[0], 'Multiplex-phix', 'sample name for run 5174 lane 4');
  
  my $names_1_7 = $bam->get_study_library_sample_names(5174, 1, 7);

  ok ( !$names_1_7->{is_phix_control}, 'run 5174 lane 1 plex 7 is not phix control data');

  my $study_names_1_7 = $names_1_7->{study};
  is( scalar @{$study_names_1_7}, 1, 'only 1 study name for run 5174 lane 1 plex 7' );
  is( $study_names_1_7->[0], 'ZF_MrSol4000', 'study name for run 5174 lane 1 plex 7');

  my $sample_names_1_7 = $names_1_7->{sample};
  is( scalar @{$sample_names_1_7}, 1, 'only 1 sample name for run 5174 lane 1 plex 7' );
  is( $sample_names_1_7->[0], '83_3_1-11', 'sample name for run 5174 lane 1 plex 7');
  
  is( $names_1_7->{library}->[0], '83_3_1-11 151125', 'library name for run 5174 lane 1 plex 7');

  my $bam_list = {q{5310_1#0.bam} => 1, q{5310_1#1.bam} => 1};
  my $bai_list = {q{5310_1#0.bai} => 1,};
  $bam->number_of_reads_list_irods()->{q{5310_1#0.bam}} = 1000;
  $bam->number_of_reads_list_irods()->{q{5310_1#1.bam}} = 2000;
  $bam->_bam_indexed_cache_irods()->{q{5310_1#0.bam}} = 1;
  $bam->_bam_indexed_cache_irods()->{q{5310_1#1.bam}} = 1;
  ok (!$bam->_check_bam_index_file($bam_list, $bai_list), 'bam index file for run 5130 lane 1 plex 1 missing');
  
  $bam->_bam_indexed_cache_irods()->{q{5310_1#3.bam}} = 1;
  $bai_list->{'5310_1#3.bai'} = 1;
  ok (!$bam->_check_bam_index_file($bam_list, $bai_list), 'bam file for run 5310 lane 1 plex 3 missing');

  $bai_list->{'5310_1#1.bai'} = 1;
  $bam_list->{'5310_1#3.bam'} = 1;
  $bam->number_of_reads_list_irods()->{q{5310_1#3.bam}} = 3000;
  ok ($bam->_check_bam_index_file($bam_list, $bai_list), 'bam index file match');

  SKIP : {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.', 4
      unless ($ENV{TOOLS_INSTALLED});

  # no reads case

  $bam_list->{'5310_1#4.bam'} = 1;
  $bam->number_of_reads_list_irods()->{q{5310_1#4.bam}} = 0;
  $bam->number_of_reads_list()->{q{5310_1#4.bam}} = 0;
  is($bam->get_number_of_reads(q{5310_1#4}), 0, 'No reads for run 5310 lane 4');
  ok ($bam->_check_bam_index_file($bam_list, $bai_list), 'bam index file still matched with one empty file');

  # unaligned case

  my $unaligned_bam_file = q{5310_1#5};
  my $unaligned_reads = 2000;
  $bam_list->{$unaligned_bam_file} = 1;
  $bam->number_of_reads_list_irods()->{$unaligned_bam_file} = $unaligned_reads;
  $bam->number_of_reads_list()->{$unaligned_bam_file} = $unaligned_reads;
  is($bam->get_number_of_reads($unaligned_bam_file), $unaligned_reads, "$unaligned_reads reads for run 5310 lane 5");
  ok ($bam->_check_bam_index_file($bam_list, $bai_list), 'bam index file still matched with one empty file and one unaligned file');
  } # end skip TOOLS_INSTALLED

  my $ref1 = '/lustre/scratch101/sanger/references/Homo_sapiens/CGP_GRCh37.NCBI.allchr_MT/all/bwa/Homo_sapiens.GRCh37.NCBI.allchr_MT.fa';
  my $bwa_pg_line1 = q{@}.qq{PG	ID:bwa_aln	PN:bwa aln	CL:/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa aln -t 2 $ref1 /nfs/sf20/ILorHSany_sf20/analysis/101101_HS3_05450_A_205NTABXX/Data/Intensities/PB_basecalls_20101108-105238/no_cal/archive/5450_3_1.fastq > /tmp/MUYklDpWoh/1.sai;/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa aln  -t 2 $ref1 /nfs/sf20/ILorHSany_sf20/analysis/101101_HS3_05450_A_205NTABXX/Data/Intensities/PB_basecalls_20101108-105238/no_cal/archive/5450_3_2.fastq > /tmp/MUYklDpWoh/2.sai;};

  is($bam->_get_ref_from_bwa_pg($bwa_pg_line1), $ref1, 'correct reference');

  my $ref2 = q{/nfs/repository/d0031/references/Salmonella_enterica/Paratyphi_A_AKU_12601/all/bwa/S_enterica_paratyphi_AKU_12601.fasta};
  my $bwa_pg_line2 = q{@}.qq{PG	ID:bwa	PN:bwa	PP:split_fastq_by_tag	VN:0.5.8c (r1536)	CL:-q 15 -t 2 $ref2};
  is($bam->_get_ref_from_bwa_pg($bwa_pg_line2), $ref2, 'correct reference');

  my $ref3 = '/nfs/srpipe_references/references/Homo_sapiens/GRCh37_53/all/bwa/Homo_sapiens.GRCh37.dna.all.fa';
  my $bwa_pg_line3 = '@PG	ID:bwa_aln	PN:bwa aln	PP:split_fastq_by_tag	VN:0.5.8c (r1536)	DS:bwa alignment from fastq files	CL:/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa aln -q 15 -t 2 '. $ref3 . ' /nfs/sf19/ILorHSany_sf19/analysis/110311_HS2_05947_A_B02YHABXX/Data/Intensities/PB_basecalls_20110318-161034/no_cal/archive/lane8/5947_8_1#2.fastq > 1.sai';
  is($bam->_get_ref_from_bwa_pg($bwa_pg_line3), $ref3, 'correct reference');

  
  my $samtools;
  eval {
    $samtools = $bam->samtools_cmd();
  };

  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
      1 unless ($ENV{TOOLS_INSTALLED});
    my $hash_ref = $bam->get_bam_reference('t/data/sequence/5431_6#3.bam');
    is( $hash_ref->{reference}, $ref2, 'correct reference from bam file');
  }
  
  $bam->alignment_filter_in_plex(1);
  is( $bam-> _bam_file_name({id_run=>5174, position=>1, tag_index=>0}), '5174_1#0.bam', 'lane 1 tag_index 0 target bam file name');
}

{
  my $bam = npg_common::irods::run::Bam->new(  id_run => 5174, 
     archive_path => q[t/data/sequence/100818_IL32_05174/Latest_Summary/archive], 
     positions => [1],     
  );
  is( scalar @{$bam->file_list()}, 2, "correct number of bam files for lane 1");

  $bam->clear_file_list(); 
  $bam->positions([1,4]);  
  is( scalar @{$bam->file_list()}, 3, "correct number of bam files for lane 1 and 4");
  
  $bam->clear_file_list();
  $bam->positions([]);
  is( scalar @{$bam->file_list()}, 6, "correct number of bam files for all lanes");
  
  $bam->clear_file_list();
  $bam->positions([1, 2, 4]);
  is( scalar @{$bam->file_list()}, 6, "correct number of bam files for lane 1 2 and 4");
}

{
  my $bam = npg_common::irods::run::Bam->new( id_run => 5174,
     archive_path => q[t/data/sequence/100818_IL32_05174/Latest_Summary/archive],
     positions => [1],
  );
  my $bam_file = $bam->file_list()->[0];
  is($bam->get_number_of_reads_qc($bam_file), 2963516, 'correct total number reads for 5174 lane 1 from qc json');

  is($bam->get_number_of_reads($bam_file), 2963516, 'correct total number reads for 5174 lane 1');

  is($bam->number_of_reads_list()->{'5174_1'},  2963516, 'correct total number reads for 5174 lane 1 from cache');

  SKIP : {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.', 3 unless ($ENV{TOOLS_INSTALLED});
    is($bam->_bam_no_alignment($bam_file), 1, 'bam file 5174_1 without alignment as no SQ line');
  
    my $irods_bam_list = {'/seq/7301/7301_8#8.bam' => 1, '/seq/7301/7301_8#8_phix.bam' => 1};
    my $staging_bam_list = ['/staging/7301_8#8.bam'];
    ok( !$bam->check_bam_list_md5($irods_bam_list, $staging_bam_list), 'number of files are different on irods and staging');
    
    my $bam = npg_common::irods::run::Bam->new(   id_run => 5174,
       archive_path => q[t/data/sequence/100818_IL32_05174/Latest_Summary/archive],
       samtools_cmd => q[t/bin/aligners/samtools/current/samtools],
     );
     my $bam_file = $bam->file_list()->[0];

     my $bam_mark_duplicate = $bam->_bam_to_mark_duplicate($bam_file);
     my $samtools_cmd = abs_path( q[t/bin/aligners/samtools/current/samtools] );
     is ($bam_mark_duplicate->samtools_cmd(), $samtools_cmd,
      q[samtools not on the path and propagated correctly to BAM_MarkDuplicate]);
  }

  is($bam->alt_process(),'','process is null');
  is($bam->collection(),'/seq/5174','collection is default');
}

{
  my $bam = npg_common::irods::run::Bam->new( id_run => 5174, alt_process=>'raspberry');
  is($bam->alt_process(),'raspberry','process is set');
  is($bam->collection(),'/seq/5174/raspberry','collection is set');
}

{
  my $bam = npg_common::irods::run::Bam->new( id_run => 5174, collection=>'/seq/js10');
  is($bam->alt_process(),'','process is clear');
  is($bam->collection(),'/seq/js10','collection is set');
}

{
  my $bam = npg_common::irods::run::Bam->new( id_run => 10371,
     archive_path => q[t/data/sequence/100818_IL32_10371/Latest_Summary/archive],
     positions => [1],
  );

  my $names = $bam->get_study_library_sample_names(10371, 1);
  
  my $study_names = $names->{study};
  is( scalar @{$study_names}, 1, 'only 1 study name for run 10371 lane 1' );
  is( $study_names->[0], 'SEQCAP_WGS_Cardiometabolic_and_infectious_traits_in_sub_Saharan_African_populations', 'study name for run 10371 lane 1');  
  is($names->{study_title}->[0], 'Whole genome association study of cardiometabolic and infectious traits in sub-Saharan African populations', 'correct study title'); 
  is($names->{study_id}->[0], 2617, 'correct study id');
  is($names->{is_xahuman_split}, 0, 'xahuman is defined at lane level');
  is($names->{is_yhuman_split}, undef, 'separate y chromosome data NOT defined at lane level');
  is($names->{human_split_type}, undef, 'human split is undefined at lane level');

}

{
   package myBam;
   use Moose;
   extends 'npg_common::irods::run::Bam';
   sub tag_list {
     return {36 => 'AAAGTCT',};
   }
   sub _recipe_store {
   }
   no Moose;
 
   package main;
   
   my $bam = myBam->new( id_run => 10371,
    archive_path => q[t/data/sequence/100818_IL32_10371/Latest_Summary/archive/lane1],
    runfolder_path => q[t/data/sequence/100818_IL32_10371],
    recalibrated_path => q[t/data/sequence/100818_IL32_10371],
    positions => [1],
    tag_index => 36,
    is_indexed => 1,
    collection => $IRODS_TEST_AREA1,
  );

  my $bam_filename1 = '10371_1#36.bam';
  my $bam_file1 = "t/data/sequence/100818_IL32_10371/Latest_Summary/archive/lane1/$bam_filename1";

  my $yhuman_bam_filename2 = '10371_1#36_yhuman.bam';
  my $yhuman_bam_file2 = "t/data/sequence/100818_IL32_10371/Latest_Summary/archive/lane1/$yhuman_bam_filename2";

  is($bam->file_list()->[1], $bam_file1, 'correct main bam file for lane 1');
  is($bam->file_list()->[2], $yhuman_bam_file2, 'correct file for lane 1 yhuman');
  is($bam->collection(), $IRODS_TEST_AREA1, 'collection is set');

  my $names = $bam->get_study_library_sample_names(10371, 1, 36);
  
  my $study_names = $names->{study};
  is( scalar @{$study_names}, 1, 'only 1 study name for run 10371 lane 1 tag 36' );
  is( $study_names->[0], 'SEQCAP_WGS_Cardiometabolic_and_infectious_traits_in_sub_Saharan_African_populations', 'study name for run 10371 lane 1');  
  is($names->{study_title}->[0], 'Whole genome association study of cardiometabolic and infectious traits in sub-Saharan African populations', 'correct study title'); 
  is($names->{study_id}->[0], 2617, 'correct study id');
  is($names->{is_xahuman_split}, 0, 'xahuman is NOT set for tag 36');
  is($names->{is_yhuman_split}, 1, 'separate y chromosome data set for tag 36');
  is($names->{human_split_type}, 'yhuman', 'split is yhuman for tag 36');

  SKIP: {
  unless ($have_irods_execs) {
    skip 'unable to access iRODS executables', 24;
  } elsif (!$test_area_created) {
    skip 'unable to create iRODS test area (try kinit to log in to iRODS)', 24;
  } elsif (!defined $ENV{TOOLS_INSTALLED} || (defined $ENV{TOOLS_INSTALLED} && $ENV{TOOLS_INSTALLED} == 0)) {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.', 24;
  }

  lives_ok { $bam->process() } 'Loading files to irods succeeds';

  my $file1 = "$IRODS_TEST_AREA1/$bam_filename1";
  my $file2 = "$IRODS_TEST_AREA1/$yhuman_bam_filename2";

  my $file1_ils = `ils -A $file1`;
  my $file2_ils = `ils -A $file2`;

  like($file2_ils, qr/ss_2617#seq:read object/msx, 'yhuman bam file is readable by study based group' );
  unlike($file2_ils, qr/public#seq:read object/msx, 'and yhuman bam file is not public readable' );
  my $meta = npg_common::irods::Loader->new(
                  file => $file2, 
                  collection => $IRODS_TEST_AREA1,
             )->_check_meta_data($file2);
  ok(exists $meta->{'target'}->{'1'}, 'yhuman bam file target meta attr is 1');
  ok(exists $meta->{'alignment_filter'}->{'yhuman'}, 'yhuman bam file alignment_filter meta attr is "yhuman"');
  ok(exists $meta->{alignment}, 'yhuman bam file has aligned value in irods meta data');
  ok(exists $meta->{alignment}->{1}, 'yhuman bam file is aligned in irods meta data');

  like($file1_ils, qr/ss_2617#seq:read object/msx, 'main bam file is readable by study based group' );
  unlike($file1_ils, qr/public#seq:read object/msx, 'and main bam file is not public readable' );
  $meta = npg_common::irods::Loader->new(
                  file => $file1, 
                  collection => $IRODS_TEST_AREA1,
             )->_check_meta_data($file1);
  ok(exists $meta->{'target'}->{'1'}, 'main bam file target meta attr is 1');
  ok(!exists $meta->{'alignment_filter'}->{'yhuman'}, 'alignment_filter meta attr is not set for a main bam file');
  ok(exists $meta->{alignment}, 'main bam file has aligned value in irods meta data');
  ok(exists $meta->{alignment}->{1}, 'main bam file is aligned in irods meta data');

  is($bam->bam_fully_archived(), 1, 'bam_fully_archived reports bam and bai in staging match the files in irods');
  is($bam->bam_fully_archived_for_deletion(), 1, 'bam_fully_archived_for_deletion reports bam and bai in staging match the files in irods');

  my $expected_file1_reads = 5;
  is($bam->get_number_of_reads_irods_meta($bam_filename1), $expected_file1_reads, "$bam_filename1 has $expected_file1_reads reads in irods metadata");
  is($bam->get_number_of_reads($bam_filename1),  $expected_file1_reads, "$bam_filename1 has $expected_file1_reads reads in cache");

  my $expected_file2_reads = 5;
  is($bam->get_number_of_reads_irods_meta($yhuman_bam_filename2), $expected_file2_reads, "$yhuman_bam_filename2 has $expected_file2_reads reads in irods metadata");
  is($bam->get_number_of_reads($yhuman_bam_filename2), $expected_file2_reads, "$yhuman_bam_filename2 has $expected_file2_reads reads in cache");

  my $index_file1 = $file1;
  $index_file1 =~ s/bam$/bai/mxs;

  my $index_file1_rm = system(qq{irm $index_file1});

  is($index_file1_rm, 0, "Successfully removed $index_file1 from irods");
  is($bam->bam_fully_archived(), 0, 'bam_fully_archived reports bam and bai in staging DO NOT match the files in irods');
  is($bam->bam_fully_archived_for_deletion(), 0, 'bam_fully_archived_for_deletion reports bam and bai in staging DO NOT match the files in irods');

  $meta = npg_common::irods::Loader->new(
                  file => $file1, 
                  collection => $IRODS_TEST_AREA1,
             )->_check_meta_data($file1);
  ok(exists $meta->{alignment}, 'main bam file has aligned value in irods meta data');
  ok(exists $meta->{alignment}->{1}, 'main bam file is aligned in irods meta data');

  } # end of skip
} 

sub exist_irods_executables {
  return 0 unless `which ienv`;
  return 0 unless `which imkdir`;
  return 1;
}
sub create_irods_test_area {
  system("imkdir $IRODS_TEST_AREA1") == 0 or return 0;
  return 1;
}

END {
  if($test_area_created) {
    eval {system("irm -r $IRODS_TEST_AREA1")};
  }
}

1;
__END__
