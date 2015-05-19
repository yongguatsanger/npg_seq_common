use strict;
use warnings;
use Test::More tests => 32;
use Test::Exception;
use Test::Deep;

use_ok('npg_common::run::file_finder');

{
    my $finder;
    lives_ok { $finder = npg_common::run::file_finder->new(id_run => 22, position => 1); }
      q{create npg_common::run::finder object ok};
    isa_ok($finder, q{npg_common::run::file_finder});
}

{
    throws_ok { npg_common::run::file_finder->new(
      id_run => 22, position => 1, with_t_file => 1, tag_index => 22) }
         qr/tag_index and with_t_file attributes cannot be both set/, 
         'error when attempting to set both with_t_file and tag_index';
    throws_ok {npg_common::run::file_finder->new(position => 12, id_run => 11)} 
           qr/Validation\ failed\ for\ \'NpgTrackingLaneNumber\'/, 
           'error on passing to the constructor invalid int as a position';
    throws_ok {npg_common::run::file_finder->new(position => 'dada', id_run => 11)} 
           qr/Validation\ failed\ for\ \'NpgTrackingLaneNumber\'/, 
           'error on passing to the constructor position as string';
    throws_ok {npg_common::run::file_finder->new(position => 1.2, id_run => 11)} 
           qr/Validation\ failed\ for\ \'NpgTrackingLaneNumber\'/, 
           'error on passing to the constructor position as a float';
}

{
    my $finder = npg_common::run::file_finder->new(id_run => 22, position => 2);
    is ($finder->position, 2, 'position set');
    is ($finder->file_extension, 'fastq', 'default file extension');
    is ($finder->with_t_file, 0, 'no _t file by default');
    is ($finder->tag_index, undef, 'tag_index undefined by default');
    is ($finder->tag_label, q[], 'empty string as a tag label');
}

{
    my $f = npg_common::run::file_finder->new(position  => 1, id_run    => 1234,);

    is ($f->create_filename(q[fastq]), '1234_1.fastq', 'generate filename, no args');

    $f->file_extension(q[]);
    is ($f->create_filename(), '1234_1', 'generate filename, no args, no ext');
   
    $f->file_extension(q[fastqcheck]);
    is ($f->create_filename(q[fastqcheck], 2), '1234_1_2.fastqcheck', 'generate filename, end 2');
    is ($f->create_filename(q[fastqcheck], 1), '1234_1_1.fastqcheck', 'generate filename, end 1');
    is ($f->create_filename(q[fastqcheck], q[t]), '1234_1_t.fastqcheck', 'generate filename, end t');
    is ($f->create_filename(q[fastqcheck]), '1234_1.fastqcheck', 'generate filename, single new-style');

    throws_ok { $f->create_filename(q[fastqcheck], q[22]) } qr/Unrecognised end string 22/, 
            'error for an end that is not 1, 2 or t';
}

{
    my $f = npg_common::run::file_finder->new( position  => 1,
                                               id_run    => 1234,
                                               file_extension => q[fastqcheck],
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               db_lookup    => 0,
                                             );

    my $forward = q[t/data/fuse/mpsafs/runs/1234/1234_1.fastqcheck];
    my $expected = { forward => $forward, };
    cmp_deeply ($f->files(), $expected, 'one fastqcheck input file found');

    $f =    npg_common::run::file_finder->new( position  => 2,
                                               id_run    => 1234,
                                               file_extension => q[fastqcheck],
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               db_lookup    => 0,
                                             );
    $forward = q[t/data/fuse/mpsafs/runs/1234/1234_2_1.fastqcheck];
    my $reverse = q[t/data/fuse/mpsafs/runs/1234/1234_2_2.fastqcheck];
    $expected = { forward => $forward, reverse => $reverse, };
    cmp_deeply ($f->files(), $expected, 'two fastqcheck input files found');

    $f =    npg_common::run::file_finder->new( position  => 1,
                                               id_run    => 1234,
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               db_lookup    => 0,
                                               with_t_file => 1,
                                             );
    $forward = q[t/data/fuse/mpsafs/runs/1234/1234_1_1.fastq];
    $reverse = q[t/data/fuse/mpsafs/runs/1234/1234_1_2.fastq];
    my $t = q[t/data/fuse/mpsafs/runs/1234/1234_1_t.fastq];
    $expected = { forward => $forward, reverse => $reverse, tags => $t, };
    cmp_deeply ($f->files(), $expected, 'three fastq input files found');
}

{
    my $f = npg_common::run::file_finder->new( position  => 1,
                                               id_run    => 2549,
                                               tag_index => 33,
                                               archive_path => q[t/data/fuse/mpsafs/runs/2549],
                                               db_lookup => 0
                                              );
    is ($f->create_filename(q[fastq]),  '2549_1#33.fastq', 'generate filename, no args, tag_index');
    is ($f->create_filename(q[fastq], 1), '2549_1_1#33.fastq', 'generate filename, end 1, tag_index');

    $f->file_extension(q[fastqcheck]);
    is ($f->create_filename(q[fastqcheck], 2), '2549_1_2#33.fastqcheck', 'generate filename, end 2, tag_index');
    is ($f->create_filename(q[fastqcheck]), '2549_1#33.fastqcheck', 'generate filename, no args, tag_index');  
}

{
    my $f = npg_common::run::file_finder->new( position  => 1,
                                               id_run    => 1234,
                                               tag_index => 33,
                                               file_extension => q[fastqcheck],
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               lane_archive_lookup => 0,
                                               db_lookup => 0
                                              );
    my $forward = q[t/data/fuse/mpsafs/runs/1234/1234_1_1#33.fastqcheck];
    my $reverse = q[t/data/fuse/mpsafs/runs/1234/1234_1_2#33.fastqcheck];
    my $expected = { forward => $forward, reverse => $reverse, };
    cmp_deeply ($f->files(), $expected, 'two fastqcheck input files found, tag_index');
}

{
    my $f = npg_common::run::file_finder->new( position  => 1,
                                               id_run    => 1234,
                                               tag_index => 33,
                                               lane_archive_lookup => 1,
                                               file_extension => q[fastq],
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               db_lookup => 0);
    is(scalar keys %{$f->files}, 0, 'no files if they are not in the lane archive');
    
    $f = npg_common::run::file_finder->new(    position  => 3,
                                               id_run    => 1234,
                                               tag_index => 33,
                                               lane_archive_lookup => 1,
                                               file_extension => q[fastq],
                                               archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                               db_lookup => 0
                                          );
    my $forward = q[t/data/fuse/mpsafs/runs/1234/lane3/1234_3_1#33.fastq];
    my $reverse = q[t/data/fuse/mpsafs/runs/1234/lane3/1234_3_2#33.fastq];
    my $expected = { forward => $forward, reverse => $reverse, };
    cmp_deeply ($f->files(), $expected, 'two fastqcheck input files found, tag_index, lane archive');
}

{
    my $f     = npg_common::run::file_finder->new(    position  => 3,
                                                      id_run    => 1234,
                                                      file_extension => q[fastqcheck],
                                                      archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                                      db_lookup => 0);
    my $forward = q[t/data/fuse/mpsafs/runs/1234/1234_3_1.fastqcheck];
    my $expected = { forward => $forward};
    cmp_deeply ($f->files(), $expected,  'one fastqcheck input files found; with _1 to identify the end');
   
    $f     = npg_common::run::file_finder->new( position  => 4,
                                                id_run    => 1234,
                                                archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                                db_lookup => 0);
    $forward = q[t/data/fuse/mpsafs/runs/1234/1234_4.fastq];
    $expected = { forward => $forward};
    cmp_deeply ($f->files(), $expected,  'one fastqcheck input files found; no _1 to identify the end');   

    $f     = npg_common::run::file_finder->new(  position  => 8,
                                                 id_run    => 1234,
                                                 archive_path => q[t/data/fuse/mpsafs/runs/1234],
                                                 db_lookup => 0);
    is (scalar (keys %{$f->files}), 0, 'no input files found');
}

1;
