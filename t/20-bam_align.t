#########
# Author:        John O'Brien
# Created:       12/01/2011
#
package test_bam_align;

use strict;
use warnings;
use Test::More tests => 127;
use Test::Deep;
use Test::Exception;
use Test::Trap;
use Test::Warn;
use File::Slurp;
use English qw{-no_match_vars};
use Cwd qw(getcwd abs_path);

sub wipe_test_scratch {
    my ($dir) = @_;
    foreach ( glob qq{$dir/*} ) { unlink $_; }
    write_file("$dir/This_dir_is_wiped_during_tests");
    return;
}

use_ok('npg_common::bam_align');


my $SOLEXA_BIN = qw[t/bin];
my $JAVA_CMD   = qw[java];

my $scratch = q{t/data/scratch};    # This is used for tmp as well.

{
    local $ENV{CLASSPATH} = 't/bin/aligners/picard/current';
    local $ENV{PATH} = join q[:], 't/bin', $ENV{PATH}; 

    my $base_test;
    lives_ok { $base_test = npg_common::bam_align->new(
                                 repository => q[t],
                                 out_bam => 'out.bam',
                                 bwa_cmd => q[t/bin/bwa],
                               );
             } 'Create object';
    isa_ok( $base_test, 'npg_common::bam_align' );

    my $wtsi_collection = abs_path qq[$SOLEXA_BIN/aligners];

    is ( $base_test->bwa_cmd(), join(q[/], getcwd(), q[t/bin/bwa] ),
            'abs path to bwa'
    );
    like(
        $base_test->bam_markduplicates(),
        qr{$wtsi_collection/picard/picard[-]tools[-].+/MarkDuplicates[.]jar}mxs,
        'Return the default MarkDuplicates'
    );

    my $default_scratch = $base_test->scratch();

    like(
	$base_test->read_1_fastq(),
	qr{^ $default_scratch \S+ read_1[.]fastq $}msx,
	'Return path to first read fastq'
    );
    like(
	$base_test->read_2_fastq(),
	qr{^ $default_scratch \S+ read_2[.]fastq $}msx,
	'Return path to second read fastq'
    );
    like(
	$base_test->single_fastq(),
	qr{^ $default_scratch \S+ single[.]fastq $}msx,
	'Return path to single read fastq'
    );
    like(
	$base_test->old_header_file(),
	qr{^ $default_scratch \S+ header[.]temp $}msx,
	'Return path to header file'
    );
    like(
	$base_test->pg_records(),
	qr{^ $default_scratch \S+ pg_records[.]temp $}msx,
	'Return path to header file'
    );
    like(
	$base_test->bam_tags_file(),
	qr{^ $default_scratch \S+ bam_tags[.]temp $}msx,
	'Return path to BAM tags file'
    );
    like(
	$base_test->new_sam_file(),
	qr{^ $default_scratch \S+ new_alignment[.]sam $}msx,
	'Return path to new SAM file'
    );

    is( $base_test->debug(), 0, 'Debug is off by default' );

    is(
            $base_test->sortsam_command(),
            undef,
            'Attribute sortsam_command is undef'
    );
    is(
            $base_test->samview_command(),
            undef,
            'Attribute samview_command is undef'
    );

    my $x = $base_test->input_pipe();
    ok( $base_test->sortsam_command(), 'But now sortsam_command is set' );
    ok( $base_test->samview_command(), 'And samview_command is set' );
}

my $idx_base = 't/data/references/Homo_sapiens/NCBI36/all/bwa/ncbi36.fa';

my $init = {
    species         => 'Human',
    index_base      => $idx_base,
    scratch         => $scratch,
    temp_dir        => $scratch,
    bwa_cmd         => '/bin/true',
    java_cmd        => '/bin/true',
    bam_markduplicates  => '/bin/true',
    ref_repository  => 't/data/references',
    debug           => 1,
    out_bam         => $scratch.q{/out.bam},
    repository      => q[t/data],
};


# Test the subroutines first.
{
    my $first = q{@} . qq{PG\tID:first\tblah\n};
    my $final = q{@} . qq{PG\tID:lAst\tBlah\n};

    throws_ok { npg_common::bam_align::validate_pg_arg( {} ); }
              qr/Arrayref arguments only/ms,
              'Croak on incorrect input to validate_pg_arg';

    my @lines = npg_common::bam_align::validate_pg_arg( [] );
    is( scalar @lines, 0, 'Empty array gets an empty array' );

    @lines = npg_common::bam_align::validate_pg_arg( [ $first, $final ] );
    cmp_deeply(
        \@lines,
        [ $first, $final ],
        'Parse array-ref arg correctly'
    );

    throws_ok { npg_common::bam_align::make_unique_pg_id(); }
              qr/Base id argument required/ms,
              'Insist on base id argument';

    my $id = npg_common::bam_align::make_unique_pg_id(
                [ $first, $final ], 'last'
    );
    is( $id, 'last_2', 'Version is one more than the count found' );

    my $middle = q{@} . qq{PG\tID:lAst_010\tBlah\n};

    $id = npg_common::bam_align::make_unique_pg_id(
                [ $first, $middle, $final ], 'last'
    );
    is( $id, 'last_11', 'Version is one more than the previous highest' );

    $id = npg_common::bam_align::make_unique_pg_id(
                [ $first, $first ], 'last'
    );
    is( $id, 'last', 'No version of the base is not found' );


    is( npg_common::bam_align::prev_command_field( [ $first, $final ] ),
        qq{\tPP:lAst}, 'Return correct value' );

    is( npg_common::bam_align::prev_command_field( [] ),
        q{}, 'Return a null string for empty input' );

    $final =~ s/ID:/broken/msx;
    throws_ok
        { npg_common::bam_align::prev_command_field( [ $first, $final ] ); }
        qr/Error inferring PP field from $final/ms,
        'Croak if the last ID is not found in the argument';

    throws_ok { npg_common::bam_align::reverse_complement('x'); }
              qr/x doesn't look like a nucleic acid sequence/ms,
              'Croak unless the string is a nucleic acid sequence';

    is(
        npg_common::bam_align::reverse_complement(q{AaAcCuC.GNGgTtTU}),
        q{AAaAcCNC.GaGgTtT},
        'Reverse complement of sequence is correct'
    );


    my @r = trap { npg_common::bam_align::timer_log('BlahBlahBlah') };
    like( $trap->stderr(), qr/[:][ ]BlahBlahBlah$/msx, 'Logger writes to STDERR' );
}


# Test processing of the BAM input.
{
    my $test = npg_common::bam_align->new($init);

    lives_ok { $test->_input_pipe('cat t/data/good.sam'); } 'Mock STDIN';
    is( $test->scratch(), $scratch, 'Test scratch hasn\'t changed' );

    ( -e $test->old_header_file() ) && ( unlink $test->old_header_file() );
    ( -e $test->read_1_fastq() )    && ( unlink $test->read_1_fastq() );
    ( -e $test->read_2_fastq() )    && ( unlink $test->read_2_fastq() );

    lives_ok { $test->read_stdin(); } 'Read BAM method runs without error';

    my $header_contents = read_file( $test->old_header_file() );
    is(
        $header_contents,
        q{@} . qq{CO\tBlah, blah, blah.\n},
        'Header saved as expected'
    );
    is( $test->species(), 'Human', 'Organism ignored in header' );

    # Test the contents when we're testing the method that writes them.
    ok( -e $test->read_1_fastq(), 'First read fastq created' );
    ok( -e $test->read_2_fastq(), 'Second read fastq created' );
}


{
    my $test = npg_common::bam_align->new($init);

    $test->clear_species();
    $test->_input_pipe('cat t/data/good.sam');
    $test->read_stdin();

    is( $test->species(), 'Caterpillar', 'Organism extracted from header' );
}


{
    my $test = npg_common::bam_align->new($init);

    ( -e $test->single_fastq() ) && ( unlink $test->single_fastq() );
    $test->_input_pipe('cat t/data/mixed.sam');
    $test->clear_is_paired();
    throws_ok
        { $test->read_stdin(); }
        qr{
            Mixed[ ]single/paired[ ]BAM[ ]file[ ]
            [(]HS1_5634:2:68:21066:200418[#]6/2[)]
        }msx,
        'Croak if BAM has both single- and paired-end reads';

    ok( -e $test->single_fastq(), 'Single read fastq created' );

    wipe_test_scratch( $test->scratch() );
}

{
    my $test = npg_common::bam_align->new($init);
    throws_ok { $test->parse_sam(); }
              qr/Pipe handle not supplied/ms,
              'Insist on pipe handle argument';

    throws_ok { $test->parse_sam( 'token' ); }
              qr/Header filehandle not supplied/ms,
              'Insist on header filehandle argument';
}

# Test parsing of alignment lines.
{
    my $test = npg_common::bam_align->new($init);
    throws_ok { $test->parse_alignment(); }
              qr/Alignment line not supplied/ms,
              'Croak if no argument passed to alignment parsing method';

    my $test_line = qq{zero\t1\ttwo\tthree\tfour\tfive\tsix\tseven\teight\t}
                  . qq{CCgtATcaG\tten\televen\ttwelve\tthirteen\n};

    is(
        scalar $test->parse_alignment($test_line),
        12,
        'Return all the fields in the line' );

    $test_line =~ s/1/16/msx;
    my @test_fields = $test->parse_alignment($test_line);
    is( $test_fields[9],  'CtgATacGG', 'Reverse and complement sequence' );
    is( $test_fields[10], 'net',       'Reverse quality' );
    is( $test_fields[11], qq{eleven\ttwelve\tthirteen}, 'Join all the tags' )
}


# Test addition of read pair suffix.
{
    my $test = npg_common::bam_align->new($init);
    throws_ok { $test->add_paired_suffix(); }
                qr/Name argument required/ms,
                'Insist on query name as first argument to suffix method';

    throws_ok { $test->add_paired_suffix('qname'); }
                qr/Flag argument required/ms,
                'Insist on flag as second argument';

    warning_is { $test->add_paired_suffix( 'qname/2', 64 ); }
               q{Contradiction in '/2' tag and read 1 flag: qname/2},
               'Warn if read has /2 but flag says /1';

    warning_is { $test->add_paired_suffix( 'qname/1', 128 ); }
               q{Contradiction in '/1' tag and read 2 flag: qname/1},
               'Warn if read has /1 but flag says /2';

    is(
        $test->add_paired_suffix( 'qname/1', 64 ),
        'qname/1',
        'No change where /1 is already present'
    );

    is(
        $test->add_paired_suffix( 'qname/2', 128 ),
        'qname/2',
        'No change where /2 is already present'
    );

    is(
        $test->add_paired_suffix( 'qname', 64 ),
        'qname/1',
        'Add /1 when required'
    );

    is(
        $test->add_paired_suffix( 'qname', 128 ),
        'qname/2',
        'Add /2 when required'
    );

    is(
        $test->add_paired_suffix( 'qname', 1 ),
        'qname',
        'No change to qname if we don\'t know the pair-end'
    );

}


# Test tag string construction.
{
    my $test = npg_common::bam_align->new($init);
    is(
        $test->make_tag_string(),
        q{},
        'Tag string method returns an empty string with no argument'
    );

    throws_ok { $test->make_tag_string( [ 'XM:i:1', 'SQ:H:0x0220' ] ); }
              qr/This tool does not support the deprecated 'SQ' tag/ms,
              'Croak on the SQ tag';

    my @unwanted_tags = qw(XC AS NH IH CS CQ CM CC CP MD HI SM AM MQ PQ UQ);
    @unwanted_tags = map { $_ . q[:test] } @unwanted_tags;
    is(
        $test->make_tag_string( \@unwanted_tags ),
        q{},
        'Strip alignment-specific tags'
    );

    is(
        $test->make_tag_string( [ 'ZZ:one', 'GG:two', 'CC:dd', 'PL:four' ] ),
        qq{GG:two\tPL:four},
        'Output as expected'
    );
}


# Test read pair consistency check.
{
    my $test = npg_common::bam_align->new($init);
    throws_ok { $test->consistent_read_pairs(); }
              qr/First fields arrayref not supplied/ms,
              'First argument to read pair consistency check missing';

    throws_ok { $test->consistent_read_pairs( [ 'r1', 64 ] ); }
              qr/Second fields arrayref not supplied/ms,
              'Second argument to read pair consistency check missing';

    throws_ok { $test->consistent_read_pairs('not_an_arrayref'); }
              qr/First fields arrayref not supplied/ms,
              'Arguments must be array refs';

    throws_ok { $test->consistent_read_pairs( ['r'], 'not_an_arrayref' ); }
              qr/Second fields arrayref not supplied/ms,
              'Both arguments must be array refs';

    throws_ok { $test->consistent_read_pairs( [ 'r1', 64 ], [ 'r5', 128 ] ); }
              qr/Mismatch between paired read names: r1, r5/ms,
              'Croak on mismatched read names';

    lives_ok {
                $test->consistent_read_pairs(
                    [ 'name', 64 ], [ 'name', 128 ]
                );
             }
             'Successful consistency check';

    my $first_read_first;
    lives_ok
        {
            $first_read_first = $test->consistent_read_pairs(
                [ 'name/1', 64 ], [ 'name/1', 128 ]
            );
        }
        'Trailing /1, /2 markers don\'t interfere';
    is($first_read_first, 1, 'first read first');

    my $second_read_first;
    lives_ok
        {
            $second_read_first = $test->consistent_read_pairs(
                [ 'name/1', 128 ], [ 'name/1', 64 ]
            );
        }
        'Trailing /1, /2 markers don\'t interfere';
    is($second_read_first, 0, 'second read first');

    throws_ok { $test->consistent_read_pairs( [ 'r', 64 ], [ 'r', 127 ] ); }
              qr/Error in paired-read flags: r 64,127/ms,
              'Croak on inconsistent pair flags';

    throws_ok { $test->consistent_read_pairs( [ 'r', 63 ], [ 'r', 128 ] ); }
              qr/Error in paired-read flags: r 63,128/ms,
              'Either can be wrong';
}


# Test fastq output.
{
    my $test = npg_common::bam_align->new($init);
    my $fastq = $test->single_fastq();

    my $field =
        [
            'read_name', 20, 'ref_name', 2, '33141336', q{*}, q{=},
            '33141336', 0, 'AGTTGTGTC', 'AB>BBABCB',
            qq{RG:Z:1\tQT:Z:HIIIIIII\tRT:Z:ATCACGTT}
        ];

    throws_ok { $test->print_read_data(); }
              qr/Fields arrayref not supplied/ms,
              'Insist on argument to print read method';

    throws_ok { $test->print_read_data($fastq); }
              qr/Fields arrayref not supplied/ms,
              'Insist on an array ref argument';

    throws_ok { $test->print_read_data($field); }
              qr/Fastq filehandle not supplied/ms,
              'Insist on a fastq filehandle argument';

    ( -e $fastq ) && ( unlink $fastq );
    ( -e $test->bam_tags_file() ) && ( unlink $test->bam_tags_file() );
    lives_ok { $test->print_read_data( $field, $test->single_fh() ); }
             'Print read data';

    $field->[1] = 4;    # Change the strand.
    $test->print_read_data( $field, $test->single_fh() );

    close $test->single_fh();
    close $test->bam_tags_fh();

    ok( -e $fastq, 'Created fastq output' );
    my $contents = read_file($fastq);
    is(
        $contents,
        ( q{@} . qq{read_name\nAGTTGTGTC\n+\nAB>BBABCB\n} ) x 2,
        'Fastq contents as expected'
    );
    unlink $fastq;

    ok( -e $test->bam_tags_file(), 'Created BAM tag output' );
    $contents = read_file( $test->bam_tags_file() );
    is(
        $contents,
          qq{read_name\t1\tRG:Z:1\tQT:Z:HIIIIIII\tRT:Z:ATCACGTT\n}
        . qq{read_name\t0\tRG:Z:1\tQT:Z:HIIIIIII\tRT:Z:ATCACGTT\n},
        'BAM tag file contents as expected'
    );
    unlink $test->bam_tags_file();
}


# Test the method that adds the pipe commands to the @PG records.
{
    my @copy = %$init;
    my %copy_hash = @copy; 
    $copy_hash{samtools_cmd} = q[t/data/aligners/bin/samtools];
    my $test = npg_common::bam_align->new(%copy_hash);

    wipe_test_scratch($scratch);
    write_file( $test->pg_records() );
    $test->sortsam_command(q{sortsam command here});
    $test->samview_command(q{samview command here});

    lives_ok { $test->add_pipe_to_pg_records(); }
             'No errors adding pipe PG records';

    my $update = read_file( $test->pg_records() );
    my $expect = q{@} . qq{PG\tID:sort_sam\tPN:samtools\tVN:0.1.11 (r851)\t}
               . qq{CL:sortsam command here\n}
               . q{@} . qq{PG\tID:samtools_view\tPN:samtools\tPP:sort_sam\t}
               . qq{VN:0.1.11 (r851)\tCL:samview command here\n};

    is( $update, $expect, 'Update verified' );

    wipe_test_scratch($scratch);
}


# Test method that finds the reference.
{
    my $test  = npg_common::bam_align->new($init);
    my $repos = 't/data/references';

    ok( $test->has_index_base(), 'Test prerequisite' );

    my $old_base    = $test->index_base();
    my $old_species = $test->species();

    lives_ok { $test->look_for_reference(); }
             'We already know the index base';
    is( $old_base, $test->index_base(), 'And we haven\'t changed it' );


    $test->clear_species();
    $test->clear_index_base();
    throws_ok
        { $test->look_for_reference(); }
        qr/Organism not passed as option nor determined from BAM header/ms,
        q[Croak if we haven't been able to identify the species];
    $test->index_base($old_base);
    $test->species($old_species);


    my $human_repos = "$repos/Homo_sapiens/NCBI36/all";

    local $ENV{CLASSPATH} = 't/bin/aligners/picard/current';
    my $common_name_test = npg_common::bam_align->new(
            {
                species        => 'Human',
                strain         => 'NCBI36',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    lives_ok { $common_name_test->look_for_reference(); }
             'Can look up species by common name';

    is(
        $common_name_test->index_base(),
        "$human_repos/bwa/ncbi36.fa",
        'And we correctly construct the index files\' basename'
    );

    my $sci_name_test = npg_common::bam_align->new(
            {
                species        => 'Homo sapiens',
                strain         => 'NCBI36',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    lives_ok { $sci_name_test->look_for_reference(); }
             'Can also use the scientific name';

    is(
        $sci_name_test->index_base(),
        "$human_repos/bwa/ncbi36.fa",
        'Again we correctly construct the index files\' basename'
    );

    my $no_strain_test = npg_common::bam_align->new(
            {
                species        => 'Human',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    lives_ok { $no_strain_test->look_for_reference(); }
             'Can look up species without a specified strain';

    is(
        $no_strain_test->index_base(),
        "$human_repos/bwa/ncbi36.fa",
        'And the index files\' basename is correct still'
    );


    # The next test requires a fresh object.
    my $two_fastas_test = npg_common::bam_align->new(
            {
                species        => 'Human',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    write_file("$human_repos/fasta/temp.fasta");

    throws_ok { $two_fastas_test->look_for_reference(); }
              qr/Failure in determining index basename from /ms,
              'Baulk if more than one fasta file is found';


    unlink "$human_repos/fasta/temp.fasta";


    my $missing_species_test = npg_common::bam_align->new(
            {
                species        => 'Dragon',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    throws_ok { $missing_species_test->look_for_reference(); }
              qr/Dragon is not represented in the repository/ms,
              'Croak if we can\'t find the species requested';

    my $missing_strain_test = npg_common::bam_align->new(
            {
                species        => 'Human',
                strain         => 'made_up_strain',
                repository     => 't/data',
                out_bam        => 'out.bam',
            }
    );

    throws_ok { $missing_strain_test->look_for_reference(); }
              qr/Human has no made_up_strain variant in the repository/ms,
              'Croak if we can\'t find the strain requested';
}


# Test conversion from fastq path to sai path
{
    my $hold = $init->{temp_dir};
    $init->{temp_dir} = 'somewhere/else';
    $init->{repository} = 't/data';
    my $test = npg_common::bam_align->new($init);

    my $sai;
    lives_ok { $sai = $test->get_sai_name_from_fastq("$scratch/test.fastq") }
             'Converts fastq path to sai path';
    is( $sai, 'somewhere/else/test.sai', 'Conversion is correct' );

    $init->{temp_dir} = $hold;
}


# Test the method that calls bwa aln.
{
    my $test = npg_common::bam_align->new($init);
    $init->{repository} = 't/data';
    wipe_test_scratch($scratch);
    write_file( $test->pg_records() );

    lives_ok { $test->bwa_aln(); } 'Call bwa aln';
    my @list = glob "$scratch/*";
    is( scalar @list, 2, 'Does nothing if there are no fastq files' );

    write_file( $test->read_1_fastq() );
    write_file( $test->read_2_fastq() );
    write_file( $test->single_fastq() );

    $test->bwa_aln();

    @list = glob "$scratch/*.sai";
    cmp_bag(
            \@list,
            [
                qq{$scratch/read_1.sai},
                qq{$scratch/read_2.sai},
                qq{$scratch/single.sai}
            ],
            'Generate sai files for every fastq found'
    );
   wipe_test_scratch($scratch);
}

# Test the method that calls bwa samse.
{
    my $test = npg_common::bam_align->new($init);
    write_file( $test->pg_records() );
    write_file( $scratch . q{/single.sai} );
    
    my $command_out_fh = $test->bwa_samse();
 
    isa_ok($command_out_fh, 'GLOB', 'file handler return');
    
    my @lines = <$command_out_fh>;
    is(scalar @lines, 0, 'no lines in test bwa samse output stream');
    close $command_out_fh;

    my $contents = read_file( $test->pg_records());    
    my $command = qq{/bin/true samse $idx_base} . qq{ $scratch/single.sai } . $test->single_fastq();    
    like($contents, qr/[@]PG\tID:bwa_samse\tPN:bwa\tPP:bam_aln\tVN:.+\tCL:$command/, 'correct bwa samse command stored in pg list');
 
    wipe_test_scratch($scratch);
}


# Test the method that calls bwa sampe.
{
    my $test = npg_common::bam_align->new($init);
    write_file( $test->pg_records() );
    write_file( $scratch . q{/read_1.sai} );
    write_file( $scratch . q{/read_2.sai} );

    my $command_out_fh = $test->bwa_sampe();
    
    isa_ok($command_out_fh, 'GLOB', 'file handler return');

    my @lines = <$command_out_fh>;
    is(scalar @lines, 0, 'no lines in test bwa sampe output stream');
    close $command_out_fh;

    my $contents = read_file( $test->pg_records());    
    my $command = qq{/bin/true sampe $idx_base} 
                  . qq{ $scratch/read_1.sai $scratch/read_2.sai }
                  . $test->read_1_fastq()
                  . q{ } . $test->read_2_fastq();

    like($contents, qr/[@]PG\tID:bwa_sampe\tPN:bwa\tPP:bam_aln\tVN:.+\tCL:$command/, 'correct bwa sampe command stored in pg list');
            
    wipe_test_scratch($scratch);
}

# Test head reconstruction.
{
    my $test = npg_common::bam_align->new($init);
    wipe_test_scratch($scratch);

    # Just to silence warnings.
    $test->sortsam_command(q{sortsam command here});
    $test->samview_command(q{samview command here});
    write_file( $test->pg_records(), q{@} . qq{PG\taligner pg line here\n} );

    my $old_head = q{@} . qq{HD\tVN:old\n}
                 . q{@} . qq{SQ\tSN:old\tLN:1\n}
                 . q{@} . qq{RG\tID:old\tSM:old\n}
                 . q{@} . qq{PG\tID:old\n}
                 . q{@} . qq{CO\told_comment\n};

    my $expect = q{@} . qq{RG\tID:old\tSM:old\n}
               . q{@} . qq{PG\taligner pg line here\n}
               . q{@} . qq{CO\told_comment\n}
               . q{@} . qq{CO\tnpg_common::bam_align (v. }
               . $npg_common::bam_align::VERSION
               . qq{) used to realign from BAM\n};

    throws_ok { $test->make_new_header(); }
              qr/Old alignment header not found/ms,
              'Need the old header to make a new header';

    write_file( $test->old_header_file(), $old_head );

    lives_ok { $test->make_new_header(); } 'Make a new header';
    ok( -e $test->new_sam_file(), 'The file is created' );

    my $merged_content = read_file( $test->new_sam_file() ); 
    is( $merged_content, $expect, 'And has the correct contents' );

    local $ENV{CLASSPATH} = 't/bin/aligners/picard/current';
    my $samfail = npg_common::bam_align->new(
                    {
                        samtools_cmd => '/bin/false',
                        scratch  => $scratch,
                        temp_dir => $scratch,
                        out_bam  => 'out.bam',
                        repository => 't/data',
                    }
    );
    wipe_test_scratch($scratch);

    my $comm_init = $init;
    $comm_init->{comment} = 'Expect to see this comment';
    $test = npg_common::bam_align->new($comm_init);
    $test->sortsam_command(q{sortsam command here});
    $test->samview_command(q{samview command here});
    write_file( $test->pg_records(), q{@} . qq{PG\taligner pg line here\n} );
    write_file( $test->old_header_file(), $old_head );

    lives_ok { $test->make_new_header() } 'Constructer accepts comments';

    $expect .= q{@} . qq{CO\tExpect to see this comment\n};
    $merged_content = read_file( $test->new_sam_file() );
    is( $merged_content, $expect, 'User comment saved to header' );

    wipe_test_scratch($scratch);
}


# Test main loop for tag rejoining.
{
    my $test = npg_common::bam_align->new($init);

    throws_ok { $test->rejoin_tags(); } qr/New alignment not found/ms,
              'Need the new alignment to re-attach the tags';

    my $test_sam = $scratch.q{/test.sam};
    write_file( $test_sam );

    open my $sam_out_fh, '<', $test_sam;
    throws_ok { $test->rejoin_tags($sam_out_fh); } qr/Old BAM tags file not found/ms,
              'Need the list of old tags to re-attach them';
    close $sam_out_fh;
    wipe_test_scratch($scratch);
}

SKIP: {
  skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.', 13 unless ($ENV{'TOOLS_INSTALLED'});

  {
    my $test;
    lives_ok { $test = npg_common::bam_align->new($init)} 'creating test object';
    write_file( $test->pg_records() );
    write_file( $test->old_header_file() );
    write_file( $test->bam_tags_file() );
    write_file( $test->new_sam_file() );

    my $test_sam = $scratch.q{/test.sam};
    write_file( $test_sam );
    open my $sam_out_fh, '<', $test_sam or die "Failed to open $test_sam";

    lives_ok { $test->rejoin_tags($sam_out_fh);} 'Call the rejoin tags method';

    wipe_test_scratch($scratch);
  }

  # Test re-attachment of BAM tags.
  {
    my $align = qq{abc\t80\tdef\t3\t5\t12I\tghi\t40\t100\tACGT\t5682\tAS:i:3};
    my $tag   = qq{def/1\t1\tAS:i:5\tU2:Z:3456\tMG:i:4\tE2:Z:GGTA\tOQ:Z:2468};

    throws_ok { npg_common::bam_align::re_tag(); }
              qr/BAM alignment line not supplied/ms,
              'Need an alignment to re-tag';
    throws_ok { npg_common::bam_align::re_tag($align); }
              qr/Old tag list not supplied/ms,
              'Need tags to re-tag';

    throws_ok { npg_common::bam_align::re_tag( $align, $tag ); }
              qr/Re-aligned reads and old tag list not synced/ms,
              'Tags are from the wrong read';

    $tag =~ s/^def/abc/msx;
    is(
        npg_common::bam_align::re_tag( "$align\n", "$tag\n" ),
        $align . qq{\tU2:Z:3456\tMG:i:4\tE2:Z:GGTA\tOQ:Z:2468\n},
        'Re-tag when the strand hasn\'t changed'
    );

    $align =~ s{\t80\t}{\t64\t}msx;

    is(
        npg_common::bam_align::re_tag( "$align\n", "$tag\n" ),
        $align . qq{\tU2:Z:6543\tMG:i:4\tE2:Z:TACC\tOQ:Z:8642\n},
        'Re-tag when the strand has changed'
    );

    $tag .= qq{\tbadly:formatted:tag};
    throws_ok { npg_common::bam_align::re_tag( "$align\n", "$tag\n" ); }
              qr/Badly formed tag: badly:formatted:tag/ms,
              'Croak on badly formed tags';
  }

  # Test the final output.
  {
    my $test = npg_common::bam_align->new($init);
    write_file($test->pg_records(), q{@}.qq{PG\tID:samtools_view\n});

    my $out_bam_fh;

    lives_ok{ $out_bam_fh = $test->output_bam(); } 'output bam pipe';
    isa_ok( $out_bam_fh, 'GLOB', 'output bam file handler returned' );

    my $pg_contents = read_file( $test->pg_records());    

    like($pg_contents, qr/[@]PG\tID:samtools_view_2\tPN:samtools\tPP:samtools_view.+[@]PG\tID:samtools_fixmate\tPN:samtools\tPP:samtools_view_2.+[@]PG\tID:sort_sam\tPN:samtools\tPP:samtools_fixmate.+/mxs, 'correct samtools sam2bam , fixmate and sort commands stored in pg list');

    close $out_bam_fh;
    wipe_test_scratch($scratch);
  }

  {
    my $test = npg_common::bam_align->new($init);
    my $out_bam = $scratch.q{/out.bam};
    $test->out_bam($out_bam);
    my $temp_bam = $test->temp_bam();
    is( $temp_bam, $scratch.q{/duplicates_unmarked.bam}, 'correct temp bam file');
    like ( $test->picard_markduplicates_command(), qr/bammarkduplicates I=t\/data\/scratch\/duplicates_unmarked\.bam O=\/dev\/stdout tmpfile=t\/data\/scratch\/ M=t\/data\/scratch\/out\.bam\.metrics/, 'correct picard markduplicates command');
    wipe_test_scratch($scratch);
  }
} # end skip

1;
