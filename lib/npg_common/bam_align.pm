#########
# Author:        jo3

package npg_common::bam_align;

use Moose;
use MooseX::StrictConstructor;
use Carp;
use English qw(-no_match_vars);
use File::Basename;
use File::Spec::Functions qw(catfile catdir splitdir);
use File::Temp qw( tempfile tempdir );
use IO::File;
use IPC::Open3;
use List::Util qw(max);
use File::Slurp;
use autodie qw(:all);
use Parallel::ForkManager;
use Cwd qw(cwd abs_path);

use npg_common::sequence::BAM_MarkDuplicate;

with qw/
    MooseX::Getopt
    npg_tracking::data::reference::list
    npg_common::roles::software_location
/;

our $VERSION = '0';

## no critic (Documentation::RequirePodAtEnd)
sub run {
    my ($self) = @_;

    $Carp::Verbose = $self->debug(); ## no critic (Policy::Variables::ProhibitPackageVars )

    timer_log('Start');
    $self->read_stdin();
    $self->add_pipe_to_pg_records();

    timer_log('Input read');
    $self->look_for_reference();

    timer_log('Reference found');
    my $new_sam_out_fh = $self->call_aligner();

    timer_log('Re-aligned');
    $self->rejoin_tags($new_sam_out_fh);

    timer_log('Tags re-attached');
    $self->markduplicates();

    timer_log('Duplicates marked');

    timer_log('Stop');
    return;
}


sub timer_log {
    my ($step) = @_;
    my $time = localtime;
    print {*STDERR} "$time: $step\n" or croak $OS_ERROR;
    return;
}

Readonly::Scalar my $MARK_DUPLICATES_JAR => 'MarkDuplicates.jar';

Readonly::Scalar my $DEFAULT_BWA_OPS   => '-q 15 -t 2';
Readonly::Scalar my $PAIRED_READ_FLAG  => 0x0001;
Readonly::Scalar my $QUERY_STRAND_FLAG => 0x0010;
Readonly::Scalar my $IS_READ_1_FLAG    => 0x0040;
Readonly::Scalar my $IS_READ_2_FLAG    => 0x0080;
Readonly::Scalar my $QNAME_INDEX       =>  0;
Readonly::Scalar my $FLAG_INDEX        =>  1;
Readonly::Scalar my $SEQ_INDEX         =>  9;
Readonly::Scalar my $QUAL_INDEX        => 10;
Readonly::Scalar my $FIRST_TAG_INDEX   => 11;
Readonly::Scalar my $DEBUG             =>  0;
Readonly::Scalar my $EXIT_CODE_SHIFT   =>  8;
Readonly::Scalar my $ALIGNER_NAME      => 'bwa';

# We don't want the roles' attributes to be settable via the command line.
# But do allow alternate reference repositories and aligner collections to be
# specified.
foreach my $attr ( npg_tracking::data::reference::list->meta->get_attribute_list()) {
    next if $attr =~ m/^ (?: ref_repository ) $/msx;
    has q{+} . $attr => ( metaclass => 'NoGetopt', );
}

has '+ref_repository' => (
    documentation => 'An alternative reference collection',
);

has 'debug' => (
    is            => 'ro',
    isa           => 'Int',
    default       => $DEBUG,
    documentation => 'For timestamps, verbose carping, '
                   . 'and to keep interim files',
);

=head2 position

Lane number. An integer from 1 to 8 inclusive.

=cut
has 'position'    => (isa       => 'Int',
                      is        => 'rw',
                      required  => 0,
                     );

=head2 id_run

Run id for the lane to be checked.

=cut
has 'id_run'      => (
                       isa            => 'Int',
                       is             => 'rw',
                       required       => 0,
		      );

has 'species' => (
    is            => 'rw',
    isa           => 'Str',
    predicate     => 'has_species',
    clearer       => 'clear_species',
    documentation => 'Name of the reference species',
);


has 'strain' => (
    is            => 'rw',
    isa           => 'Str',
    predicate     => 'has_strain',
    documentation => 'Name of the reference version or organism strain',
);

# This is the reference index argument for the aligner, i.e. everything up to
# the file extension.
has 'index_base' => (
    is            => 'rw',
    isa           => 'Str',
    predicate     => 'has_index_base',
    clearer       => 'clear_index_base',
    documentation => 'Base name (incl. path) of the aligner index files',
);

has 'comment' => (
    is            => 'ro',
    isa           => 'Str',
    predicate     => 'has_comment',
    documentation => 'A comment to add to the output bam file header',
);

has 'out_bam' => (
    is            => 'rw',
    isa           => 'Str',
    required      => 1,
    documentation => 'final output bam file name with path',
);

has '_temp_bam' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
    reader      => 'temp_bam',
    documentation => 'final output bam file name with path',
);

sub _build__temp_bam {
   my $self = shift;
   return catfile( $self->scratch, q{duplicates_unmarked.bam} );
}

=head2 tag_index

Tag index. An integer from 0 to 10000. Zero means that no tag has been matched.

=cut
has 'tag_index'   => (isa        => 'Int',
                      is         => 'rw',
                      required   => 0,
                     );

=head2 alignment_filter

to show the bam file is human part or nonhuman, or phix part
default 'all', which means not split

=cut
has 'alignment_filter' => (isa        => 'Str',
                           is         => 'rw',
                           required   => 0,
                          );


# Picard's SortSam, and this whole class, has a tendency to fill /tmp areas.
# Default to /tmp if we have no choice, but try to use a scratch area.
has 'scratch' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
    documentation => 'The working directory for interim files',
);

sub _build_scratch {
    my ($self) = @_;
    return $self->temp_dir;
}


# sai creation is awfully slow on scratch, so do this on /tmp.
has 'temp_dir' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
    documentation => 'A working directory for sai creation',
);

sub _build_temp_dir {
    my ($self) = @_;
    my $cleanup = $self->debug() ? 0 : 1;
    return tempdir( CLEANUP => $cleanup );
}

has 'bam_markduplicates' => (
    is            => 'ro',
    isa           => 'NpgCommonResolvedPathJarFile',
    default       => $MARK_DUPLICATES_JAR,
    coerce        => 1,
    documentation => 'Path to Picard MarkDuplicates jar file',
);

has '_picard_markduplicates_command' => (
    is            => 'rw',
    isa           => 'Str',
    lazy_build    => 1,
    accessor      => 'picard_markduplicates_command',
    documentation => 'Picard Markduplicates command line - for the bam header',
);

sub _build__picard_markduplicates_command {
    my ($self) = @_;

    my $command = $self->picard_markduplicates_wrapper()->mark_duplicate_cmd();

    return $command;
}

has 'picard_markduplicates_wrapper' => (
    is            => 'rw',
    isa           => 'npg_common::sequence::BAM_MarkDuplicate',
    lazy_build    => 1,
    documentation => 'Picard Markduplicates wrapper object',
);

sub _build_picard_markduplicates_wrapper {
    my ($self) = @_;

    my @path = splitdir $self->out_bam;
    my $filename = pop @path;
    my $in =  @path ? catdir @path : cwd;
    my $qc_out = catfile($in, q[qc]);
    if (!-d $qc_out) {
        $qc_out = $in;
    }

    my $wrapper = npg_common::sequence::BAM_MarkDuplicate->new(
          input_bam    => $self->temp_bam(),
          output_bam   => $self->out_bam(),
          metrics_file => $self->out_bam() . q{.metrics},
          metrics_json =>  catfile( $qc_out , $filename . q{_flagstats.json} ),
          jar_file     => $self->bam_markduplicates(),
          temp_dir     => $self->temp_dir(),
          not_strip_bam_tag => 1,
          reference => $self->index_base(),
          $self->resolved_paths,
    );

    if( $self->id_run() ){
       $wrapper->id_run( $self->id_run() );
    }
    if( $self->position() ){
       $wrapper->position( $self->position() );
    }
    if( defined $self->tag_index() ){
       $wrapper->tag_index( $self->tag_index() );
    }
    if( $self->alignment_filter() ){
       $wrapper->human_split( $self->alignment_filter());
    }

    return $wrapper;
}

has '_sortsam_command' => (
    is            => 'rw',
    isa           => 'Str',
    accessor      => 'sortsam_command',
    documentation => 'Picard or smatools SortSam command line - for the bam header',
);

has '_samview_command' => (
    is            => 'rw',
    isa           => 'Str',
    accessor      => 'samview_command',
    documentation => 'Samtools view command line - for the bam header',
);

has '_input_pipe' => (
    is            => 'rw',
    reader        => 'input_pipe',
    isa           => 'Str',
    lazy_build    => 1,
    documentation => 'Pipe the BAM file through this command',
);

sub _build__input_pipe {
    my ($self) = @_;

    my $sam_args  =  q{view -o - -h -};
    my $samview   = sprintf '%s %s', $self->samtools_cmd(), $sam_args;

    my ($fh, $temp_sort_out) = tempfile( DIR => $self->scratch() );
    $sam_args  =  qq{sort -n - -o $temp_sort_out};
    my $samsort = sprintf '%s %s', $self->samtools_cmd(), $sam_args;

    $self->sortsam_command($samsort);
    $self->samview_command($samview);

    return qq{$samsort | $samview};
}


has '_is_paired' => (
    is            => 'rw',
    isa           => 'Bool',
    accessor      => 'is_paired',
    predicate     => 'has_is_paired',
    clearer       => 'clear_is_paired',    # Useful for tests.
    documentation => 'Whether the aligned reads are paired or single',
);


# Paths for the re-constructed fastq files, and filehandles for them.
foreach my $attr ( qw( read_1 read_2 single ) ) {
    has q{_} . $attr . q{_fastq} => (
        is         => 'ro',
        isa        => 'Str',
        reader     => $attr . q{_fastq},
        lazy_build => 1,
    );

    has q{_} . $attr . q{_fh} => (
        is         => 'ro',
        isa        => 'IO::File',
        reader     => $attr . q{_fh},
        lazy_build => 1,
    );
}

sub _build__read_1_fastq {
    my ($self) = @_;

    return catfile( $self->scratch(), 'read_1.fastq' );
}

sub _build__read_2_fastq {
    my ($self) = @_;

    return catfile( $self->scratch(), 'read_2.fastq' );
}

# This could really just point to 'read_1.fastq'. Leave it in case we ever
# decide to tackle mixed paired/single BAM files.
sub _build__single_fastq {
    my ($self) = @_;

    return catfile( $self->scratch(), 'single.fastq' );
}

sub _build__read_1_fh {
    my ($self) = @_;

    return IO::File->new( $self->read_1_fastq(), q{>} );
}

sub _build__read_2_fh {
    my ($self) = @_;

    return IO::File->new( $self->read_2_fastq(), q{>} );
}

sub _build__single_fh {
    my ($self) = @_;

    return IO::File->new( $self->single_fastq(), q{>} );
}


# The @PG records will not be written here
has '_old_header_file' => (
    is            => 'ro',
    isa           => 'Str',
    reader        => 'old_header_file',
    lazy_build    => 1,
    documentation => 'Write the header lines from the input bam here',
);

sub _build__old_header_file {
    my ($self) = @_;

    return catfile( $self->scratch(), 'old_header.temp' );
}


# As the script progresses we generate new @PG records. We have to make sure
# all ID fields are unique and that the PP field in the following records are
# chained to them. This is simpler if we store the @PG records seperately.
has '_pg_records' => (
    is            => 'ro',
    isa           => 'Str',
    reader        => 'pg_records',
    lazy_build    => 1,
    documentation => 'Write the header PG records here',
);

sub _build__pg_records {
    my ($self) = @_;

    return catfile( $self->scratch(), 'pg_records.temp' );
}


has '_new_sam_file' => (
    is            => 'ro',
    isa           => 'Str',
    reader        => 'new_sam_file',
    lazy_build    => 1,
    documentation => 'The name of the new alignment\'s SAM file',
);

sub _build__new_sam_file {
    my ($self) = @_;

    return catfile( $self->scratch(), 'new_alignment.sam' );
}


has '_bam_tags_file' => (
    is            => 'ro',
    isa           => 'Str',
    reader        => 'bam_tags_file',
    lazy_build    => 1,
    documentation => 'A temporary file that holds the old alignment BAM tags',
);

sub _build__bam_tags_file {
    my ($self) = @_;

    return catfile( $self->scratch(), 'bam_tags.temp' );
}

has '_bam_tags_fh' => (
    is         => 'ro',
    isa        => 'IO::File',
    reader     => 'bam_tags_fh',
    lazy_build => 1,
);

sub _build__bam_tags_fh {
    my ($self) = @_;

    return IO::File->new( $self->bam_tags_file(), q{>} );
}

has '_samtools_version' => (
    is         => 'ro',
    isa        => 'Str',
    lazy_build => 1,
);
sub _build__samtools_version {
    my ($self) = @_;
    return $self->current_version($self->samtools_cmd) || q[not known];
}

has '_bwa_version' => (
    is         => 'ro',
    isa        => 'Str',
    lazy_build => 1,
);
sub _build__bwa_version {
    my ($self) = @_;
    return $self->current_version($self->bwa_cmd) || q[not known];
}


sub read_stdin {
    my ($self) = @_;

    timer_log($self->input_pipe());

    #TODO: check input stream is in bam format
    #Picard doesn't complain it, the script will die later
    #
    #OR set -o pipefail in shell before pipe bam to this scrip
    #and then the step generating bam stream will die if any problem for the input bam
    #
    #user may give a input bam file directly, may say no need to sort

    open my $sampipe, q{-|}, $self->input_pipe();
    open my $head_fh, q{>},  $self->old_header_file();
    open my $pg_fh,   q{>},  $self->pg_records();

    $self->parse_sam( $sampipe, $head_fh, $pg_fh );

    close $sampipe; # catch error with pipefail URGENT
    close $head_fh;
    close $pg_fh;

    timer_log('Old header file generated: ' . $self->old_header_file());
    timer_log('Old PG records file generated: ' . $self->pg_records());

    if ( $self->is_paired() ) {
        $self->read_1_fh->close();
        $self->read_2_fh->close();
        timer_log('Fastq file generated: '.$self->read_1_fastq . q{ and }.$self->read_2_fastq);
    }
    else {
        $self->single_fh->close();
        timer_log('Fastq file generated: '. $self->single_fastq);
    }
    $self->bam_tags_fh->close();
    timer_log('old bam tags saved: '. $self->bam_tags_file);

    return;
}


sub parse_sam {
    my ( $self, $pipe, $head, $pg ) = @_;

    croak 'Pipe handle not supplied'           if !defined $pipe;
    croak 'Header filehandle not supplied'     if !defined $head;
    croak 'PG records filehandle not supplied' if !defined $pg;

    while (<$pipe>) {
        if (m/^ [@]SQ /imsx) {               # Find the species name.
            #TODO: AS should be checked as well to determine reference
            # However, SQ is not necessary to be unique
            next if $self->has_species();

            if (m/ .* \s+ SP:( [^\t\n]+ ) /msx) {
                my $species = $1;
                $species =~ s/\s+/_/gmsx;
                $self->species($species);
            }
            next;
        }

        ###########################################
        #these just uk10k sample and study name changing
        #Readonly::Scalar my $NEW_STUDY_NAME    => 'UK10K COHORT TWINSUK';
        #Readonly::Scalar my $NEW_SAMPLE_NAME   => 'UK10K TW1537988';
        #if (m/^ [@]RG /msx) {    # Save these.
        #    my $line = $_;
        #    $line =~ s/DS:[^\t\n]+/DS:Study $NEW_STUDY_NAME/gmsx;
        #    $line =~ s/SM:[^\t\n]+/SM:$NEW_SAMPLE_NAME/gmsx;
        #    print {$head} $line or croak $OS_ERROR;
        #    next;
        #}
        ###########################################

        if (m/^ [@] (?: RG | CO ) /msx) {    # Save these.
            print {$head} $_ or croak $OS_ERROR;
            next;
        }

        if (m/^ [@] PG /msx) {               # Save these seperately.
            print {$pg} $_ or croak $OS_ERROR;
            next;
        }

        next if m/^ [@] HD /msx;             # Don't save these.

        # If we get this far we're reading alignments.
        my @first_read = $self->parse_alignment($_);

        # Don't deal with mixed paired and single reads.
        my $this_read_is_paired =
            $first_read[$FLAG_INDEX] & $PAIRED_READ_FLAG;

        if ( $self->has_is_paired() ) {
            croak "Mixed single/paired BAM file ($first_read[$QNAME_INDEX])"
                if $self->is_paired() ne $this_read_is_paired;
        }
        else {
            $self->is_paired($this_read_is_paired);
        }

        my @second_read;
        if ( $self->is_paired() ) {
            $_           = <$pipe>;
            @second_read = $self->parse_alignment($_);

            my $first_read_first = $self->consistent_read_pairs( \@first_read, \@second_read );
            my $actual_read_ref1;
            my $actual_read_ref2;
            if($first_read_first){
               $actual_read_ref1 = \@first_read;
               $actual_read_ref2 = \@second_read;
            }else{
               $actual_read_ref2 = \@first_read;
               $actual_read_ref1 = \@second_read;
            }
            $self->print_read_data( $actual_read_ref1,  $self->read_1_fh() );
            $self->print_read_data( $actual_read_ref2 , $self->read_2_fh() );
        }
        else {
            $self->print_read_data( \@first_read,  $self->single_fh() );
        }
    }

    return;
}


sub parse_alignment {
    my ( $self, $line ) = @_;

    croak 'Alignment line not supplied' if !defined $line;
    chomp $line;

    my @fields = split m/\t/msx, $line;
    $fields[$FLAG_INDEX] += 0;    # Numeric for bitwise comparisons.

    $fields[$QNAME_INDEX] =
        $self->add_paired_suffix( @fields[ $QNAME_INDEX, $FLAG_INDEX ] );

    if ( $fields[$FLAG_INDEX] & $QUERY_STRAND_FLAG ) {
        $fields[$QUAL_INDEX] = reverse $fields[$QUAL_INDEX];
        $fields[$SEQ_INDEX]  = reverse_complement( $fields[$SEQ_INDEX] );
    }

    # Concatenate all the tags into one string so that we can retrieve the
    # whole lot from one array index. Get rid of alignment-specific tags too.
    $fields[$FIRST_TAG_INDEX] =
        $self->make_tag_string( [ @fields[ $FIRST_TAG_INDEX .. $#fields ] ] );

    delete @fields[ ( $FIRST_TAG_INDEX + 1 ) .. $#fields ];

    return @fields;
}


sub add_paired_suffix {
    my ( $self, $name, $bit_string ) = @_;

    croak 'Name argument required' if !defined $name;
    croak 'Flag argument required' if !defined $bit_string;

    # Make sure paired reads end in /[12]
    if ( $bit_string & $IS_READ_1_FLAG ) {

        if ( $name =~ m{/2$}msx ) {
            carp "Contradiction in '/2' tag and read 1 flag: $name";
            return $name;
        }

        ( $name =~ m{/1$}msx ) || ( $name .= q{/1} );
    }

    if ( $bit_string & $IS_READ_2_FLAG ) {

        if ( $name =~ m{/1$}msx ) {
            carp "Contradiction in '/1' tag and read 2 flag: $name";
            return $name;
        }

        ( $name =~ m{/2$}msx ) || ( $name .= q{/2} );
    }

    # What about the pathological, but technically possible, case of both
    # $IS_READ_1_FLAG and $IS_READ_2_FLAG being set?
    # Also it's possible for $PAIRED_READ_FLAG to be set and both those to be
    # unset. In this case $name will pass through the method unaltered. This
    # is correct behaviour.
    return $name;
}


sub make_tag_string {
    my ( $self, $tag_ref ) = @_;

    my $string = q{};
    return $string if !defined $tag_ref;

    foreach my $tag ( @{$tag_ref} ) {

        # This has disappeared altogether from v1.3 of the SAM specification.
        # croak rather than just ignore it because we don't handle this currently,
        # we may need change the ori
        croak 'This tool does not support the deprecated \'SQ\' tag'
            if $tag =~ m/^SQ:/msx;

        # Skip tags that are alignment-specific
        # Skip all tags starting with X, Y, Z, reserved fields for end users
# ***
# This regex needs a careful review
# ***
        ##no critic (ProhibitComplexRegexes)
        next if $tag =~ m{^
                            (?:
                                AS | N[HM] | [NI]H | H[012I] | C[SQMCP]
                                |
                                [XYZ]. | M[DQ] | [SA]M | [BMPU]Q
                            ):
                         }msx;
        ##use critic
        $string .= "\t$tag";
    }

    $string =~ s/^\t//msx;

    return $string;
}


sub consistent_read_pairs {
    my ( $self, $fields1, $fields2 ) = @_;

    croak 'First fields arrayref not supplied'  if ref $fields1 ne 'ARRAY';
    croak 'Second fields arrayref not supplied' if ref $fields2 ne 'ARRAY';

    my $name1 = $fields1->[$QNAME_INDEX];
    my $name2 = $fields2->[$QNAME_INDEX];

    $name1 =~ s{/[12]$}{}msx;
    $name2 =~ s{/[12]$}{}msx;

    croak "Mismatch between paired read names: $name1, $name2"
        if $name1 ne $name2;

    # Make sure that the two reads are flagged as one pair.
    my $both_flags = ( $IS_READ_1_FLAG | $IS_READ_2_FLAG );
    my $pair_flag1 = $fields1->[$FLAG_INDEX] & $both_flags;
    my $pair_flag2 = $fields2->[$FLAG_INDEX] & $both_flags;

    croak "Error in paired-read flags: $name1 "
        . $fields1->[$FLAG_INDEX] . q{,} . $fields2->[$FLAG_INDEX]
        if $pair_flag1 + $pair_flag2 != $IS_READ_1_FLAG + $IS_READ_2_FLAG;

    my $first_read_first = 1;
    if($pair_flag1 == $IS_READ_2_FLAG){
       $first_read_first = 0;
    }

    return $first_read_first;
}


sub print_read_data {
    my ( $self, $fields, $fastq_fh ) = @_;

    croak 'Fields arrayref not supplied'  if ref $fields ne 'ARRAY';
    croak 'Fastq filehandle not supplied' if !defined $fastq_fh;
    my $tag_fh = $self->bam_tags_fh();

    # We may have to reverse/complement some tags if the strand changes in the
    # new alignment so store the strand also.
    my $strand = ( $fields->[$FLAG_INDEX] & $QUERY_STRAND_FLAG ) ? 1 : 0;

    printf {$fastq_fh} "@%s\n%s\n+\n%s\n",
        @{$fields}[ $QNAME_INDEX, $SEQ_INDEX, $QUAL_INDEX ]
        or croak $OS_ERROR;

    printf {$tag_fh} "%s\t%s\t%s\n",
        $fields->[$QNAME_INDEX], $strand, $fields->[$FIRST_TAG_INDEX]
        or croak $OS_ERROR;

    return;
}


sub add_pipe_to_pg_records {
    my ($self) = @_;

    my @pg_header = read_file( $self->pg_records() );

    my $pg_id     = make_unique_pg_id( \@pg_header, 'sort_sam' ) || 'sort_sam';
    my $pp_field  = prev_command_field( \@pg_header );

    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
                . q{VN:} . $self->_samtools_version . qq{\t}
                . q{CL:} . $self->sortsam_command() . qq{\n};

    push @pg_header, $pg_line;

    $pg_id    = make_unique_pg_id( \@pg_header, 'samtools_view' );
    $pp_field = prev_command_field( \@pg_header );

    $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
             . q{VN:} . $self->_samtools_version . qq{\t}
             . q{CL:} . $self->samview_command() . qq{\n};

    write_file( $self->pg_records(), @pg_header, $pg_line );

    return;
}


sub look_for_reference {
    my ($self) = @_;

    return if $self->has_index_base();

    croak 'Organism not passed as option nor determined from BAM header'
        if !$self->has_species();


    # We could just build a path here /repos/species/strain/all/aligner/fasta
    # and test for its existence, but then we could only report failure and
    # not the specific reason why.

    my $repos_species = $self->species();
    $repos_species =~ s/\s+/_/gmsx;

    # Get the repository report
    my $report = $self->report();

    # Check if the species is a common name such as Human, Mouse, etc.
    if ( $report =~ m/^ ( [^:]+ ) [^\n]+ $repos_species $/imsx ) {
        $repos_species = $1;
    }

    # Allow for case insensitive argument (e.g. --species homo_sapiens)
    if ( $report =~ m/^ ($repos_species) : /imsx ) {
        $repos_species = $1;
    }
    else {
        croak $self->species() . ' is not represented in the repository';
    }


    # If the strain/version argument has been specified, look for it.
    my $repos_strain = q{};
    if ( $self->has_strain() ) {
        $repos_strain = $self->strain();

        if ( $report =~ m/^ $repos_species : ($repos_strain) , /imsx ) {
            $repos_strain = $1;
        }
        else {
            croak $self->species()
                . " has no $repos_strain variant in the repository";
        }
    }

    # Otherwise look for the default variant.
    elsif ( $report =~ m/^ $repos_species : ( [^,]+ ) , 1/msx ) {
        $repos_strain = $1;
    }

    # This can't be triggered. npg_tracking::data::reference::list
    # will complain about this first.
    croak 'No default strain in the repository for ' . $self->species()
        if $repos_strain eq q{};


    # The aligner files all have the full fasta filename as their basename.
    my $fasta_path = catfile( $self->ref_repository(), $repos_species,
                              $repos_strain, 'all', 'fasta' );


    my @fasta_files = glob "$fasta_path/*.*";
    @fasta_files = grep {m/[.] f (?: ast | n )?  a $/imsx} @fasta_files;
    croak "Failure in determining index basename from $fasta_path"
        if scalar @fasta_files != 1;


    my $index_base = $fasta_files[0];
    $index_base =~ s/fasta/$ALIGNER_NAME/emsx;

    my @index_files = glob $index_base . q{*};
    croak "No aligner index files found: $index_base*"
        if scalar @index_files == 0;

    $self->index_base($index_base);

    return;
}


sub call_aligner {
    my ($self) = @_;
    $self->bwa_aln();
    return $self->is_paired() ? $self->bwa_sampe() : $self->bwa_samse();
}


sub get_sai_name_from_fastq {
    my ( $self, $base ) = @_;

    my $scratch = $self->scratch();
    my $tmp     = $self->temp_dir();
    $base =~ s/fastq/sai/msx;
    $base =~ s/^ $scratch /$tmp/emsx;

    return $base;
}


sub bwa_aln {
    my ($self) = @_;

    my $command   = $self->bwa_cmd() . q{ aln } . $DEFAULT_BWA_OPS .q{ }. $self->index_base();

    my @fastq_files = ( $self->read_1_fastq(), $self->read_2_fastq(), $self->single_fastq() );
    my @fastq_files_exists =  grep { -e } @fastq_files;

    my $pm = Parallel::ForkManager->new( scalar @fastq_files_exists );
    $pm->run_on_finish(
      sub { my ($pid, $exit_code, $ident) = @_;
         if($exit_code){
            croak "Failed to execute bwa aln command: PID $pid fastq $ident";
         }
       }
    );

    my $count = 0;
    foreach my $fastq  ( @fastq_files_exists )
    {
        $count++;

        my $sai = $self->get_sai_name_from_fastq($fastq);

        my $aln_command = qq{$command $fastq > $sai};

        timer_log($aln_command);

        # Sort out the ID and PP for the bam header @PG line.
        my @pg_header = read_file( $self->pg_records() );
        my $pg_id = 'bwa_aln';
        if ( $fastq =~ m/ _ ([12]) [.] fastq $/msx ) {
            $pg_id .= q{_read} . qq{$1};
        }

        $pg_id = make_unique_pg_id( \@pg_header, $pg_id );

        my $pp_field = prev_command_field( \@pg_header )
                        || qq{PP:samtools_view\t};

        my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
                    . q{VN:} . $self->_bwa_version . qq{\t}
                    . q{CL:} . $aln_command . qq{\n};

        write_file( $self->pg_records(), @pg_header, $pg_line );

        $pm->start($count) and next;

        my $exit_code;
        {
          no autodie;
          $exit_code = system $aln_command;
        }
        $exit_code = $exit_code >> $EXIT_CODE_SHIFT;
        $pm->finish($exit_code  );
    }

    $pm->wait_all_children;

    return;
}


sub bwa_samse {
    my ($self) = @_;

    my $command = $self->bwa_cmd() . q{ samse } . $self->index_base();

    my $sai     = $self->get_sai_name_from_fastq( $self->single_fastq() );

    $command .= qq{ $sai } . $self->single_fastq();

    timer_log($command);

    open my $samse_out_fh, q{-|}, $command; ## no critic (InputOutput::RequireBriefOpen);

    # Sort out the ID and PP for the bam header @PG line.
    my @pg_header = read_file( $self->pg_records() );
    my $pg_id     = make_unique_pg_id( \@pg_header, 'bwa_samse' );
    my $pp_field  = prev_command_field( \@pg_header ) || qq{\tPP:bam_aln};

    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
                . q{VN:} . $self->_bwa_version . qq{\t}
                . q{CL:} . $command . qq{\n};

    write_file( $self->pg_records(), @pg_header, $pg_line );

    return $samse_out_fh;
}


sub bwa_sampe {
    my ($self) = @_;

    my $command = $self->bwa_cmd() . q{ sampe } . $self->index_base();

    my @sai = (
                $self->get_sai_name_from_fastq( $self->read_1_fastq() ),
                $self->get_sai_name_from_fastq( $self->read_2_fastq() ),
    );

    $command .= q{ } . join( q{ }, @sai )
              . q{ } . $self->read_1_fastq()
              . q{ } . $self->read_2_fastq();

    timer_log($command);

    open my $sampe_out_fh, q{-|}, $command; ## no critic (InputOutput::RequireBriefOpen)

    # Sort out the ID and PP for the bam header @PG line.
    my @pg_header = read_file( $self->pg_records() );
    my $pg_id     = make_unique_pg_id( \@pg_header, 'bwa_sampe' );
    my $pp_field  = prev_command_field( \@pg_header ) || qq{\tPP:bam_aln};

    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
                . q{VN:} . $self->_bwa_version . qq{\t}
                . q{CL:} . $command . qq{\n};

    write_file( $self->pg_records(), @pg_header, $pg_line );

    return $sampe_out_fh;
}


sub make_new_header {
    my ($self) = @_;

    croak 'Old alignment header not found' if !-e $self->old_header_file();

    my @old_lines  = read_file( $self->old_header_file() );
    my @pg_records = read_file( $self->pg_records() );

    ## no critic (InputOutput::RequireBriefOpen)
    open my $head, q{>}, $self->new_sam_file();

    # Try to find and use a picard dict file, otherwise use the new @SQ lines later.
    my $picard_dict = $self->index_base();

    $picard_dict =~ s{/ $ALIGNER_NAME /}{/picard/}msx;
    $picard_dict .= q{.dict};

    if ( -e $picard_dict ) {
        my $dict_contents = read_file($picard_dict);
        $dict_contents =~ s/^[@]HD[^\n]+\n//msx;
        print {$head} $dict_contents or croak $OS_ERROR;
    }

    # Add the old @RG lines - there are no new ones.
    print {$head} ( grep {m/^ [@]RG /msx} @old_lines ) or croak $OS_ERROR;

    # We've been building up the @PG lines independently.
    print {$head} @pg_records or croak $OS_ERROR;

    # Add old @CO lines followed by new ones.
    print {$head} ( grep {m/^ [@]CO /msx} @old_lines ) or croak $OS_ERROR;

    # Create a comment for this script.
    my $auto_comment = q{@} . qq{CO\t} . __PACKAGE__
                     . qq{ (v. $VERSION) used to realign from BAM\n};

    print {$head} $auto_comment or croak $OS_ERROR;

    # Finish with user comments.
    if ( $self->has_comment() ) {
        print {$head} q{@}, qq{CO\t}, $self->comment(), qq{\n}
            or croak $OS_ERROR;
    }

    close $head;

    $self->debug() || ( unlink $self->old_header_file() );

    ## use critic
    return;
}

sub rejoin_tags { ##no critic (Subroutines/ProhibitExcessComplexity)
    my ($self, $new_sam_out_fh) = @_;

    #sam stream should be given, probably from bwa sampe or samse output
    croak 'New alignment not found'     if !$new_sam_out_fh;

    #old tags from a file for each alignment
    croak 'Old BAM tags file not found' if !-e $self->bam_tags_file();
    open my $tags,  q{<},  $self->bam_tags_file(); ## no critic (InputOutput::RequireBriefOpen)

    #input stream for samtools to convert sam to bam or further 
    my $bam_out_fh = $self->output_bam();

    #make new head
    $self->make_new_header();

    #write old header to the bam output
    my $new_header_with_sq = 0;
    open my $sam_header_fh,  q{<}, $self->new_sam_file();
    while( my $line = <$sam_header_fh> ){
       if($line =~ /^[@]SQ/mxs){
          $new_header_with_sq = 1;
       }
       print {$bam_out_fh} $line or croak $OS_ERROR;
    }
    close $sam_header_fh;

    my $line = <$new_sam_out_fh>;
    while ( $line && $line =~ /^[@]/mxs ) {
        #keep HD, CO header etc
        #only keep SQ if no SQ added already
        if( ( $line =~ /^[@]SQ/mxs && !$new_header_with_sq ) || ( $line !~ /^[@]SQ/mxs && $line !~ /^[@]PG/mxs )){
           print {$bam_out_fh} $line or croak $OS_ERROR;
        }
        $line = <$new_sam_out_fh>
    }

    # Assume that there is only one alignment per read. This could change if
    # the aligner is not bwa, or if non-default options are used.
    my $alignment = $line;

    while ( $alignment ) {
        my $tag_list = <$tags>;
        print {$bam_out_fh} re_tag( $alignment, $tag_list ) or croak $OS_ERROR;
        $alignment = <$new_sam_out_fh>;
    }

    close $tags;
    close $new_sam_out_fh;
# occassionally the next close breaks during tests; needs investigating !!
    close $bam_out_fh;

    $self->debug() || ( unlink $self->bam_tags_file(), $self->new_sam_file() );

    if ( $self->is_paired() ) {
        $self->debug() || unlink $self->read_1_fastq();
        $self->debug() || unlink $self->read_2_fastq();
    }
    else {
        $self->debug() || unlink $self->single_fastq();
    }

    return;
}


# Do the heavy lifting for the above method.
sub re_tag {
    my ( $sam_line, $tag_line ) = @_;

    croak 'BAM alignment line not supplied' if !defined $sam_line;
    croak 'Old tag list not supplied'       if !defined $tag_line;

    chomp $sam_line;
    chomp $tag_line;
    my @new_align = split m/\t/msx, $sam_line;
    my @old_tags  = split m/\t/msx, $tag_line;

    $new_align[$FLAG_INDEX] += 0;

    # Make sure we're dealing with the same read.
    # This is fragile if both $IS_READ_1_FLAG and $IS_READ_2_FLAG are set.
    my $bam_name = $new_align[0];
    ( $new_align[$FLAG_INDEX] & $IS_READ_1_FLAG ) && ( $bam_name .= q{/1} );
    ( $new_align[$FLAG_INDEX] & $IS_READ_2_FLAG ) && ( $bam_name .= q{/2} );

    croak 'Re-aligned reads and old tag list not synced'
        if $bam_name ne $old_tags[0];


    # Is the read on the same strand as before?
    my $old_strand = $old_tags[1] + 0;
    my $new_strand = ( $new_align[$FLAG_INDEX] & $QUERY_STRAND_FLAG ) ? 1 : 0;
    my $strand_changed = ( $old_strand != $new_strand ) ? 1 : 0;

    my %new_tags = map { substr( $_, 0, 2 ) => 1 }
                       @new_align[ $FIRST_TAG_INDEX .. $#new_align ];

    foreach ( @old_tags[ 2 .. $#old_tags ] ) {

        # Parse the tag into its three parts. The delimiter, ':', is a legal
        # character in the value, so we can't use split m/:/msx.
        ## no critic (RegularExpressions::ProhibitCaptureWithoutTest)
        croak "Badly formed tag: $_"
            if $_ !~ m{
                        ^ ( [[:alpha:]] [[:alpha:][:digit:]] )
                        :
                        ([AifZH])
                        :
                        ([^\t\n\r]+) $
                      }msx;
        my ( $tag, $vtype, $value ) = ( $1, $2, $3 );

        ## use critic
        # Newly generated tags have priority over old ones.
        next if defined $new_tags{$tag};

# Review here

        # Some flags need to be reversed and maybe complemented also. But
        # won't catch user-defined tags [XYZ][A-Z0-9] that should be reversed.
        if ($strand_changed) {
            ( $tag =~ m/^ OQ | U2 $/msx ) && ( $value = reverse $value );
            ( $tag eq 'E2' ) && ( $value = reverse_complement($value) );
        }

        $sam_line .= qq{\t} . ( join q{:}, $tag, $vtype, $value );
    }

    return qq{$sam_line\n};
}


sub output_bam {
    my ($self) = @_;

    my $sam2bam_command = $self->samtools_cmd() . q{ view -bhS -};

    my $bam_fixmate_command = $self->samtools_cmd() . q{ fixmate - -};

    my($filename, $directories, $suffix) = fileparse($self->temp_bam(), '.bam');
    my $sort_bam_command = $self->samtools_cmd() . q{ sort - } . catfile($directories, $filename) ;

    my $command = $sam2bam_command . q{ | } . $bam_fixmate_command . q{ | } . $sort_bam_command ;

    timer_log($command);

    my @pg_header = read_file( $self->pg_records() );

    my $pg_id    = make_unique_pg_id( \@pg_header, 'samtools_view' ) || 'samtools_view' ;
    my $pp_field = prev_command_field( \@pg_header );
    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
             . q{VN:} . $self->_samtools_version . qq{\t}
             . q{CL:} . $sam2bam_command . qq{\n};
    push @pg_header, $pg_line;

    $pg_id     = make_unique_pg_id( \@pg_header, 'samtools_fixmate' ) || 'samtools_fixmate';
    $pp_field  = prev_command_field( \@pg_header );
    $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
                . q{VN:} . $self->_samtools_version . qq{\t}
                . q{CL:} . $bam_fixmate_command . qq{\n};
    push @pg_header, $pg_line;

    $pg_id     = make_unique_pg_id( \@pg_header, 'sort_sam' ) || 'sort_sam';
    $pp_field  = prev_command_field( \@pg_header );
    $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
                . q{VN:} . $self->_samtools_version . qq{\t}
                . q{CL:} . $sort_bam_command . qq{\n};
    push @pg_header, $pg_line;

    $pg_id     = make_unique_pg_id( \@pg_header, 'picard_markduplicates' ) || 'picard_markduplicates';
    $pp_field  = prev_command_field( \@pg_header );
    my $mark_duplicates_version = $self->current_version($self->bam_markduplicates) || q[not known];
    $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:Picard MarkDuplicates} . qq{$pp_field\t}
                . q{VN:} . $mark_duplicates_version . qq{\t}
                . q{CL:} . $self->picard_markduplicates_command . qq{\n};

    write_file( $self->pg_records(), @pg_header, $pg_line );

    open my $bam_out_fh, q[|-], qq[/bin/bash -c "set -o pipefail; $command"];

    return $bam_out_fh;
}

sub markduplicates {
   my ($self) = @_;

   if( ! -e $self->temp_bam()){
       croak 'unmarked bam file ' . $self->temp_bam() . ' not found';
   }

   my $command = $self->picard_markduplicates_command();

   timer_log($command);

   $self->picard_markduplicates_wrapper()->process();

   $self->debug() || ( unlink $self->temp_bam() );

   return;
}


sub reverse_complement {
    my ($sequence) = @_;

    croak "$sequence doesn't look like a nucleic acid sequence"
        if $sequence =~ m/ [^ACGTUN.] /imsx;

    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTUacgtu/TGCAAtgcaa/;

    return $sequence;
}


sub prev_command_field {
    my ($pg_lines) = @_;

    my @lines = validate_pg_arg($pg_lines);
    return q{} if scalar @lines == 0;

    my $last_line = $lines[-1];

    my $pp_name;
    if ( $last_line =~ m/ ID: ( [^\t\n]+ ) /msx ) {
        $pp_name = $1;
    }
    else {
        croak "Error inferring PP field from $last_line" if !$pp_name;
    }

    return qq{\tPP:$pp_name};
}


sub make_unique_pg_id {
    my ( $pg_lines, $base_id ) = @_;

    croak 'Base id argument required' if !$base_id;
    my @lines = validate_pg_arg($pg_lines);
    return $base_id if !scalar @lines;

    my ( $max_version, $count ) = ( 0, 0 );
    foreach my $pg_line (@lines) {
        next if $pg_line !~ m/ ID:$base_id (?: _ (\d+) )? [\t\n] /imsx;

        # We could just count how often we get here, but it's better to
        # increment one greater than the maximum version number we find.
        $count++;

        if ( defined $1 ) {
            $max_version = max( $max_version, $1 + 0 );
        }
    }

    $max_version = max( $count, $max_version );

    # The version part will never be '_1' with this code.
    return $max_version ? $base_id . q{_} . ( $max_version + 1 ) : $base_id;
}


sub validate_pg_arg {
    my ($pg_lines) = @_;

    croak 'Arrayref arguments only' if ref $pg_lines ne 'ARRAY';

    return @{$pg_lines};
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
__END__

=head1 NAME

npg_common::bam_align - align the sequence in a BAM file to a reference.

=head1 VERSION

=head1 SYNOPSIS

use npg_common::bam_align;

my $bam_aln = npg_common::bam_align->new_with_options();
$bam_aln->run();

=head1 DESCRIPTION

Provide the methods to support the useful_modules script, bam_align.

=head1 SUBROUTINES/METHODS

=head2 run

Main method. Control calling of other methods.

=head2 timer_log

Print a timestamp and text string to report how long is spent on various
stages.

=head2 read_stdin

Read the BAM from STDIN, piped through a sort-by-name method and through a
BAM-to-SAM conversion. Save the header to a temporary file. Write the sequence
data to fastq(s). If we don't already know the species try to work it out
from the header.

The file handles for the fastq and bam_tags files are closed in this method.
This causes errors if print_read_data is called independently of read_stdin.

=head2 parse_sam

Do the heavy lifting for read_stdin.

=head2 parse_alignment

Take SAM file line for an alignment and return the fields as an array,
processing individual fields where required - make sure the flag is an
integer, append '/1', or '/2' to the read name if it is paired, reverse the
quality field and reverse-complement the sequence field if the read is on the
reverse strand, remove alignment-specific tags.

=head2 add_paired_suffix

Take a query name and bitwise flag, determine if the read is paired_end and,
unless it's already present, append '/1' or '/2' to the query name. Check for
inconsistency between the flag and the suffix where it is already present.
Return the query name whether modified or not.

=head2 make_tag_string

Accept an array ref of alignment BAM tags. Remove any that match hard-coded
tag names, and return the remainder as a single, tab-delimited, string.

=head2 consistent_read_pairs

Take two array refs containing the fields of two alignment lines and make sure
they are consistent with being matching read pairs. Don't care about which one
is first and which second, just that they complement each other in that sense.

=head2 print_read_data

Take the processed fields list, as an array ref, and a handle for an output
file and print the read details in fastq format. Print the tag details to
another file.

See read_stdin for a potential gotcha. Need to rethink this.

=head2 add_pipe_to_pg_records

Add bam header PG records for the pipe commands to the temporary file that we
are using to store these. Making sure that ID fields are unique and that PP
fields chain to the ID field of the previous record.

=head2 look_for_reference

The aligner will need to know the base name to use for its indexed reference
files. Following WTSI practice, this consists of the path to the files, which
will depend on the aligner name, concatenated with the name of the original
fastq file, including its extension.

=head2 call_aligner

Call the methods required for the particular aligner chosen. This method will
have to be edited for each aligner that the class should handle - as well as
writing the aligner-specific methods themselves.

=head2 get_sai_name_from_fastq

Given a path to a fastq file, convert it to the path to the corresponding sai
file.

=head2 bwa_aln

Run the bwa aln command, add a record to the PG header file.

=head2 bwa_samse

Run the bwa samse command, add a record to the PG header file.

=head2 bwa_sampe

Run the bwa sampe command, add a record to the PG header file.

=head2 make_new_header

Merge the original header and the new header. Use the new HD and SQ lines, the
old RG lines, the PG header that we've assembled separately, and lastly the old CO lines and new CO lines.

=head2 rejoin_tags

Extract the alignments from the interim BAM file. Scroll down through the
alignments and the file with the saved old tags line by line, passing them to
another method for processing. Append the returned string to the new header
file to create a SAM file for output.

=head2 re_tag

A subroutine that takes two lines passed from rejoin_tags. Croak if they do
not correspond to the same read. Retain newly generated BAM tags in each
case, but append any remaining old tags to the new alignment line - reversing
and complementing the tag value where required. Return a new alignment line.

We can only process user-defined tags with values that are derived from base
position in the read, by hard-coding them.

=head2 output_bam

Convert the SAM file to BAM. The output goes to a file.

=head2 markduplicates

run Picard Markduplicates command

=head2 reverse_complement

This is a subroutine for internal use. Reverse the input string and base-pair
complement it. If the string doesn't look like nucleic acid sequence, croak.

=head2 pg_header_lines

We can't rely on the aligner inserting its own @PG lines in the interim
alignment. This method takes an array ref holding whatever PG lines were in
the interim alignment and adds appropriate lines for the input parsing, and
for the creation of the interim alignment. It merges them together and returns
a string of concatenated @PG lines that can be directly appended to the new
header.

=head2 prev_command_field

This subroutine extracts the value of the ID field from the last line of a
set of bam header @PG records. The only argument being a concatenated string
or array-ref of those records.

=head2 make_unique_pg_id

This subroutine takes a bam header's @PG records and a string as its input.
It searches the ID: fields of the @PG records for instances of the string, and returns a value that will be unique in the bam header, appending an
incremented version number if needed. E.g. an input string of 'bwa', might be
returned as 'bwa_1', or 'bwa_2' depending on the ID: fields present in the
@PG records.

The subroutine doesn't assume that the previous @PG records are consistent or
correctly chained, and simply counts how often the input string is used as the
ID field either bare or with m/_\d+/ appended.

The assumption of an underscore between the basename and the version number
may be a bug.

=head2 validate_pg_arg

A utility subroutine to take a single argument, make sure it's an array-ref, return the de-referenced array.

Return an empty array on an empty array-ref or a null string.

The motivation for the subroutine is so that we can pass the set of @PG
records from a bam header for some parsing operation without fussing about
what format they're in. The parsing subroutines/methods will call this one to
validate the argument they recieve.


=head1 DIAGNOSTICS

=head1 DEPENDENCIES

=over

=item npg_common::roles::software_location

=back

=head1 CONFIGURATION AND ENVIRONMENT

Many of the default file paths are based on practices in the WTSI.

=head1 INCOMPATIBILITIES

Developed with reference to the SAM format specification v1.3-r882

=head1 BUGS AND LIMITATIONS

The class will croak if the input BAM file is derived from a mix of single-
and paired-end runs. It might break if the BAM file holds more than one
alignment per read - it sorts the input by read name and expects paired reads
in consecutive lines.

A partial genome reference (e.g. chromosome 1 only) can only be specified with the index_base option. If the class has to look up the reference
repository itself it will only consider full genomes (i.e. 'all'
subdirectories).

Each aligner that the class can deal with has its own methods, so it's
important that the class know the name of the aligner.

If aligner_path is specified, but not aligner, the last element of the path
will be taken as the aligner name. Should this be an alias or non-standard
name for an aligner that the class recognises, the aligner attribute should
be passed as well.

There is no generic method for calling the aligner. Code must be written for
each one that the class should work with. Although this is fairly straight-
forward.

=head1 AUTHOR

John O'Brien, E<lt>jo3@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by John O'Brien

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
