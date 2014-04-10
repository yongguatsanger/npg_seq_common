#########
# Author:        Marina Gourtovaia
# Created:       03 April 2009
#

package npg_common::extractor::fastq;

use strict;
use warnings;
use Carp;
use English qw(-no_match_vars);
use Math::Round qw(round);
use Exporter qw(import);
use IO::Tee;
use Fcntl ':seek';
use File::Spec;
use File::stat;
use File::Basename;

use npg_common::fastqcheck;

use Readonly;
our $VERSION = '0';

our @EXPORT_OK = qw(
                    first_read_length
                    read_count
                    generate_equally_spaced_reads
                    get_equally_spaced_reads
                    get_equally_spaced_reads_single
                    split_reads
                    generate_cache
                    retrieve_from_cache
                    to_fasta
                   );

Readonly::Scalar our $LANES_PER_READ  => 4;
Readonly::Scalar our $NUM_ARGS  => 5;
Readonly::Scalar our $CACHE_DIR_PREFIX  => q[.npg_cache];
Readonly::Scalar our $CACHE_DIR_ENV_VAR => q[NPG_FASTQ_CACHE];
Readonly::Scalar our $SAMPLE_SIZE       => 10_000;
Readonly::Scalar our $EXT => q[fastq];
Readonly::Scalar our $CACHED_FILE_DELIM => q[.];
Readonly::Scalar our $INIT_VALUE => -1;

sub _jump_by {
    my ($sources, $count) = @_;
    foreach my $fh ( @{$sources} ) {
        my $i = 0;
        while ($i < $count ) {
            <$fh> or return 0;
            $i++;
        }
    }
    return 1;
}

sub _is_pair {

    my ($l1, $l2) = @_;
    $l1 =~ s/\/\d?//smx;
    $l2 =~ s/\/[\d|t]?//smx;
    return ($l1 eq $l2);
}


sub _line_count_by_counting {

    my $fname = shift;

    local ($ENV{PATH}) = $ENV{PATH} =~ /(.*)/smx;
    my $exe = qq[wc -l $fname];
    open my $fh, q[-|], $exe or croak $ERRNO;
    my $line_count = <$fh>;
    close $fh or carp $ERRNO;
    ($line_count) = $line_count =~ /(\d+)/smx;
    if (defined $line_count) {
        $line_count = int $line_count/$LANES_PER_READ;
    }
    return $line_count;
}


sub _cache_dir {
    my ($path, $sample_size) = @_;
    my $dir =  $ENV{$CACHE_DIR_ENV_VAR};
    if ($dir) { return $dir; }
    if (!$sample_size) { $sample_size = $SAMPLE_SIZE; }
    return File::Spec->catfile($path, join q[_], $CACHE_DIR_PREFIX, $sample_size);
}


sub _link2cache {
    my ($filename, $new_file, $sample_size) = @_;

    my $target;
    ##no critic (RequireCheckingReturnValueOfEval)
    eval {
        $target = retrieve_from_cache($filename, $sample_size);
    };
    ##use critic;
    if ($target) {
        my $result =  eval { symlink $target, $new_file; 1 };
        if ($result) {
	    my ($actual_sample_size) = $target =~ m/ [.] (\d+) $/smx;
            return $actual_sample_size;
	}
    }
    return;
}


sub _close_file_handles {
    my @handles = @_;
    foreach my $fh (@handles) {
        if ($fh) {close $fh or croak $ERRNO;}
    }
    return;
}

sub first_read_length {
    my ($fname, $zero4empty_file) = @_;

    open my $fh, q[<], $fname or croak qq[Cannot open $fname for reading];
    my $err = qq[Cannot close filehandle to $fname];
    my $lane = <$fh>;
    if (!$lane) {
        close $fh or croak $err;
        if ($zero4empty_file) { return 0; }
        croak qq[First line empty in $fname];
    }
    $lane = <$fh>;
    if (!$lane) {
        close $fh or croak $err;
        croak qq[Second line empty in $fname];
    }
    close $fh or croak $err;
    chomp $lane;
    return length $lane;
}


sub read_count {

    my $fname = shift;
    my $fqname = $fname . q[check];
    my $num_reads;
    if (-e $fqname) {
        eval {
            $num_reads = npg_common::fastqcheck->new(fastqcheck_path => $fqname)->num_reads();
            1;
	} or do {
            if ($EVAL_ERROR) {
                $num_reads = _line_count_by_counting($fname);
	    }
        };
    } else {
        $num_reads = _line_count_by_counting($fname);
    }

    if (!defined $num_reads) {
        croak "Failed to get a line count for file $fname";
    }

    return $num_reads;
}


sub generate_equally_spaced_reads {  ##no critic (ProhibitExcessComplexity)
    my ($from, $to, $num_reads) = @_;

    if (!$num_reads) {
        croak q[Number of reads should be defined];
    }

    my $num_files = scalar @{$from};
    if (!$num_files) {
        croak q[Input file array empty];
    }
    if ($num_files != scalar @{$to}) {
        croak q[Array size of from and to files should be the same];
    }

    my $actual_read_num = $INIT_VALUE;
    my $previous_actual_read_num = $INIT_VALUE;
    my $count = 0;
    while ($count < $num_files) {
        $actual_read_num =  _link2cache($from->[$count],$to->[$count],$num_reads);
        if (!$actual_read_num) {
            last;
        }
        if ($count != 0 && $actual_read_num != $previous_actual_read_num) {
            $actual_read_num = undef;
            last;
        }
        $previous_actual_read_num = $actual_read_num;
        $count++;
    }
    if (defined $actual_read_num) { return $actual_read_num; }

    my $read_count = read_count($from->[0]);

    if ($read_count < $num_reads) {
        carp q[File ] . $from->[0]  . qq[ is shorter ($read_count reads) than requested ($num_reads).];
    }

    $num_reads = $read_count <= $num_reads ? $read_count : $num_reads;
    if ($num_reads == 0) {
        return $num_reads;
    }

    ## no critic (RequireBriefOpen)
    my @sources = ();
    my @destinations = ();
    $count = 0;
    while ($count < $num_files) {
        my $file = $from->[$count];
        open my $source, q[<], $file or croak qq[Cannot open $file for reading.];
        push @sources, $source;
        $file = $to->[$count];
        open my $dest, q[>], $file or croak qq[Cannot open $file for writing.];
        push @destinations, $dest;
        $count++;
    }

    my $i = 0;
    my $jump_by = 0;
    my $line = q[];
    my $comparison = q[];
    my $end_loop = $LANES_PER_READ - 1;
    my $scale_factor = $read_count / $num_reads;
    my $read_index = 0;
    my $old_read_index = 0;

    while ($i < $num_reads) {
        $old_read_index = $read_index;
        $read_index = round($i * $scale_factor);
        $jump_by = $read_index - $old_read_index - 1;
        if ($jump_by > 0) {
            _jump_by(\@sources, $jump_by * $LANES_PER_READ ) or last;
        }

        for my $j ( 0 .. $end_loop ) {
            my $hcount = 0;
            while ($hcount < $num_files) {
                my $sh = $sources[$hcount];
                $line = <$sh>;
                if (!$line) {
                    push @sources, @destinations;
                    _close_file_handles(@sources);
                    croak q[File ] . $from->[$hcount] . qq[ is shorter than reported $read_count reads.];
		}
                if ($j == 0) {
                    if( $hcount == 0) {
                        $comparison = $line;
                    }
                    elsif (!_is_pair($comparison, $line)) {
                            push @sources, @destinations;
                            _close_file_handles(@sources);
                            croak q[Reads are out of order in ] . $from->[0] . q[ and ] .
                                     $from->[$hcount] . q[, read No ] . ($i+1);
		    }
		}
                my $dh = $destinations[$hcount];
                print {$dh} $line or croak $ERRNO;
                $hcount++;
           }
       }
        $i++;
    }

    push @sources, @destinations;
    _close_file_handles(@sources);

    if ($i < $num_reads) {
        croak qq[One of the input files is shorter than reported $read_count reads.];
    }

    return $num_reads;
}


sub generate_cache {
    my ($path, $from, $sample_size) = @_;

    if (!$path) {
        croak q[Directory in which to generate the cache is not given];
    }
    if (!@{$from}) {
        croak q[Input file array empty];
    }
    if (!$sample_size) { $sample_size = $SAMPLE_SIZE; }
    my $cache = _cache_dir($path, $sample_size);
    if (!-d $cache) {
        ##no critic (RequireCheckingReturnValueOfEval)
	eval { mkdir $cache; }
        ##use critic
    }
    if (!-d $cache) {
        croak qq[Cache directory $cache does not exist];
    }

    ## no critic (ProhibitEscapedMetacharacters)
    if($from->[0] =~ /\.bam$/smx) {
		foreach my $file (@{$from}) {
			my ($outbase,$dir,$suffix) = fileparse $file, '.bam';
			my $cmd = "fastq_summ -s $sample_size -k -o ${outbase} -d $cache $file";
			if(system $cmd) {
				croak qq[fastq_summ failed for file: $file, cmd: $cmd: $CHILD_ERROR];
			}
		}
     } else {
		my @to = ();
		foreach my $file (@{$from}) {
			my ($filename,$dir,$suffix) = fileparse $file;
			push @to, File::Spec->catfile($cache, $filename)
		}
		my $actual_sample_size = generate_equally_spaced_reads($from, \@to, $sample_size);

		foreach my $file (@to) {
			rename $file, join $CACHED_FILE_DELIM, $file, $actual_sample_size;
		}
    }

    return $cache;
}


sub retrieve_from_cache {
    my ($file_path, $sample_size) = @_;

    my $source_stats = stat $file_path or croak qq[Failed to get file stats for $file_path: $ERRNO];
    my ($filename,$path,$suffix) = fileparse $file_path;
    my $cache = _cache_dir($path, $sample_size);
    if (!-d $cache) {
        croak qq[Cache directory $cache does not exist];
    }

    my @found = glob File::Spec->catfile($cache, join $CACHED_FILE_DELIM, $filename, q[*]);
    if (scalar @found == 1) {
        my $target = $found[0];
        my $stats;
        ##no critic (RequireCheckingReturnValueOfEval)
        eval { $stats = stat $target; };
        ##use critic
        if ($stats && $stats->mtime >= $source_stats->mtime) {
            return $target;
        }
    }
    return;
}


sub split_reads { ##no critic (ProhibitExcessComplexity)
    my ($fq, $num_bases, $new_fqs) = @_;

    if (!$fq) {
        croak q[Input file name should be given];
    }
    if (!$num_bases) {
        croak q[Read lengths for the target files should be given as an array reference];
    }

    ## no critic (RequireBriefOpen ProhibitMagicNumbers ProhibitTwoArgOpen ProhibitDeepNests)
    open my $source, q[<], $fq or croak qq[Cannot open $fq for reading.];

    my $num_bases1;
    my $num_bases2;

    if (!@{$num_bases}) {
        my $second_line = <$source>;
        if ($second_line) {
            $second_line = <$source>;

            if ($second_line) {
                my $read_length = (length $second_line) - 1;
                if ($read_length % 2 != 0) {
                    croak "Odd number of bases in $fq that should be split in halves";
                }
                $num_bases1 = $read_length / 2;
            } else {
                croak qq[Only one line in $fq];
            }
        } else {
            $num_bases1 = 1;
	}
        seek $source, 0, SEEK_SET;
        $num_bases2 = $num_bases1;
    } else {
        if ($num_bases->[0] <= 0 || (scalar @{$num_bases} > 1 && $num_bases->[1] <= 0)) {
	    croak q[Target read length should be positive];
        }
        $num_bases1 = int $num_bases->[0];
        $num_bases2 = scalar @{$num_bases} > 1 ? int $num_bases->[1] : 0;
    }

    my $total_wanted = $num_bases1 + $num_bases2;

    my $dest1;
    my $dest2;

    my $fq1_fh;
    my $fq2_fh;
    my @destinations = ();

    my $fqe1 = $new_fqs && @{$new_fqs} ? $new_fqs->[0] : $num_bases1 . q[_1.fastq];
    open $fq1_fh, q[>], $fqe1 or croak qq[Cannot open $fqe1 for writing.];
    push @destinations, $fq1_fh;
    $dest1 = IO::Tee->new(@destinations);

    if ($num_bases2) {
        @destinations = ();
        my $fqe2 = $new_fqs && @{$new_fqs} > 1 ? $new_fqs->[1] : $num_bases2 . q[_2.fastq];
        open $fq2_fh, q[>], $fqe2 or croak qq[Cannot open $fqe2 for writing.];
        push @destinations, $fq2_fh;
        $dest2 = IO::Tee->new(@destinations);
    }

    my $count = 0;
    while (my $line = <$source>) {
        if ($line eq qq[\n]) { next; }
        my $remainder = $count % $LANES_PER_READ;
        if ($remainder == 1 || $remainder == 3) {
	    if (length $line >= $total_wanted + 1) {
                print {$dest1} substr $line, 0, $num_bases1 or croak $ERRNO;
                print {$dest1} qq[\n] or croak $ERRNO;
                if ($num_bases2) {
                    print {$dest2} substr $line, $num_bases1, $num_bases2 or croak $ERRNO;
                    print {$dest2} qq[\n] or croak $ERRNO;
	        }
	    } else {
                my $l = $count+1;
                croak qq[Line number $l in $fq is too short];
	    }
        } elsif ($remainder == 2) {
            print {$dest1} $line or croak $ERRNO;
            if ($num_bases2) {
                print {$dest2} $line or croak $ERRNO;
	    }
	} else {
            if (!$num_bases2) {
                print {$dest1} $line or croak $ERRNO;
	    } else {
                my $l = (length $line) - 1;
                my $s = substr $line, 0, $l;
                print {$dest1} $s . '/1' . "\n" or croak $ERRNO;
                print {$dest2} $s . '/2' . "\n" or croak $ERRNO;
	    }
	}
        $count++;
    }

    close $source or croak $ERRNO;
    close $fq1_fh or croak $ERRNO;
    $dest1->close;

    if (defined $dest2) {
        close $fq2_fh or croak $ERRNO;
        $dest2->close;
    }

    return;
}

sub to_fasta {
    my ($fastq, $fasta)   = @_;
    ## no critic (InputOutput::RequireBriefOpen)
    my $fa_fh;
    if ($fasta) {
        open $fa_fh, q[>], $fasta or croak qq[Cannot open $fasta for writing];
    } else {
        open $fa_fh, q[>&STDOUT] or croak q[Failed to open stdout for writing];
    }

    open my $fq_fh, q{<}, $fastq or croak qq[Cannot open $fastq for reading $ERRNO];

    while (my $line = <$fq_fh>) {
        my ($first_word) = $line =~ m/^ (\S+) /gmsx;
        print {$fa_fh} ">$first_word\n" or croak 'cannot print to the pipe';
        $line = <$fq_fh>;
        if ($line) {
            print {$fa_fh} $line or croak 'cannot print to the pipe';
	} else {
            croak q[File ended earlier than expected];
	}
        _jump_by([$fq_fh],2);
    }
    close $fq_fh or croak qq[Cannot close a handle to $fastq $ERRNO];
    close $fa_fh or croak q[Cannot close an output filehandle];
    ## use critic

    return;
}

1;

__END__

=head1 NAME

npg_common::extractor::fastq

=head1 VERSION


=head1 SYNOPSIS
This module is for extracting parts of fastq files.

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 first_read_length
  my $l = first_read_length($my_path);
  my $zero4empty_file = 1;
  first_read_length($empty, $zero4empty_file);  #returns 0
  my $zero4empty_file = 0;
  first_read_length($empty, $zero4empty_file);  #throws an error

=head2 read_count - returns a number of reads in a fastq file.
If a fastqcheck file is present alongside the first fastq file,
the read count is taken from there. The fallback is a slow
system call to wc.

  my $count = line_count($my_path);

=head2 generate_equally_spaced_reads - extract reads from an arbitrary number of fastq files;
names of the files should be given as an array refs
  generate_equally_spaced_reads($sources, $destinations, $num_reads),
$num_reads is the number of reads that should be in the output files.

=head2 split_reads - depending on input, either trims the number of bases to the requested number
or, if two numbers are given, one output file has the trimmed reads and another has the bases that
start after the first part and extend for as long as required. There should be enough bases in a read.
 If a ref to an empty array of length is given as an argument, the
source file is split in two halves.

For a read like this:

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNANNNNNNNNNANAAANNANNNNNNNNGNANNNNN
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<&&*&&&&&&&&&<&<<;&&;&&&&&&&&(&<&&&&&

  split_read("input.fastq", [37]);

gives output

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<

written to file 37_1.fastq in the current directory

and

  split_read("input.fastq", [37,37]);

gives output

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<

and

@IL14_1008:1:1:470:276/2
ANNANNNNNNNNNANAAANNANNNNNNNNGNANNNNN
+
<&&*&&&&&&&&&<&<<;&&;&&&&&&&&(&<&&&&&

written to files 37_1.fastq and 37_2.fastq in the current directory.

Files to output can be given as a third argument
  split_read("input.fastq", [37,37], ['out1', 'out2']);


=head2 generate_cache
Call

 generate_cache($path, $files, $sample_size);

to generate the cache of extracts in the default location, i.e. in the subdirectory
of the input directory path, where $files is a ref to an array of file names.
To generate the cache in any other location, set $ENV{NPG_FASTQ_CACHE} variable.
If the sample size is not given, it defaults to 10000.

=head2 retrieve_from_cache

Call

 retrieve_from_cache($file_path, $sample_size);

to get the path to the cached file or undef if the cached file does not exist;
To retrieve from a non-standard cache location, set $ENV{NPG_FASTQ_CACHE} variable.
If the sample size is not given, it defaults to 10000.

=head2 to_fasta
A fastq to fasta converter
 
 to_fasta($fastq_in); #writes to stdout
 to_fasta($fastq_in, $fasta_out);

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item Carp

=item English -no_match_vars

=item Readonly

=item Exporter

=item Fcntl

=item Math::Round qw(round);

=item npg_common::fastqcheck

=item File::Spec

=item File::stat

=item File::Basename

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author: mg8 $

=head1 LICENSE AND COPYRIGHT

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


