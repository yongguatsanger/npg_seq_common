#!/usr/bin/env perl

###########################################################################################################
# seqchksum_merge.pl:
#  combines the output of bamseqchksum (crc32prod) files. When used with a set of bam/cram/sam files
#  split from a common source, the result should be identical with a seqchksum of the original.
#  Combination methods can be specified by column using the command line flags:
#    -a : accumulate (add values)
#    -c : chksum (apply crc32prod check sum, i.e. combine by multiplication modulo the prime number 2^31-1)
#    -m : match (value for a given row/col position across input files is constant)
#    -n : no check, value in the initial input file for the column is copied to the output; comment
#               rows (where the row in the initial input file begins with '#') will also behave this way
###########################################################################################################

use strict;
use warnings;
use autodie;
use Carp;
use Getopt::Std;
use Readonly;
use English qw(-no_match_vars);
use File::Basename;
use bignum;

our $VERSION = '0';

Readonly::Scalar my $MINUS_ONE => -1;
Readonly::Scalar my $CRC32PROD_MOD => 0x7fffffff; # 2^31-1

Readonly::Scalar my $NOOP => 0;
Readonly::Scalar my $ACC => 1;
Readonly::Scalar my $CONST => 2;
Readonly::Scalar my $CHKSUM => 3;

my %fnc_list = (
	$NOOP => \&noop_foldin,
	$ACC => \&addup_foldin,
	$CONST => \&const_foldin,
	$CHKSUM => \&chksum_foldin,
);

my %flag_fnc = (
	a => $ACC,
	c => $CHKSUM,
	m => $CONST,
	n => $NOOP,
);
my %opts;
getopts('a:c:m:n:h', \%opts);

# default function map matches the layout of bamseqchksum output (version 0.0.183 currently)
my @fnc_map = ();
my @default_fnc_map = (
	$CONST,
	$CONST,
	$ACC,
	$CONST,
	$CHKSUM,
	$CHKSUM,
	$CHKSUM,
	$CHKSUM,
);

if($opts{h}) {
	croak basename($PROGRAM_NAME), qq{ [-a <accumulate_fld1,accumulate_fld2,...] [-c <chksum_fld1,chksum_fld2,...] [-m <const_fld1,const_fld2,...] [-n <noop_fld1,noop_fld2,...] [-h]\n};
}

# using values from the command-line flags, initialise the fnc_map vector specifying each column's merge function
for my $flag (keys %flag_fnc) {
	set_column_fnc(\@fnc_map, $opts{$flag}, $flag_fnc{$flag});
}

my $outrows; # merged results accumulated here
my $column_count;

# process the input files
for my $fn (@ARGV) {
	open my $f, q[<], $fn or croak qq[Error: Failed to open $fn for input];
	my @inrows = <$f>;
	close $f or croak qq[Error: Failed to open $fn for input];
	chomp @inrows;

	if(not $outrows) { # determine operating parameters, finish setting up fnc_map from command line flags
		($outrows, $column_count) = initialise_outrows(\@inrows, \@fnc_map, \@default_fnc_map, $fn);
	}
	else { # combine the incoming results with what has come before
		combine_inrows($outrows, \@inrows, $column_count, $fn);
	}
}

# convert any chksum numbers back to hex strings
$outrows = chksums_to_hex_strings($outrows, \@fnc_map);

# print results
for my $row (@{$outrows}) {
	if(ref $row eq q[ARRAY]) {
		$row = join "\t", @{$row};
	}

	print $row, "\n" or carp $OS_ERROR;
}

exit;

sub set_column_fnc {
	my ($fnc_map, $fld_list, $fnc) = @_;

	if(defined $fld_list) {
		for my $i (split /,/smx, $fld_list) {
			if(not defined $fnc_map[$i-1]) {
				$fnc_map[$i-1] = $fnc;
			}
			else {
				croak q[Warning: Field ], $i-1, q[ has multiple functions defined];
			}
		}
	}

	return $fnc_map;
}

sub count_tab_columns {
	my ($rows, $fn) = @_;
	my $count = $MINUS_ONE;

	for my $row (@{$rows}) {
		next if($row =~ /^\#/smx or $row eq q[]);

		if(my $n = $row =~ tr/\t/\t/) {
			if($count >= 0 and $count != $n) {
				croak q[Inconsistent column count in non-comment input rows, file: ], $fn;
			}
			else {
				$count = $n;
			}
		}
	}

	return $count + 1;
}

sub initialise_outrows {
	my ($inrows, $fnc_map, $default_fnc_map, $fn) = @_;
	my $init_rows;
	my $col_count;

	$col_count = count_tab_columns($inrows, $fn);
	if($col_count <= 0) {
		carp q[Only comment lines or empty in initial input file ], $fn, q[, nothing to do];
		exit;
	}

	# ensure all columns have a processing function (default or NOOP)
	for my $col (0..$col_count-1) {
		if(not defined $fnc_map->[$col]) {
			$fnc_map->[$col] = $default_fnc_map->[$col];
			$fnc_map->[$col] ||= $NOOP;
		}
	}

	# convert non-comment rows into arrays
	$init_rows = [ map { /^\#/smx ?  $_ :[ split /\t/smx ]} @{$inrows} ];
	if(not @{$init_rows}) {
		croak q[No input in initial input file ], $fn;
	}

	# convert hex strings to decimals in columns containing checksum values
	for my $rowidx (0..$#{$init_rows}) {
		## no critic (ControlStructures::ProhibitUnlessBlocks ControlStructures::ProhibitPostfixControls)
		next unless(ref $init_rows->[$rowidx]); # skip comment rows
		for my $colidx (0..$col_count-1) {
			if($fnc_map->[$colidx] == $CHKSUM) {
				$init_rows->[$rowidx]->[$colidx] = hex($init_rows->[$rowidx]->[$colidx]);
			}
		}
	}

	return ($init_rows, $col_count);
}

sub combine_inrows {
	my ($init_rows, $inrows, $col_count, $fn) = @_;

	# first check layout consistency with initial input
	if(@{$inrows} != @{$init_rows}) {
		croak q[Error: number of rows is not consistent between input files (first file has ], scalar @{$init_rows}, qq[ rows, $fn has ], scalar @{$inrows}, q[)];
	}
	if((my $c=count_tab_columns($inrows, $fn)) != $col_count) {
		croak q[Error: column count inconsistent between initial input and ], $fn, q[ ( ], $col_count, q[ != ], $c , q[)];
	}

	for my $rowidx (0..$#{$init_rows}) {
		my $outrow = $init_rows->[$rowidx];

		## no critic (ControlStructures::ProhibitUnlessBlocks ControlStructures::ProhibitPostfixControls)
		next unless(ref $outrow); # skip comment rows

		unless(ref $outrow eq q[ARRAY]) { croak q[Error: non-comment row should be an ARRAY ref, not ], ref $outrow; }

		# combine the incoming values using the method specified by the function map for the column
		my @inrow = (split /\t/smx, $inrows->[$rowidx]);
		for my $colidx (0..$col_count-1) {
			if(not defined ($outrow->[$colidx] = $fnc_list{$fnc_map[$colidx]}->($outrow->[$colidx], $inrow[$colidx]))) {
				croak q[Error: Failed to merge value for column ], $colidx, q[ from file ],$fn;
			}
		}
	}

	return;
}

# foldin functions
sub noop_foldin {
	my ($outval, $inval) = @_;

	return $outval;
}

sub addup_foldin {
	my ($outval, $inval) = @_;

	return $outval + $inval;
}

sub const_foldin {
	my ($outval, $inval) = @_;

	if($outval ne $inval) {
		carp q[Value not constant (], $outval, q[ ne ], $inval, q[)];

		return;
	}

	return $outval;
}

sub chksum_foldin {
	my ($outval, $inval) = @_;

	$outval *= hex $inval;
	$outval %= $CRC32PROD_MOD; # 2^31-1

	return $outval;
}

# convert chksum numbers back to hex strings
sub chksums_to_hex_strings {
	my ($out, $fnc_map) = @_;

	for my $rowidx (0..$#{$outrows}) {
		## no critic (ControlStructures::ProhibitUnlessBlocks ControlStructures::ProhibitPostfixControls)
		next unless(ref $outrows->[$rowidx]); # skip comment rows
		for my $colidx (0..$column_count-1) {
			if($fnc_map->[$colidx] == $CHKSUM) {
				$outrows->[$rowidx]->[$colidx] = sprintf '%x', $outrows->[$rowidx]->[$colidx];
			}
		}
	}

	return $out;
}

