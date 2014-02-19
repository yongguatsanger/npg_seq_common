#############
# $Id$
# Created By: ajb
# Mast Maintained By: $Author$
# Created On: 2009-10-29
# Last Changed On: $Date$
# $HeadURL$

package npg_common::roles::log;
use strict;   # here to satisfy
use warnings; # webpublish
use Moose::Role;
use POSIX qw(strftime);
use Carp;
use Cwd;
use English qw{-no_match_vars};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mxs; $r; };

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::roles::log

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  package MyPackage;
  use Moose;
  ...
  with qw{npg_common::roles::log};

  my $MyPackage = MyPackage->new({
    filename => $sFileName,
  });

  my $MyPackage = MyPackage->new({
    filename => $sFileName,
    require_memory_log => 1,
  });


=head1 DESCRIPTION

This role imports the features which enable logging to STDERR, including setting STDERR to be a file.

=head1 SUBROUTINES/METHODS

=head2 log - method to call to send log message to STDERR (and store in the memory log if require_memory_log was set to true)

  $oMyPackage->log(q{log message});

The log message will be prepended by the current localtime. This will not croak if it cannot log!

=cut

sub log {  ## no critic (Subroutines::ProhibitBuiltinHomonyms)
  my ($self, @strs) = @_;

  my $log_string = (strftime q{[%Y-%m-%dT%H:%M:%S] }, localtime) . (join q{ }, @strs) . qq{\n};
  if (!$self->has_log_file()) {
    $self->log_file();
  }

  print {*STDERR} $log_string or carp qq{Problem trying to print to STDERR \n\t $log_string};

  if ($self->require_memory_log()) {
    my $mem_log = $self->_memory_log() || q{};
    $mem_log .= $log_string;
    $self->_memory_log($mem_log);
  }
  return 1;
}


=head2 require_memory_log - Boolean which needs to be set on creation. If so, starts to store any log string into memory, appending as it goes.

=cut

has q{require_memory_log} => (isa => q{Bool}, is => q{ro},
                              documentation => q{a memory log is required},);

# the private store of the memory log
has q{_memory_log}     => (isa => q{Str},  is => q{rw}); # logs into memory, which may be wanted later

=head2 memory_log - method to retrieve the memory_log if set. Is always guaranteed to return an empty string, so if you are allowing
the user to set require_memory_log, you can always retrieve a string successfully

  my $sMemoryLog = $oMyPackage->memory_log();

Note: the only way to get something into the memory_log is to 'log' it, this is a read only method

=cut

sub memory_log {
  my ($self) = @_;
  if ($self->_memory_log()) {
    return $self->_memory_log();
  }
  return q{};
}

=head2 log_file_name - on construction only, providing this will log to this file, rather than default STDERR (STDERR becomes this file!)

=head2 log_file_path - if a path is provided, this will be where the file is written to, default is to use current working directory

  Both of these can be set with the following, after construction, however, if you have already called the log function before setting them, you cannot start logging to that file
  You will get NO warnings about this - it is up to you to ensure you have your log_file_name and log_file_path in before you start calling log

=head2 set_log_file_name

=head2 set_log_file_path

  $class->set_log_file_name($sLogFileName);
  $class->set_log_file_path($sLogFilePath);

  Both also have predicates, to check if the are set - boolean return

=head2 has_log_file_name

=head2 has_log_file_path

  $bHasLogFileName = $class->has_log_file_name();
  $bHasLogFilePath = $class->has_log_file_path();

=cut

has q{log_file_name} => (isa => q{Str}, is => q{ro}, writer => q{set_log_file_name}, predicate => q{has_log_file_name}, documentation => q{Give a name for the log_file to use, if not provided, no log_file will be created, and log will be to STDERR});
has q{log_file_path} => (isa => q{Str}, is => q{ro}, writer => q{set_log_file_path}, lazy_build => 1, documentation => q{Provide a directory path for the log_file to written to - default is current working directory [.]});
has q{log_file}      => (isa => q{Str}, is => q{ro}, lazy_build => 1, documentation => q{Do not set on construction, use log_file_name and log_file_path});

sub _build_log_file_path {
  my ($self) = @_;
  return getcwd();
}

sub _build_log_file {
  my ($self) = @_;

  if ($self->has_log_file_name()) {
    my $log_file = $self->log_file_path() . q{/} . $self->log_file_name();
    ##no critic (RequireBracedFileHandleWithPrint)
    print STDERR "log file: $log_file\n" or croak;
    open STDERR, q{>>}, $log_file or croak qq{Unable to open log file as STDERR:$EVAL_ERROR}; # writes STDERR (so logging) to a log file in the run_folder
    ##use critic
    return $log_file;
  }
  return q{no log file};
}

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item Readonly

=item POSIX

=item Cwd

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2009 Andy Brown (ajb@sanger.ac.uk)

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
