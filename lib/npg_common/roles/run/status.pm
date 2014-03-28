#############
# $Id$
# Created By: ajb
# Last Maintained By: $Author$
# Created On: 2010-12-06
# Last Changed On: $Date$
# $HeadURL$

package npg_common::roles::run::status;
use Moose::Role;
use Carp;
use English qw{--no_match_vars};
use npg::api::run;
use npg::api::request;
use npg::api::run_status_dict;


our $VERSION = '0';

with qw{npg_common::roles::log};

## no critic (Documentation::RequirePodAtEnd)
=head1 NAME

npg_common::roles::run::status

=head1 VERSION

$Revision$

=head1 SYNOPSIS

  package My::Moose::Package;
  use Moose;
  with qw{npg_common::roles::run::status};


=head1 DESCRIPTION

This role imports some methods for obtaining and updating the run status of a run, via the npg::api::run
and npg::api::run_status_dict modules

=head1 SUBROUTINES/METHODS

=head2 update_run_status

Takes a run id and a status, and updates this run to this status.

  $oClass->update_run_status( {
    id_run => $iIdRun,
    status_desc => $sStatus
  } );

It is clever enough to look at $self->id_run() for an id_run, however, it won't look in any run objects,
and will then croak if you haven't now provided both.

status_desc must be provided in the args given

If you have a a run object already created ( in $self->run() ), it will use this (and it's util object),
else it will try for a util object in $self->apiutil(), or finally will try to obtain the objects, allowing
them to create their own util objects

=cut

sub update_run_status {
  my ($self, $args) = @_;

  local $ENV{npg::api::request->cache_dir_var_name()} = q();
  local $ENV{npg::api::request->save2cache_dir_var_name()} = 0;

  if ( ! $args->{id_run} && $self->can( q{id_run} ) ) {
    $args->{id_run} = $self->id_run();
  }

  if ( ! $args->{id_run} || ! $args->{status_desc} ) {
    croak q{no id_run and/or status_desc provided};
  }

  my $id_run = $args->{id_run};

  $self->log(qq{Updating Run $args->{id_run} to '$args->{status_desc}'});

  if ( $self->can( q{run} ) && ( ref $self->run() ) && $self->run()->id_run() == $id_run ) {
    $args->{run} = $self->run();
    $args->{util} = $self->run()->util();
    $self->_obtain_run_and_all_other_objects_for_run_statuses( $args );
  }

  if ( ! $args->{run} ) {
    $args->{util} = $self->can( q{apiutil} ) ? $self->apiutil() : undef;
    $self->_obtain_run_and_all_other_objects_for_run_statuses( $args );
  }

  my $current_run_status = $args->{current_run_status};
  my $rsd = $args->{rsd};
  my $run = $args->{run};

  if ( ! $run || ! $rsd || ! $current_run_status ) {
    croak q{Unable to obtain a necessary objects for updating the status of } . $id_run;
  }

  if ( $current_run_status->id_run_status_dict() == $rsd->id_run_status_dict() ) {
    $self->log( qq{Leaving $id_run at '$args->{status_desc}'/} . $rsd->id_run_status_dict() );
    return 1;
  }

  $self->log( qq{Updating $id_run to '$args->{status_desc}'/} . $rsd->id_run_status_dict() );
  my $new_run_status = npg::api::run_status->new({
		util               => $run->util(),
		id_run             => $id_run,
		id_run_status_dict => $rsd->id_run_status_dict(),
	});

  return $new_run_status->create();
}

=head2 run_statuses

Returns a hashref of run statuses from the npg tracking application, with a set of descriptive keys
which can remain constant

=cut

Readonly::Hash our %STATUSES => (
  6 => q{STATUS_ANALYSIS_PENDING},
  7 => q{STATUS_ANALYSIS_IN_PROGRESS},
  9 => q{STATUS_ANALYSIS_COMPLETE},
  12 => q{STATUS_RUN_ARCHIVED},
  14 => q{STATUS_ANALYSIS_PRELIM},
  15 => q{STATUS_ANALYSIS_PRELIM_COMPLETE},
  17 => q{STATUS_ARCHIVAL_PENDING},
  18 => q{STATUS_ARCHIVAL_IN_PROGRESS},
  19 => q{STATUS_QC_REVIEW_PENDING},
  20 => q{STATUS_QC_COMPLETE},
  24 => q{STATUS_SECONDARY_ANALYSIS_IN_PROGRESS},
);

sub run_statuses {
  my ( $self, $args ) = @_;

  if ( $self->_has_run_statuses() ) {
    return $self->_run_statuses();
  }

  my $rsds = $self->_obtain_rsd_object( $args )->run_status_dicts();

  my $href = {};

  foreach my $rsd ( @{ $rsds } ) {
    if ( $STATUSES{ $rsd->id_run_status_dict } ) {
      $href->{ $STATUSES{ $rsd->id_run_status_dict } } = $rsd->description();
    }
  }

  $self->_run_statuses( $href );

  return $href;
}

# cache the run statuses
has q{_run_statuses} => (
  isa => q{HashRef[Str]},
  is  => q{rw},
  predicate => q{_has_run_statuses},
);

# does the evaluation and obtaining of the objects needed
sub _obtain_run_and_all_other_objects_for_run_statuses {
  my ( $self, $args ) = @_;

  my $run = $args->{run};
  my $rsd = $args->{rsd};
  my $util = $args->{util};
  my $id_run = $args->{id_run};
  my $status_desc = $args->{status_desc};

  eval {
    # if run was not already provided
    # use any util provided
    if ( ! $args->{run} ) {
      $args->{run} = npg::api::run->new( {
        id_run  => $id_run,
        ( $util ? ( util => $util ) : ()),
      } );
    }
    # if run_status_dict not already provided
    if ( ! $args->{rsd} ) {
      $args->{rsd} = $self->_obtain_rsd_object( {
        status_desc => $status_desc,
        util        => $args->{run}->util(),
	    } );
    }
    # force a read on run to get the current run status, which will test if the object is ok
    $args->{current_run_status} = $args->{run}->current_run_status();
  } or do {
    $self->log( $EVAL_ERROR );
    undef $args->{run};
    undef $args->{rsd};
    undef $args->{util};
  };

  return 1;
}

sub _obtain_rsd_object {
  my ( $self, $args ) = @_;

  my $init = {};
  $args->{status_desc} && ( $init->{description} = $args->{status_desc} );
  $args->{util}        && ( $init->{util}        = $args->{util} );

  return npg::api::run_status_dict->new( $init );
}

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item English -no_match_vars

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Andy Brown (ajb@sanger.ac.uk)

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
