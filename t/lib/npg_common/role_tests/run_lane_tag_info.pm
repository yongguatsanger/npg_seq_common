#############
# Created By: ajb
# Created On: 2010-06-12

package npg_common::role_tests::run_lane_tag_info;
use Moose;
use Carp;
use English qw{-no_match_vars};
use Readonly;
use File::Spec::Functions qw(splitdir);
use List::Util qw(first);

with qw{npg_tracking::illumina::run::short_info npg_tracking::illumina::run::folder};
with qw{npg_common::roles::run::lane::tag_info};

our $VERSION = '0';

has q{verbose} => ( isa => q{Int}, is => q{ro} );

sub _build_run_folder {
  my ($self) = @_;
  my $path = $self->runfolder_path();
  return first {$_ ne q()} reverse splitdir($path);
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
__END__

=head1 NAME

npg_common::role_tests::run_lane_tag_info

=head1 VERSION


=head1 SYNOPSIS

=head1 DESCRIPTION

This is a test class object for testing the role npg_common::roles::run::lane::tag_info, whilst you may want to use
some of this code, it should not be considered production as is

=head1 SUBROUTINES/METHODS

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item Carp

=item English -no_match_vars

=item Readonly

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
