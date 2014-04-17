#############
# Created By: ajb
# Created On: 2010-12-10

package npg_common::role_tests::run_status_has_run_method;
use Moose;

extends qw{npg_common::role_tests::run_status_no_extra_methods};

use Readonly;
our $VERSION = '0';

has q{run} => (
  isa => q{npg::api::run},
  is => q{ro},
  lazy_build => 1,
);

sub _build_run {
  my ( $self ) = @_;
  return npg::api::run->new();
}

sub id_run {
  return 1234;
}

=head1 NAME

=head1 VERSION

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=cut

no Moose;
__PACKAGE__->meta->make_immutable;
1;
__END__

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
