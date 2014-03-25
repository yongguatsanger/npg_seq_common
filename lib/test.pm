#########
# Author:        kt6
# Maintainer:    $Author$
# Created:       2009-04-17
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package test;

use strict;
use warnings;

# could have this
our $VERSION = '0';

# or this
use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

1;
__END__

=head1 NAME

npg_common::test

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

A test package for substitution of git version

=cut
