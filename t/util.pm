package t::util;

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Class::Std;
use HTML::PullParser;
use YAML qw(LoadFile);
use MIME::Parser;
use MIME::Lite;
use Test::More;
use Cwd;
use File::Temp qw{ tempdir };

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mx; $r; };

Readonly::Scalar our $TEMP_DIR => q{/tmp};

{

  my %emails_of       :ATTR( :get<emails>,        :set<emails>       );
  my %parsed_email_of :ATTR( :get<parsed_email>,  :set<parsed_email> );

  sub rendered {
    my ($self, $tt_name) = @_;
    local $RS = undef;
    open my $fh, q(<), $tt_name or croak "Error opening $tt_name: $ERRNO";
    my $content = <$fh>;
    close $fh or croak "Error closing $tt_name: $ERRNO";
    return $content;
  }

  sub test_rendered {
    my ($self, $chunk1, $chunk2) = @_;

    if(!$chunk1) {
      diag q(No chunk1 in test_rendered);
    }

    if(!$chunk2) {
      diag q(No chunk2 in test_rendered);
    }

    if($chunk2 =~ m{^t/}mx) {
      $chunk2 = $self->rendered($chunk2);

      if(!length $chunk2) {
        diag("Zero-sized $chunk2. Expected something like\n$chunk1");
      }
    }

    my $chunk1els = $self->parse_html_to_get_expected($chunk1);
    my $chunk2els = $self->parse_html_to_get_expected($chunk2);
    my $pass      = $self->match_tags($chunk2els, $chunk1els);

    if($pass) {
      return 1;

    } else {
      my ($fn) = $chunk2 =~ m{([^/]+)$}mx;
      $fn     .= q(-);
      open my $fh1, q(>), "/tmp/${fn}chunk1" or croak "Error opening /tmp/${fn}chunk1";
      open my $fh2, q(>), "/tmp/${fn}chunk2" or croak "Error opening /tmp/${fn}chunk2";
      print $fh1 $chunk1;
      print $fh2 $chunk2;
      close $fh1 or croak "Error closing /tmp/${fn}chunk1";
      close $fh2 or croak "Error closing /tmp/${fn}chunk2";
      diag("diff /tmp/${fn}chunk1 /tmp/${fn}chunk2");
    }

    return;
  }

  sub parse_html_to_get_expected {
    my ($self, $html) = @_;
    my $p;
    my $array = [];

    if ($html =~ m{^t/}xms) {
      $p = HTML::PullParser->new(
			       file  => $html,
			       start => '"S", tagname, @attr',
			       end   => '"E", tagname',
			      );
    } else {
      $p = HTML::PullParser->new(
			       doc   => $html,
			       start => '"S", tagname, @attr',
			       end   => '"E", tagname',
			      );
    }

    my $count = 1;
    while (my $token = $p->get_token()) {
      my $tag = q{};
      for (@{$token}) {
        $_ =~ s/\d{4}-\d{2}-\d{2}/date/xms;
        $_ =~ s/\d{2}:\d{2}:\d{2}/time/xms;
        $tag .= " $_";
      }
      push @{$array}, [$count, $tag];
      $count++;
    }

    return $array;
  }

  sub match_tags {
    my ($self, $expected, $rendered) = @_;
    my $fail = 0;
    my $a;

    for my $tag (@{$expected}) {
      my @temp = @{$rendered};
      my $match = 0;
      for ($a= 0; $a < @temp;) {
        my $rendered_tag = shift @{$rendered};
        if ($tag->[1] eq $rendered_tag->[1]) {
          $match++;
          $a = scalar @temp;
        } else {
          $a++;
        }
      }

      if (!$match) {
        diag("Failed to match '$tag->[1]'");
        return 0;
      }
    }

    return 1;
  }

  ###########
  # for catching emails, so firstly they don't get sent from within a test
  # and secondly you could then parse the caught email
  #

  sub catch_email {
    my ($self) = @_;
    if (!$self->get_emails()) {
      $self->set_emails([]);
    }
    my $sub = sub {
                my $msg = shift;
        	      push @{$self->get_emails()}, $msg->as_string;
                return;
	        };
    MIME::Lite->send('sub', $sub);
    return;
  }

  ##########
  # for parsing emails to get information from them, probably caught emails
  #

  sub parse_email {
    my ($self, $email) = @_;
    if (!$email) {
      $email = shift @{$self->get_emails()};
    }

    my $parser = MIME::Parser->new();
    $parser->output_to_core(1);
    my $entity = $parser->parse_data($email);
    my $ref    = {
		    body          => $entity->bodyhandle->as_string()     || undef,
		    subject       => $entity->head->get('Subject', 0)     || undef,
		    to            => $entity->head->get('To',0)           || undef,
		    cc            => $entity->head->get('Cc',0)           || undef,
		    bcc           => $entity->head->get('Bcc',0)          || undef,
		    from          => $entity->head->get('From',0)         || undef,
		    precendence   => $entity->head->get('Precedence',0)   || undef,
		    content_type  => $entity->head->get('Content-Type',0) || undef,
	    };

    return $ref;
  }
  
  my %temp_dir_of :ATTR( :get<temp_dir>, :set<temp_dir> );
  sub temp_directory {
    my ( $self ) = @_;

    if ( $self->get_temp_dir() ) {
      return $self->get_temp_dir();
    }

    my $tempdir = tempdir(
      DIR => $TEMP_DIR,
      CLEANUP => 1,
    );

    $self->set_temp_dir( $tempdir );

    note $tempdir;
    return $tempdir;
  }
}

1;
