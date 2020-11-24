#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;

my $url = 'http://genome.ucsc.edu/cgi-bin/hgCustom';
my $usage = "Usage: perl $0 <bed_file> [hgsid]\n";
my $bed_in = shift;
my $hgsid = shift;


my $ua = LWP::UserAgent->new( timeout => 10 );
$ua->env_proxy;

warn "Submitting reads ...\n";
my %form = ( 'mainForm' => '',#'mainForm',
             'clade'   => 'mammal',
             'hgsid'  => $hgsid,
             'org'   => 'Human',
             'db'   => 'hg38',
             'hgt.customFile' => [ $bed_in ],
             'Submit' => 'Submit',
           );
my $response2 = $ua->post( $url,
                        Content_Type => 'multipart/form-data',
                        Content => \%form ) ;
if ( $response2->is_success ) {
    my $content2 = $response2->content;
}
else {
    die "request 2 failed:", $response2->status_line,"\n";
}


