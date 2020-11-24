#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;
use PDF::API2::Simple;

### Settings ###
my $url = 'http://genome.ucsc.edu/cgi-bin/hgCustom';
my $track_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks';
my $first_row_x = 130.0;   # horizontal position where first row of chromosomes start 
my $second_row_x = 1191.8; # horizontal position where second row of chromosomes start
my $units_per_mb = 3.7;  # horizontal scale, units per megabase
my $chr1_y = 546.5;         # vertical position where first chromosomes start (chr1, chr13)
my $chr_width = 46.36;     # vertical scale, units per one chromosome
my $haplotype_shift = 32.0; # shift down so many units when on second haplotype '-'

my $usage = "Usage: perl $0 <bed_file> <pdf_file> <male/female> [hgsid]\n";
my $bed_in = shift or die $usage;
my $pdf_in = shift or die $usage;
my $sex = shift or die $usage;
die "$sex should be either male or female (case-sensitive)" unless $sex eq 'male' or $sex eq 'female';
my $hgsid = shift;


my $ua = LWP::UserAgent->new( timeout => 10 );
$ua->env_proxy;

if ( !defined( $hgsid ) ) {

    warn "Fetching session id...\n";
    my $response = $ua->get( $url );
    if ( $response->is_success ) {
        my $content = $response->content;
        ( $hgsid ) = $content =~ m/hgsid\=([^\"]+)/;
        if ( $hgsid ) {
            warn "Session ID: $hgsid\n";
        } else {
            die "Did not find hgsid in response";
        }
    }
    else {
        die "request 1 failed:", $response->status_line,"\n";
    }
}

warn "Submitting track...\n";
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

warn "Insering URLs to PDF...\n";
my $pdf = PDF::API2::Simple->open(open_file => $pdf_in);
$pdf->add_font('VerdanaBold');
$pdf->add_font('Verdana');
open F, $bed_in;
while ( <F> ) {
    chomp;
    my ( $chr, $start, $end, $name, $score, $strand ) = split /\t/;
    my $middle = int(( $start + $end )/2);
    my $size = $end-$start+1;
    my $start_show = $middle - $size; # half inversion size before
    my $end_show = $middle + $size; # half inversio size after
    my $link = $track_url.'?db=hg38&hgsid='.$hgsid.'&position='.$chr.'%3A'.$start_show.'%2D'.$end_show;
    my ( $x, $y );
    if ( $chr =~ m/^chr(1|2|3|4|5|6|7|8|9|10|11|12)$/ ) { # Chromosomes 1-12, first row
        $x += $first_row_x + $units_per_mb * $middle / 1_000_000;
        $y = $chr1_y - $chr_width * ( $1 - 1 );
    }
    elsif ( $chr =~ m/^chr(13|14|15|16|17|18|19|20|21|22|X|Y)$/ ) { # Chromosomes 13-22,X,Y second row
        $x += $second_row_x + $units_per_mb * $middle / 1_000_000;
        if ( $1 eq 'X' ) {
            if ( $sex eq 'female' ) {
                $y = $chr1_y - $chr_width * 10 ;
            }
            else {
                $y = $chr1_y - $chr_width * 9.43;
            }
        }
        elsif ( $1 eq 'Y' ) {
            $y = $chr1_y - $chr_width * 10.43;
        }
        else {
            $y = $chr1_y - $chr_width * ( $1 - 13 ); 
        }
    }
    else {
        warn "No rule for chromosome $chr, link won\'t be placed\n";
        next;
    }
    $y -= $haplotype_shift if defined($strand) and $strand eq '-';
    $pdf->link($link, '  ', x => $x, y => $y, font_size => 5, fill_color => 'red');
}
close F;

warn "Saving PDF...\n";
my $pdf_out = $pdf_in;
$pdf_out =~ s/pdf$//i;
$pdf->save($pdf_out.''.$hgsid.'.pdf');
warn "Done!\n";


print "$hgsid"
