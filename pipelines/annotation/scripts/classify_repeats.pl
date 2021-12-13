#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $percent_id     = 95;
my $classify_table = 'repeat_family_curate.tab';
my $db_pref        = 'lib/repeats.nr_%s';
my $hits_pref      = 'lib/repeats.nr_%s.diamond_blastx.swissprot.tab';
my ( $hitsfile, $dbfile );
GetOptions(
    't|table:s' => \$classify_table,
    'p|pid:i'   => \$percent_id,
    'hits:s'    => \$hitsfile,
    'db:s'      => \$dbfile,
);

if ( !defined $hitsfile ) {
    $hitsfile = sprintf( $hits_pref, $percent_id );
}

if ( !defined $dbfile ) {
    $dbfile = sprintf( $db_pref, $percent_id );
}

my %dbhits;
open( my $fh => $hitsfile ) || die "$hitsfile: $!";
while (<$fh>) {
    my ( $q, $h, $percent_id ) = split;
    my $hitname;
    if ( $h =~ /^([^\|]+)\|([^\|]+)\|([^\|]+)/ ) {
        $hitname = $3;
    }
    $dbhits{$q}->{$hitname} = 1;
}

my %info;
open( $fh => $classify_table ) || die $!;
while (<$fh>) {
    chomp;
    my @row = split( /,/, $_ );
    $info{ $row[0] } = $row[1];
}

my $seqio = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => $dbfile
);
my $out = Bio::SeqIO->new( -format => 'fasta' );
while ( my $seq = $seqio->next_seq ) {
    my $id = $seq->display_id;
    my $skip = 0;
    if ( $id =~ /(\S+)\#Unknown/ ) {
        my $name = $1;
        for my $hit ( keys %{ $dbhits{$id} } ) {
            if ( exists $info{$hit} ) {
                if ( $info{$hit} =~ /SKIP/ ) {
                    $skip = 1;
                    warn("skipping $hit\n");
                    last;
                }
                else {
                    # overwrite name with new classification
                    $seq->display_id( sprintf( "%s#%s", $name, $info{$hit} ) );
                }
            }
        }
    }
    $out->write_seq($seq) unless $skip;
}
