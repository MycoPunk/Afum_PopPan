#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Bio::Location::Simple;

my $folder = "cluster_mapping";
opendir( F, $folder ) || die $!;
my %clusters_unique;

foreach my $file ( readdir(F) ) {
    next unless $file =~ /(\S+)\.bed$/;
    my $stem = $1;
    next if $stem =~ /clusters_mapped/;

    my $path = File::Spec->catfile( $folder, $file );
    my $mapped =
      File::Spec->catfile( $folder, sprintf( "%s_mapped_Af293.bed", $stem ) );
    if ( -e $path && -e $mapped ) {
        open( my $fh => $mapped ) || die "cannot open $mapped";
        while (<$fh>) {
            my ( $chrom, $start, $end, $target, $frame, $strand ) = split;
            $clusters_unique{$chrom}->{ sprintf( "%d_%d", $start, $end ) }++;
        }
    }
    else {
        warn("skipping $path/$file");
    }
}

mkdir("clustering");
for my $chrom ( sort keys %clusters_unique ) {
    open( my $fh => ">clustering/$chrom.bed" ) || die $!;
    for my $loc (
        sort { $a->[1] <=> $b->[1] }
        map {
            my ( $start, $end ) = split( '_', $_ );
            [ $_, $start, $end ]
        }
        keys %{ $clusters_unique{$chrom} } ) {
        print $fh join( "\t", $chrom, $loc->[1], $loc->[2] );
    }

    #	print join("\t",$chrom, split('_',$location),
    #		   $clusters_unique{$chrom}->{$location}), "\n";
    #    }
}

sub single_linkage {
    my $locations = shift @ARGV;

    my @matrix;
    my $i = 0;
    foreach ( my $i = 0 ; $i < scalar @$locations ; $i++ ) {
        foreach ( my $j = 1 ; $j < $i ; $j++ ) {
            if ( $locations->[$i]->overlaps( $locations->[$j] ) ) {
                push @matrix, $locations->[$i]->union( $locations->[$j] );
            }
        }
    }
    # unfinished
  #  for my $l (@$locations) {
#        #$matrix;
  #  }
}
