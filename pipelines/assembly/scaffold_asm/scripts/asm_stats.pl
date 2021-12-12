#!/usr/bin/perl
#
use File::Spec;
use strict;
use warnings;

my %stats;


my $read_map_stat = 'mapping_report';
my $dir = shift || 'scaffolds';
my %cols;
my @header;
my %header_seen;

opendir(DIR,$dir) || die $!;
my $first = 1;
foreach my $d ( readdir(DIR) ) {
    opendir(SUB,"$dir/$d") || die "$dir/$d: $!";
	my $stem = $d;
	foreach my $file ( readdir(SUB) ) {
	    next unless ( $file eq 'assembly.stats.txt');
	    open(my $fh => "$dir/$d/$file") || die "cannot open $dir/$d/$file: $!";
	    while(<$fh>) {
		next if /^\s+$/;
		s/^\s+//;
		chomp;
		if ( /\s*(.+)\s+=\s+(\d+(\.\d+)?)/ ) {
		    my ($name,$val) = ($1,$2);	    
		    $name =~ s/\s*$//;
		    $stats{$stem}->{$name} = $val;
#	    warn("'$name'"," '", $val,"'\n");
		    $cols{$name}++;
		    if( ! exists $header_seen{$name} ) {
			push @header, $name;
			$header_seen{$name} = 1;
		    }
		}
	    }
	    if ( $first ) { 
		push @header, qw(BUSCO_Complete BUSCO_Single BUSCO_Single BUSCO_Single
                 BUSCO_Fragmented BUSCO_Missing BUSCO_NumGenes
                 );
	    }


	    my $busco_file = File::Spec->catfile("BUSCO",sprintf("run_%s",$stem),
						 sprintf("short_summary_%s.txt",$stem));

	    if ( -f $busco_file ) {

		open(my $fh => $busco_file) || die $!;
		while(<$fh>) {	 
		    if (/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/ ) {
			$stats{$stem}->{"BUSCO_Complete"} = $1;
			$stats{$stem}->{"BUSCO_Single"} = $2;
			$stats{$stem}->{"BUSCO_Duplicate"} = $3;
			$stats{$stem}->{"BUSCO_Fragmented"} = $4;
			$stats{$stem}->{"BUSCO_Missing"} = $5;
			$stats{$stem}->{"BUSCO_NumGenes"} = $6;
		    } 
		}

	    } else {
		warn("Cannot find $busco_file");
	    }

	    my $sumstatfile = File::Spec->catfile($read_map_stat,
						  sprintf("%s.bbmap_summary.txt",$stem));
	    if ( -f $sumstatfile ) {
		open(my $fh => $sumstatfile) || die "Cannot open $sumstatfile: $!";
		my $read_dir = 0;
		my $base_count = 0;
		$stats{$stem}->{'Mapped reads'} = 0;
		while(<$fh>) {
		    if( /Read (\d+) data:/) {
			$read_dir = $1;
		    } elsif( $read_dir && /^mapped:\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/) {
			$base_count += $4;
			$stats{$stem}->{'Mapped_reads'} += $2;
		    }  elsif( /^Reads:\s+(\S+)/) {
			$stats{$stem}->{'Reads'} = $1;
		    }

		}
		$stats{$stem}->{'Average_Coverage'} =
		    sprintf("%.1f",$base_count / $stats{$stem}->{'TOTAL LENGTH'});
		if( $first )  {
		    push @header, ('Reads',
				   'Mapped_reads',			   
				   'Average_Coverage');
		}
	    }
	    $first = 0;
    }
}
print join("\t", qw(SampleID), @header), "\n";
foreach my $sp ( sort keys %stats ) {    
    print join("\t", $sp, map { $stats{$sp}->{$_} || 'NA' } @header), "\n";
}
