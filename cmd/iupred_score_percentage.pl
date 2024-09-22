#!/usr/bin/perl
use strict;
use warnings;

# Directory containing the files, change it to your directory path
my $dir = "."; # Current directory, or specify a path

# Open output file
open my $out, '>', 'output.txt' or die "Cannot open output.txt: $!";
print $out "PBid\tDisordered\tOrdered\n";

# Process each file starting with 'PB.'
foreach my $file (glob "$dir/PB.*") {
    print "Processing file: $file\n";  # Debug print
	# Open the file
    
    open my $fh, '<', $file or die "Cannot open $file: $!";


    # Variables to hold counts
    my $total = 0;
    my $disordered = 0;
    
    while (<$fh>) {
        # Skip header lines and empty lines
        next if /^#/ || /^\s*$/;

        # Split the line into columns
        my @cols = split;

        # Check the IUPRED2 score
        if ($cols[2] > 0.5) {
            $disordered++;
        }
        $total++;
    }

    # Calculate percentages
    my $disordered_percent = $total ? ($disordered / $total * 100) : 0;
    my $ordered_percent = $total ? 100 - $disordered_percent : 0;

    # Extract PBid from filename
    my ($pbid) = $file =~ /(PB\.\d+\.\d+)/;

    # Write to output
    printf $out "%s\t%.1f\t%.1f\n", $pbid, $disordered_percent, $ordered_percent;

    close $fh;
}

close $out;