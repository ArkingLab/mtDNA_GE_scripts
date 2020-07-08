#!/usr/bin/env perl

use warnings;
use strict;

my @f = `ls *_stats.txt | cut -f1 -d "_"`;

open (OUTPUT,">stats_summary.txt");
print OUTPUT "ID\tTotal\tMapped\tUnmapped\tAutosomal\tchrX\tchrY\tchrMT\trandom\tunknown\tebv\tdecoy1\tdecoy2\n";

for my $id (@f) {
	my $total_reads=0;
	my $mapped_reads=0;
	my $unmapped_reads=0;
	my $autosomal=0;
	my $x=0;
	my $y=0;
	my $mt=0;
	my $random=0;
	my $unknown=0;
	my $ebv=0;
	my $decoy1=0;
	my $decoy2=0;
	my $i=1;
	chomp($id);
	my $file = join "_", $id, "stats.txt";
	unless (open(FILE,"$file")) {die "could not open $file\n"; }
	while(<FILE>) {
		chomp;
		my $line=$_;
		my @line = split /\s+/, $line;
		if ($i<=22) {
			$autosomal = $autosomal + $line[2];
		} elsif ($i==23) {
			$x = $line[2];
		} elsif ($i==24) {
			$y=$line[2];
		} elsif ($i==25) {
			$mt = $line[2];
		} elsif ($i>=26 && $i<=67) {
			$random = $random + $line[2];
		} elsif ($i>=68 && $i<=194) {
			$unknown = $unknown + $line[2];
		} elsif ($i==195) {
			$ebv = $line[2];
		} elsif ($i>=196 && $i<=582) {
			$decoy1 = $decoy1 + $line[2];
		} elsif ($i>=583 && $i<=2580) {
			$decoy2 = $decoy2 + $line[2];
		}
		$mapped_reads = $mapped_reads + $line[2];
		$unmapped_reads = $unmapped_reads + $line[3];
		$i++;
	}
	$total_reads = $mapped_reads + $unmapped_reads;
	print OUTPUT "$id\t$total_reads\t$mapped_reads\t$unmapped_reads\t$autosomal\t$x\t$y\t$mt\t$random\t$unknown\t$ebv\t$decoy1\t$decoy2\n";
	
}
close(OUTPUT);
exit;
