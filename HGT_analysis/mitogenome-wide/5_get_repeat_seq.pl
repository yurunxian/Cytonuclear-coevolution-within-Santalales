use strict;
use warnings;

my $head = "";
my %fasta = ();
my %output = ();
my %tmp = ();
my %candidate = ();
open (FAS,"<P-tonkinensis_mito.candidate.fasta") or die();
open (IN,"<P-tonkinensis_mito.candidate.selfblastn") or die();
open (CANDIDATE,"<P-tonkinensis_mito.candidate") or die();

while (<FAS>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

while (<CANDIDATE>) {
	chomp;
	my @tmp = split /\t/,$_;
	$candidate{$tmp[0]} = $tmp[2]."-".$tmp[3];
}

while (<IN>) {
	chomp;
	my @tmp = split /\t/,$_;
	if ($tmp[0] eq $tmp[1]) {next}
	if (exists $tmp{$tmp[0]}) {next}
	if (length $fasta{$tmp[0]} eq $tmp[3] and length $fasta{$tmp[1]} eq $tmp[3]) {
		$tmp{$tmp[1]}++;
		push @{$output{$tmp[0]}},$tmp[1];
	}
}

foreach my $i (sort keys %fasta) {
	if (exists $output{$i}) {
		my @tmp = ();
		$tmp{$i}++;
#		print "$i\t",join("\t",@{$output{$i}}),"\n";
		push @tmp, $candidate{$i};
		foreach my $j (@{$output{$i}}) {push @tmp, $candidate{$j}}
		print length($fasta{$i})," (",length($fasta{$i}) * (scalar(@{$output{$i}})+1),")","\t";
		print join(", ",@tmp),"\n";

	}
	elsif (! exists $tmp{$i}) {print length($fasta{$i}),"\t",$candidate{$i},"\n"}	
}
