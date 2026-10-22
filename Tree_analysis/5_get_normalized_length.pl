use strict;
use List::Util qw(sum);
use Statistics::Descriptive;

#get branch length threshold for each node
open (IN,"<OG.Branch_length.txt") or die();
my %all_raw = ();
my %threshold = ();
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	if (/Node/) {next}
	else {
		push @{$all_raw{$tmp[1]}},$tmp[2];
	}
}
#IQT = Q3 - Q1
#threshold = Q3 + 1.5 * IQT (boxplot outliers)
foreach my $node (sort keys %all_raw) {
	my @len = sort {$a <=> $b} @{$all_raw{$node}};
	$threshold{$node}{"Up"} = $len[int(scalar(@len) * 0.75)] * 2.5 - $len[int(scalar(@len) * 0.25)] * 1.5;
}

#read the data
open (IN,"<OG.Branch_length.txt") or die();
my %branch_raw = ();
my %branch = ();
my %all = ();
my %average = ();
my %variance = ();
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	if (/Node/) {next}
	else {
		if ($tmp[2] < $threshold{$tmp[1]}{"Up"}) {push @{$branch_raw{$tmp[0]}{$tmp[1]}},$tmp[2]}
	}
}


foreach my $i (sort keys %branch_raw) {
	foreach my $j (sort keys %{$branch_raw{$i}}) {
		# for a given OG, if multiple branche lengths were obtained for a given node, average them
		$branch{$i}{$j} = sum(@{$branch_raw{$i}{$j}})/scalar(@{$branch_raw{$i}{$j}});
		push @{$all{$j}}, $branch{$i}{$j};
	}
}

foreach my $i (sort keys %all) {
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@{$all{$i}});
	$average{$i} = $stat->mean();
	$variance{$i} = $stat->standard_deviation();
}

foreach my $i (sort keys %branch) {
	foreach my $j (sort keys %{$branch{$i}}) {
		my $normalize = sprintf("%.6f",$branch{$i}{$j}/$average{$j});
		my $zscore = sprintf("%.6f",($branch{$i}{$j}-$average{$j})/$variance{$j});
		print $i,"\t",$j,"\t",$branch{$i}{$j},"\t",$normalize,"\t",$zscore,"\n";
	}
}