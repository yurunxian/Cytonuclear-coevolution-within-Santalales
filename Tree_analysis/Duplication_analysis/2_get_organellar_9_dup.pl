use strict;

open (DUP,"<lineage.all.filter.txt");
open (TARGET,"<../Orthogroup.target.1022.txt");

my %target = ();

while (<TARGET>) {
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	my $target = $tmp[1].$tmp[2];
	if ($target =~ /\S/) {$target{$tmp[0]} = $target}
}

while (<DUP>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (exists $target{$tmp[0]} and $tmp[1] eq "9") {
		print $_,"\t",$target{$tmp[0]},"\n";
	}
}