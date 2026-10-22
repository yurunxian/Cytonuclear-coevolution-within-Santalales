use strict;

open(TARGET,"<Orthogroup.target.1022.txt") or die();;
open(LEN,"<OG.Branch_length.normalized.txt") or die();

my %target = ();

while (<TARGET>) {
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\s+/,$_;
	$target{$tmp[0]} = $tmp[1].$tmp[2];
}

print "OG\tNode\tBL\tnorBL\tZBL\tTarget\n";
while (<LEN>) {
	chomp;
	my @tmp = split /\t/,$_;
	$tmp[1] =~ s/N//gi;
	if (exists $target{$tmp[0]} and $target{$tmp[0]} =~ /\S+/) {
		if ($target{$tmp[0]} eq "PTMT") {print $_,"\tPT\n"; print $_,"\tMT\n"}
		else {print $_,"\t",$target{$tmp[0]},"\n"}
	}
	else {print $_,"\tNuclear\n"}
}