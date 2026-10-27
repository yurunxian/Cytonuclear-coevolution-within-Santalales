use strict;

open (ANCHOR,"<anchorpoints.txt") or die();
open (DUP,"<Duplication_from_trees.txt") or die();
open (KS,"<../Santalum.ks.tsv");

my %KS = ();
my %duplication = ();

while (<DUP>) {
	$_ =~ s/[\r\n]+//gi;
	my @tmp = split /\t/,$_;
	if ($tmp[2] =~ /San_/ and $tmp[3] =~ /San_/) {
		my @gene1 = ();
		my @gene2 = ();
		foreach my $i (split /, /,$tmp[2]) {if ($i =~ /San_(\S+)/) {push @gene1,$1}}
		foreach my $i (split /, /,$tmp[3]) {if ($i =~ /San_(\S+)/) {push @gene2,$1}}
		foreach my $i (@gene1) {
			foreach my $j (@gene2) {
				$duplication{$i}{$j} = $tmp[1];
				$duplication{$j}{$i} = $tmp[1];
			}
		}
	}
}

while (<KS>) {
	$_ =~ s/[\r\n]+//gi;
	my @tmp = split /\t/,$_;
	$KS{$tmp[11]}{$tmp[12]} = $tmp[5]."\t".$tmp[7]."\t".$tmp[8];
	$KS{$tmp[12]}{$tmp[11]} = $tmp[5]."\t".$tmp[7]."\t".$tmp[8];
}

while (<ANCHOR>) {
	$_ =~ s/[\r\n]+//gi;
	my @tmp = split /\t/,$_;
	if (exists $duplication{$tmp[3]}{$tmp[4]}) {
		print $_,"\t",$duplication{$tmp[3]}{$tmp[4]},"\t";
	}
	else {
		print $_,"\tunclear\t";
	}
	if (exists $KS{$tmp[3]}{$tmp[4]}) {
		print $KS{$tmp[3]}{$tmp[4]},"\n";
	}
	else {
		print "\n";
	}
}