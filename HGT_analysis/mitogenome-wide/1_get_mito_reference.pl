use strict;
use warnings;

my %accession = ();
my $head = ();
my %fasta = ();

foreach my $i (glob "/data/yurx/work/New_Caledonia/work/Source/*") {
	if ($i =~ /([A-z]+)_mito.fasta/) {
		my $sp = $1;
		open (IN,"<$i") or die();
		while (<IN>) {
			$_ =~ s/[\r\n]+//gi;
			if (/>(\S+)/) {
				$head = $1;
				push @{$accession{$sp}},$head;
			}
			else {$fasta{$sp}{$head} .= $_}
		}
	}
}

foreach my $sp (sort keys %accession) {
	print "$sp\t",join(", ",@{$accession{$sp}}),"\n";
#	foreach my $seq (@{$accession{$sp}}) {
#		print ">$sp\_$seq\n",$fasta{$sp}{$seq},"\n";
#	}
}