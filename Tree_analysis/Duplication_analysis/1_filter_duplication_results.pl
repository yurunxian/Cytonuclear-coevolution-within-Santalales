use strict;
use Array::Utils qw(:all);

open(IN,"<./Duplication_analysis/lineage.all.txt") or die();

while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	my %clade1 = ();
	my %clade2 = ();
	foreach my $i (split /, /,$tmp[2]) {
		if ($i =~ /(\S+_\d+_c\d+_g\d+)_i\S+/) {$clade1{$1}++}
		else {$clade1{$i}++}
	}
	foreach my $i (split /, /,$tmp[3]) {
		if ($i =~ /(\S+_\d+_c\d+_g\d+)_i\S+/) {$clade2{$1}++}
		else {$clade2{$i}++}
	}
	my @gene1 = keys %clade1;
	my @gene2 = keys %clade2;
	my @overlap = intersect(@gene1, @gene2);
#	if (scalar @overlap > 0) {print join(", ",@overlap),"\n"}
	my %species = ();
	foreach my $gene (@gene1) {
		if (! grep {$_ eq $gene} @overlap) {$species{get_sp_name($gene)} = 1}
	}
	foreach my $gene (@gene2) {
		if (! grep {$_ eq $gene} @overlap) {$species{get_sp_name($gene)} += 1}
	}
	my @overlap_species = grep {$species{$_} > 1} keys %species;
	if (scalar @overlap_species >= 2) {
		print $_,"\n";
	}
}

sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}