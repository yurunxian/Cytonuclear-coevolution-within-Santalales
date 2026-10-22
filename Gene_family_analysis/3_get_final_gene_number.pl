use strict;

open (IN,"<RRR.txt") or die();
my @species = ("Amborella","Aristolochia","Arabidopsis","Erythropalum","Malania","Balanophora","Rhopalocnemis","Macrosolen","Tolypanthus","Scurrula","Taxillus","Hachettea","Champereia","Thesium","Pyrularia","Scleropyrum","E-cupressiformis","E-neocaledonicus","S-yasi","Osyris","Korthalsella","V-album","V-ovalifolium","A-chinense","A-pini","Daenikera","Dendrotrophe","P-rigidula","P-tonkinensis");

print "#Gene name	Arabidopsis	OXPHOS	Function\t",join("\t",@species),"\n";

while (<IN>) {
	if (/#/) {next}
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	my $name = $tmp[0];
	my %gene = ();
	my %number = ();
	print $_;
	open (GENE,"<./RRR/$name/$name.candidate");
	while (<GENE>) {
		$_ =~ s/[\r\n]+//g;
		my @tmp2 = split /\t/,$_;
		$gene{$tmp2[0]} = $tmp2[1];
	}
	open (FINAL,"<./RRR/$name/$name.trim.fas");
	while (<FINAL>) {
		if (/>(\S+)/) {
			push @{$number{$gene{$1}}}, $1;
		}
	}
	
	foreach my $sp (@species) {
#		if (exists $number{$sp}) {print "\t",scalar @{$number{$sp}}}
		if (exists $number{$sp}) {print "\t",join(", ",@{$number{$sp}})}
		else {print "\t"};
	}
	print "\n";
}