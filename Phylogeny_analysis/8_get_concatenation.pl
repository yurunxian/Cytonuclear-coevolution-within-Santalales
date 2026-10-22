use strict;

my %fasta = ();
my @order = ("Amb","Ara","Ari","Ery","Mal","Bya","Rph","Hac","Mco","Tma","Scu","Tax","Cha","The","Pyr","Scl",
	"Ecu","Ene","Sya","Osy","Ach","Api","Kor","Val","Vov","Dae","Den","Pri","Pto");

open (IN,"<filter.list") or die();

while (<IN>) {
	if (/(\S+)/) {
		my $gene = $1;
		system("cat ./BUSCO_tree/$gene/$gene.final.cds.fas.raxml.support >> final0416.coalescence.input");
		if (-s "aln/$gene.final.fas") {
			open (FAS,"<aln/$gene.final.fas") or die();
			my $head = "";
			my %tmp = ();
			my $length = ();
			while (<FAS>) {
				chomp;
				if (/>(\S+)/) {$head = $1}
				else {$tmp{$head} .= $_}
			}
			foreach my $sp (@order) {
				if (exists $tmp{$sp}) {
					$length = length $tmp{$sp};
					last;
				}
			}
			foreach my $sp (@order) {
				if (exists $tmp{$sp}) {$fasta{$sp} .= $tmp{$sp}}
				else {$fasta{$sp} .= "-" x $length}
			}
		}
	}
}

foreach my $sp (@order) {print ">$sp\n$fasta{$sp}\n"}