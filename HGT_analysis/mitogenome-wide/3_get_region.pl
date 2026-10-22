use strict;
use warnings;

open (IN,"<P-tonkinensis_mito.blastn") or die(); #identity > 70; hit > 200 bp
open (FAS,"<P-tonkinensis_mito.fasta") or die();

my %clade = ("Amyema"=>"Loranthaceae","Helixanthera"=>"Loranthaceae","Scurrula"=>"Loranthaceae","Taxillus"=>"Loranthaceae","Tolypanthus"=>"Loranthaceae",
	"Arceuthobium"=>"Santalales","Daenikera"=>"Santalales","Dendrotrophe"=>"Santalales","Viscum"=>"Santalales","Santalum"=>"Santalales","Opilia"=>"Santalales","Strombosia"=>"Santalales","Malania"=>"Santalales","Hachettea"=>"Santalales","Mystropetalon"=>"Santalales");
my $head = "";
my %fasta = ();
my %region = ();

while (<FAS>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}
while (<IN>) {
	chomp;
	my @tmp = split /\t/,$_;
	my $string = substr($fasta{$tmp[0]},$tmp[6]-1,$tmp[7]-$tmp[6]+1);
	my $species = get_species($tmp[1]);
	my @cover = $string =~ /N/g;
	if (scalar @cover eq 0) {
		substr($fasta{$tmp[0]},$tmp[6]-1,$tmp[7]-$tmp[6]+1) = "N" x ($tmp[7]-$tmp[6]+1);
		push @{$region{$tmp[0]}}, $tmp[2]."-".$tmp[6]."-".$tmp[7]."-".$species;
	}
	else {
		my @last = split ("-",$region{$tmp[0]}[-1]);
		if ($last[0] < $tmp[2]) {
			if ($tmp[6]-$last[1] < 200) {$region{$tmp[0]}[-1] = ""}
			else {$region{$tmp[0]}[-1] = $last[0]."-".$last[1]."-".($tmp[6]-1)."-".$last[3]}
			push @{$region{$tmp[0]}}, $tmp[2]."-".$tmp[6]."-".$tmp[7]."-".$species;
			substr($fasta{$tmp[0]},$tmp[6]-1,$tmp[7]-$tmp[6]+1) = "N" x ($tmp[7]-$tmp[6]+1);
		}
		else {
			if ($tmp[7]-$last[2] > 200) {
				push @{$region{$tmp[0]}}, $tmp[2]."-".($last[2]+1)."-".$tmp[7]."-".$species;
				substr($fasta{$tmp[0]},$tmp[6]-1,$tmp[7]-$tmp[6]+1) = "N" x ($tmp[7]-$tmp[6]+1);
			}
		}
	}
}

foreach my $i (sort keys %region) {
	foreach my $hit (@{$region{$i}}) {
		if ($hit ne "") {
			my @hit = split /-/,$hit;
			print $hit[3],"\t",$hit[0],"\t",$hit[1],"\t",$hit[2];
			if (exists $clade{$hit[3]}) {print "\t",$clade{$hit[3]},"\n"}
			else {print "\n"}
			
#			print "$i\t$hit\n";
		}
	}
}

sub get_species {
	my $seq = shift;
	if ($seq =~ /^(\S+?)_/) {return "$1"}
}