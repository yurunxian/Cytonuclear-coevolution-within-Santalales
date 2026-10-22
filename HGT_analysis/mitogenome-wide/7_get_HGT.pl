use strict;
use warnings;

open (IN,"<P-tonkinensis_mito.candidate.final") or die();

while (<IN>) {
	chomp;
	my @tmp = split /\t/,$_;
	my $best_hit = "";
	my %number = ();
	open(BLAST,"<./P-tonkinensis/$tmp[0].blastn.filter") or die();
	while (<BLAST>) {
		chomp;
		my @blast = split /\t/,$_;
		if ($blast[1] =~ /Amborella/) {next}
		if ($best_hit eq "") {$best_hit = $_}
		if ($blast[1] =~ /^(\S+?)_/) {$number{$1}++}
	}
	print $tmp[0],"\t",scalar keys %number,"\t",$best_hit,"\t",join(", ",sort keys %number),"\n";
}