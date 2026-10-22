use strict;

open (OG,"<ribosome.txt") or die();
my %total_gene = ();
while (<OG>) {
	chomp;
	if (/Gene name/) {next}
	else {
		my @tmp = split /\t/,$_;
		for (my $i = 1; $i < @tmp; $i++) {
			$tmp[$i] =~ s/\s+//g;
			foreach my $gene (split /,/,$tmp[$i]) {push @{$total_gene{$tmp[0]}},$gene}
		}
	}
}
foreach my $gene (sort keys %total_gene) {
	get_ctl_file($gene);
	system("bash 7_run_codeml.sh $gene");
}

sub get_ctl_file {
	my $gene = shift;
	open (PAML,">/data/yurx/work/Amphorogynaceae/dNdS/ribosome/$gene/codeml.ctl") or die();
print PAML <<END_OF_TEXT;
      seqfile = /data/yurx/work/Amphorogynaceae/dNdS/ribosome/$gene/$gene.cds.trim.fas
     treefile = /data/yurx/work/Amphorogynaceae/dNdS/ribosome/$gene/$gene.root.tree
      outfile = /data/yurx/work/Amphorogynaceae/dNdS/ribosome/$gene/$gene.codeml

        noisy = 1  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a

        model = 1  * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;

        getSE = 0
 RateAncestor = 0

    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
END_OF_TEXT
	close PAML;
}
