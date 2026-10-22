use strict;

open (IN,"<OXPHOS.txt") or die();
open (ARA,"</data/yurx/work/Amphorogynaceae/Database/cds_1204/Arabidopsis.pep");

my $head = "";
my $outdir = "./Family";
my $datadir = "./Data";
my @species = ();
my %arabidopsis = ();
my %all_cds = ();

#system("cat $datadir/*.cds > $datadir/all_cds.fa");
foreach my $i (glob "$datadir/*") {
	if ($i =~ /$datadir\/(\S+).pep/) {
		my $sp = $1;
		push @species, $sp;
		if (! -s "$datadir/$sp.pin") {system("makeblastdb -in $datadir/$sp.pep -dbtype prot -out $datadir/$sp")}
	}
}

print join("\n",@species),"\n";

while (<ARA>) {
	$_ =~ s/[\r\n]+//g;
	if (/>(\S+)/) {$head = $1; $head =~ s/\.\d//g}
	else {$arabidopsis{$head} .= $_}
}

while (<IN>) {
	if (/#/) {next}
	$_ =~ s/[\r\n]+//g;
	my @tmp = split /\t/,$_;
	my $name = $tmp[0];
	$name =~ s/\s//g;
	$tmp[1] =~ s/\s//g;
	$tmp[1] =~ s/\.\d//g;
	my @ara_gene = split /,/,$tmp[1];

	# create the analysis folder
	if (! -e "$outdir/$name") {system("mkdir $outdir/$name")}
	if (! -s "$outdir/$name/$name.ara.fasta") {
		open (OUT,">$outdir/$name/$name.ara.fasta");
		foreach my $gene (sort @ara_gene) {print OUT ">$gene\n",$arabidopsis{$gene},"\n"}
		close OUT;
	}
#	if (! -s "$outdir/$name/$name.pin") {system("makeblastdb -in $outdir/$name/$name.ara.fasta -dbtype prot -out $outdir/$name/$name")}

	# execute blastp searching
	if (! -s "$outdir/$name/$name.candidate") {
		system("mkdir $outdir/$name/blastp");
		open (OUT,">$outdir/$name/$name.candidate") or die();
		foreach my $sp (@species) {
			if (! -e "$outdir/$name/blastp/$sp.blastp") {
				system("blastp -query $outdir/$name/$name.ara.fasta -db $datadir/$sp -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $outdir/$name/blastp/$sp.blastp");
			}
			open (BLASTP,"<$outdir/$name/blastp/$sp.blastp") or die();
			my %tmp = ();
			while (<BLASTP>) {
				if (!/\S/) {next}
				my @tmp = split /\t/,$_;
				if ($tmp[2] > 35 and $tmp[3] > length{$arabidopsis{$tmp[0]}} * 0.3) {$tmp{$tmp[1]}++}
			}
			foreach my $i (sort keys %tmp) {print OUT "$i\t$sp\n"}
		}
		close OUT;
	}

	# get the CDS sequences
	if (! -s "$outdir/$name/$name.fa") {
		system("seqtk subseq $datadir/all_cds.fa $outdir/$name/$name.candidate > $outdir/$name/$name.tmp");
		open (CDS,"<$outdir/$name/$name.tmp") or die();
		open (OUT,">$outdir/$name/$name.fa") or die();
		my %CDS = ();
		my %length = ();
		while (<CDS>) {
			chomp;
			if (/>(\S+)/) {$head = $1}
			else {$CDS{$head} .= $_}
		}
		foreach my $i (sort keys %CDS) {
			if ($i =~ /(\S+_c\d+_g\d+)_i\d+\.p\d+/) {push @{$length{$1}}, $i}
			else {push @{$length{$i}}, $i}
		}
		foreach my $i (sort keys %length) {
			if (scalar @{$length{$i}} eq 1) {print OUT ">",$length{$i}[0],"\n",$CDS{$length{$i}[0]},"\n"}
			else{
				my @order = sort {length $CDS{$b} <=> length $CDS{$a}} @{$length{$i}};
				print OUT ">",$order[0],"\n",$CDS{$order[0]},"\n";
			}
		}
		system("rm $outdir/$name/$name.tmp");
		close CDS;
		close OUT;
	}

	# align the CDS sequences
	if (! -s "$outdir/$name/$name.trim.fas") {
		system("seqkit translate -x $outdir/$name/$name.fa > $outdir/$name/$name.pep.fa");
		system("mafft --thread 8 --quiet --maxiterate 1000 --localpair $outdir/$name/$name.pep.fa > $outdir/$name/$name.pep.fas");
		system("perl pal2nal.pl $outdir/$name/$name.pep.fas $outdir/$name/$name.fa -output fasta > $outdir/$name/$name.cds.fas");
		system("trimal -in $outdir/$name/$name.cds.fas -fasta -gt 0.3 -out $outdir/$name/$name.trim.fas");
		system("rm $outdir/$name/$name.pep.fa*");
		system("rm $outdir/$name/$name.fa");
	}

	# phylogenetic inference
	if (! -s "$outdir/$name/$name.tree") {
		system("~/software/raxml-ng -all --msa $outdir/$name/$name.trim.fas --model GTR+G4 --tree pars{10} --bs-trees 100 --threads 4 --redo");
		system("mv $outdir/$name/$name.trim.fas.raxml.support $outdir/$name/$name.tree");
		system("rm $outdir/$name/*raxml*");
	}
}
