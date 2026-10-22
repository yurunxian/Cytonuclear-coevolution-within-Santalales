use strict;
#use Modern::Perl;

my @species = ("Amb","Ari","Ara","Ery","Mal","Rph","Bya","Hac","Mco","Tma","Scu","Tax","Cha","The","Scl","Pyr","Ene","Ecu","Osy","San","Dae","Den","Pri","Pto","Kor","Ach","Api","Val","Vov");
my $orthogroup = $ARGV[0];
if ($orthogroup eq "") {die()}

my $outdir = "/data/yurx/work/Amphorogynaceae/Tree_0414/Tree";
my $orthogroup_dir = "Orthogroups.4sp.tsv";
my $seq_number = 0;
my $align_length = 0;
my %total_gene = ();

open (OG,"<$orthogroup_dir") or die();

while (<OG>) {
	chomp;
	if (/Orthogroup/) {next}
	else {
		my @tmp = split /\t/,$_;
		for (my $i = 1; $i < @tmp; $i++) {
			$tmp[$i] =~ s/\s+//g;
			foreach my $gene (split /,/,$tmp[$i]) {push @{$total_gene{$tmp[0]}},$gene}
		}
	}
}

#prepare the sequences
#system("rm -rf $outdir/$orthogroup/");

if (! -e "$outdir/$orthogroup/") {
	system("mkdir $outdir/$orthogroup/");
	open (GENE,">$outdir/$orthogroup/$orthogroup.raw.head") or die();
	open (HEAD,">$outdir/$orthogroup/$orthogroup.head");
	print GENE join("\n",@{$total_gene{$orthogroup}}),"\n";
	close GENE;
	system("seqtk subseq all.4sp.cds $outdir/$orthogroup/$orthogroup.raw.head > $outdir/$orthogroup/$orthogroup.raw.cds.fa");
	system("seqkit translate $outdir/$orthogroup/$orthogroup.raw.cds.fa > $outdir/$orthogroup/$orthogroup.raw.pep.fa");
	system("cd-hit -i $outdir/$orthogroup/$orthogroup.raw.pep.fa -c 0.9 -o $outdir/$orthogroup/$orthogroup.nr.pep.fa");
	open (CDS,"<$outdir/$orthogroup/$orthogroup.raw.cds.fa") or die();
	open (CDSO,">$outdir/$orthogroup/$orthogroup.cds.fa") or die();
	my %cds = ();
	my $head = "";
	while (<CDS>) {
		chomp;
		if (/>(\S+)/) {$head = $1; print HEAD $head,"\n"} else {$cds{$head} .= $_}
	}
	foreach my $seq (sort keys %cds) {
		for (my $i = 0; $i < length $cds{$seq}; $i += 3) {
			my $string = substr($cds{$seq},$i,3);
			if ($string =~ /TAG/i or $string =~ /TAA/i or $string =~ /TGA/i or $string =~ /N/i) {
				$cds{$seq} = substr($cds{$seq},0,$i);
				last;
			}
		}
	}
	foreach my $species (@species) {foreach my $key (grep {$_ =~ /^$species/} keys %cds) {print CDSO ">",$key,"\n",$cds{$key},"\n"; $seq_number++}}
	system("rm $outdir/$orthogroup/$orthogroup.raw.head");
	system("rm $outdir/$orthogroup/$orthogroup.raw.cds.fa");
	system("rm $outdir/$orthogroup/$orthogroup.raw.pep.fa");
	system("rm $outdir/$orthogroup/$orthogroup.nr.pep.fa");
	system("rm $outdir/$orthogroup/$orthogroup.nr.pep.fa.clstr");
	close PEP;
	close CDS;
	close CDSO;
	close GENE;
}

#align the sequences
if (! -s "$outdir/$orthogroup/$orthogroup.pep.fas") {
	system("seqkit translate $outdir/$orthogroup/$orthogroup.cds.fa > $outdir/$orthogroup/$orthogroup.pep.fa");
	system("mafft --thread 8 --maxiterate 1000 --quiet --localpair $outdir/$orthogroup/$orthogroup.pep.fa > $outdir/$orthogroup/$orthogroup.pep.fas");
	system("trimal -in $outdir/$orthogroup/$orthogroup.pep.fas -backtrans $outdir/$orthogroup/$orthogroup.cds.fa -ignorestopcodon -fasta -gt 0.3 -st 0.001> $outdir/$orthogroup/$orthogroup.cds.fa.tmp");
	system("sed -i \"s/-//g\" $outdir/$orthogroup/$orthogroup.cds.fa.tmp");
	system("seqkit translate $outdir/$orthogroup/$orthogroup.cds.fa.tmp > $outdir/$orthogroup/$orthogroup.pep.fa.tmp");
	system("mafft --thread 8 --maxiterate 1000 --quiet --localpair $outdir/$orthogroup/$orthogroup.pep.fa.tmp > $outdir/$orthogroup/$orthogroup.pep.fas.tmp");
	system("pal2nal.pl $outdir/$orthogroup/$orthogroup.pep.fas.tmp $outdir/$orthogroup/$orthogroup.cds.fa.tmp -output fasta -nomismatch -nostderr > $outdir/$orthogroup/$orthogroup.cds.fas");
	open(GB,"<$outdir/$orthogroup/$orthogroup.cds.fas") or die();
	open(GBO,">$outdir/$orthogroup/$orthogroup.trim.cds.fas") or die();
	my %fasta_trimal;
	my @order = ();
	my $head = "";
	while (<GB>) {
		chomp;
		if (/>(\S+)/) {$head = $1; push @order,$head} else {$fasta_trimal{$head} .= $_}
	}
	foreach my $key (@order){
		$align_length = length($fasta_trimal{$key});
		my $gap = 0;
		if (my @tmp = $fasta_trimal{$key} =~ /-/gi) {$gap = scalar(@tmp)}
		if ($gap/$align_length <= 0.7) {print GBO ">",$key,"\n",$fasta_trimal{$key},"\n"}
	}
	system("rm $outdir/$orthogroup/*tmp*");
	system("rm $outdir/$orthogroup/$orthogroup.cds.fas");
	system("rm $outdir/$orthogroup/$orthogroup.pep.fa");
	close GB;
	close GBO;
}
#if ($align_length <= 150) {die("Ops! The $orthogroup alignment is too short: $align_length.\n")}

#construct the phylogeny tree with RAxML-NG
if (! -s "$outdir/$orthogroup/$orthogroup.raw.tree") {
	system("~/software/raxml-ng -all --msa $outdir/$orthogroup/$orthogroup.trim.cds.fas --model GTR+G --tree pars{10} --bs-trees 100 --threads 8 --redo");
	system("mv $outdir/$orthogroup/$orthogroup.trim.cds.fas.raxml.support $outdir/$orthogroup/$orthogroup.raw.tree");
	system("rm $outdir/$orthogroup/*raxml*");
}