#!/usr/bin/env perl
use warnings;

$protein_database_file=$ARGV[0];
$output_dir=$ARGV[1];

open(OUT,">$output_dir/protein-gene-mapping.txt");
%gene_protein_number=();
open(IN,$protein_database_file);
while($raw=<IN>)
{
	if($raw=~/^>/&&$raw!~/^>rev_/&&$raw!~/^>Cont/)
	{
		@inf=split(/\s+|\|/,$raw);
		$gene_name=$inf[@inf-3];
		$gene_name=~s/GN=//;
		$gene_protein_number{$gene_name}++;
		print OUT "$inf[1]\t$gene_name\n";
	}	
}
close(IN);
close(OUT);


open(OUT,">$output_dir/gene-protein-statistic.txt");

foreach $key (sort keys(%gene_protein_number))
{
	print OUT "$key\t$gene_protein_number{$key}\n";
}

close(OUT);















