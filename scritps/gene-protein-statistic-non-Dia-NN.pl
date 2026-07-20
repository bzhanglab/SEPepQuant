#!/usr/bin/env perl
use warnings;

$protein_database_file=$ARGV[0];
$output_dir=$ARGV[1];

%gene_protein_number=();
open(IN,$protein_database_file);
while($raw=<IN>)
{
	if($raw=~/^>/&&$raw!~/^>rev_/&&$raw!~/^>Cont/)
	{
		@inf=split(/\s|\|/,$raw);
		$gene_protein_number{$inf[2]}++;
	}	
}
close(IN);

open(OUT,">$output_dir/gene-protein-statistic.txt");

foreach $key (sort keys(%gene_protein_number))
{
	print OUT "$key\t$gene_protein_number{$key}\n";
}

close(OUT);















