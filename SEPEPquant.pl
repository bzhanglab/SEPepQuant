#!/usr/bin/env perl

use Getopt::Long qw(:config no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use FindBin;
use lib $FindBin::RealBin;

use Env qw(@PATH);

@PATH = ($FindBin::RealBin, @PATH);

use strict;
use warnings;
my @folder_list=();

my $SEPEP_fdr=0.01;
my $help=0;
my $protein_database="";
my $quant_type="";
my $tmt_plex="";
my $Ref_Tag="";
my $output_dir="";
my $input_dir="";

if (scalar(@ARGV) == 0)
{
    pod2usage(-msg => "Usage :\n--database: protein database used for database searching\n--quant: LF or TMT\n--plex: TMT plex, required if --quant is TMT\n--RefTag: Tag of TMT reference channel, required if --quant is TMT\n--input: folder of FragPipe output\n--output: output filder\n--help: help information\n",-exitval => 2)
}

GetOptions("database=s" => \$protein_database,
    "fdr=f" => \$SEPEP_fdr,
    "quant=s" => \$quant_type,
    "plex=i" => \$tmt_plex,
    "RefTag=s" => \$Ref_Tag,
    "input=s" => \$input_dir,
    "output=s" => \$output_dir,
    "h|help" => \$help);

if ($help ==1 )
{
    pod2usage(-msg => "Usage:\n--database: protein database used for database searching\n--quant: LF or TMT\n--plex: TMT plex, required if --quant is TMT\n--RefTag: Tag of TMT reference channel, required if --quant is TMT\n--input: folder of FragPipe output\n--output: output filder\n--help: this information\n",-exitval => 2)
}
elsif (!$protein_database || !$quant_type || !$input_dir || !$output_dir) 
{
    pod2usage(-msg => "\nAll these parameters are required: --database, --quant, --input, and --output!\n--database: protein database used for database searching\n--quant: LF or TMT\n--input: folder of FragPipe output\n--output: output filder",-exitval => 2);
}

if ($tmt_plex eq "TMT" && (!$tmt_plex || !$Ref_Tag)) 
{
    pod2usage(-msg => "\n--plex and --RefTag are required for TMT data",-exitval => 2);
}

# check files and folders

if (!(-e $protein_database)) 
{
    pod2usage(-msg => "\nCan not find protein database: $protein_database",-exitval => 2);
}

if (!(-e $input_dir)) 
{
    pod2usage(-msg => "\nCan not find input folder: $input_dir",-exitval => 2);
}
else
{
    opendir(DIR, $input_dir)|| die "Cannot open directory: $input_dir!";
    @folder_list = readdir(DIR);
    closedir(DIR);
    my $total_peptide_files=0;
    foreach my $folder_name (@folder_list)
    {
        if($folder_name ne "\."&&$folder_name ne "\.\.")
        {
            if(!(-e "$input_dir/$folder_name/peptide.tsv"))
            {
                print "No peptide.tsv found in $input_dir/$folder_name\n";
                exit;
            }
            else
            {
                $total_peptide_files++;
                print "Peptide file $total_peptide_files: $input_dir/$folder_name/peptide.tsv \n";
            }
        }
    }
    if($total_peptide_files == 0)
    {
        print "The input folder $input_dir is empty.\n";
    }
}

if(-e "$output_dir")
{
    print "\nWarning! The output direction $output_dir is existing\n\n";
}
else
{
    mkdir( $output_dir ) or die "Couldn't create output_dir directory!\n";
}


# process protein database

print("Processing protein database..... \n");
system("perl $PATH[0]/scritps/gene-protein-statistic.pl $protein_database $output_dir");

# process protein database
if($quant_type eq "TMT")
{
    foreach my $folder_name (@folder_list)
    {
        if($folder_name ne "\."&&$folder_name ne "\.\.")
        {
            print "Processing peptide.tsv for $folder_name\n";
            system("Rscript $PATH[0]/scritps/process-matrix-TMT.r $input_dir/$folder_name $tmt_plex $Ref_Tag $output_dir/gene-protein-statistic.txt");
            system("Rscript $PATH[0]/scritps/identify-sepep-TMT.r $input_dir/$folder_name $tmt_plex");
            system("Rscript $PATH[0]/scritps/sepep-level-fdr-TMT.r $input_dir/$folder_name $SEPEP_fdr");
        }
    }
    system("Rscript $PATH[0]/scritps/creat-mapping-table-TMT.r $input_dir $output_dir $tmt_plex"); 
    system("Rscript $PATH[0]/scritps/quantification-TMT.r $input_dir $output_dir $tmt_plex $Ref_Tag"); 
    system("Rscript $PATH[0]/scritps/combine-and-normalizaztion-TMT.r $input_dir $output_dir");
}
elsif($quant_type eq "LF")
{
    foreach my $folder_name (@folder_list)
    {
        if($folder_name ne "\."&&$folder_name ne "\.\.")
        {
            print "Processing peptide.tsv for $folder_name\n";
            system("Rscript $PATH[0]/scritps/process-matrix-LF.r $input_dir/$folder_name $output_dir/gene-protein-statistic.txt");
            system("Rscript $PATH[0]/scritps/identify-sepep-LF.r $input_dir/$folder_name");
            system("Rscript $PATH[0]/scritps/sepep-level-fdr-LF.r $input_dir/$folder_name $SEPEP_fdr");
        }
    }
    system("Rscript $PATH[0]/scritps/creat-mapping-table-LF.r $input_dir $output_dir"); 
    system("Rscript $PATH[0]/scritps/quantification-LF.r $input_dir $output_dir"); 
    system("Rscript $PATH[0]/scritps/combine-quantification-LF.r $input_dir $output_dir");    
}



