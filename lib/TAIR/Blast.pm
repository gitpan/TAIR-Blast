package TAIR::Blast;

use warnings;
use strict;

=head1 NAME

TAIR::Blast - A module to gather automated BLAST result from TAIR (http://www.arabidopsis.org/Blast/index.jsp) 

=head1 VERSION

Version 1.00

=cut

our $VERSION = '1.00';


=head1 SYNOPSIS

This module simply automatically! BLAST any type of sequences (nucleotide, protein) with using different type of algorithm (blastp, blastn, tblastx etc.) by using TAIR Blast engine.

    #METHOD #1 for one sequence only
    use TAIR::Blast;

    my $TB = TAIR::Blast->new();
    
    #START OF PARAMETERS FOR BLAST#
	my $query = "estExt_fgenesh4_pm.C_1450005";
	my $query_seq = "ATGGGAGCTGTTGCTGCGATATGGTGGTTAACGGTGGTTGCTGCTGCATCTCGCCTTTCATTCCTGCACGCCT";
	my $algorithm = "blastn"; #other set : blastp, blastx, tblastx,tblastn
	my $maxscore; #http://www.arabidopsis.org/help/helppages/BLAST_help.jsp#alignments (OPTIONAL #default=250)
	my $blast_target_set; #http://www.arabidopsis.org/help/helppages/BLAST_help.jsp#datasets (OPTIONAL #default for blastn, tblastx, tblastn is At_transcripts and for blastp and blastx is ATH1_pep. You can alternatively select from the list according to algorithm)
	#At_transcripts = TAIR10 Transcripts (-introns, +UTRs) (DNA) 
	#ATH1_cds = TAIR10 CDS (-introns, -UTRs) (DNA) 
	#ATH1_seq = TAIR10 Genes (+introns, +UTRs) (DNA) 
	#ATH1_pep = TAIR10 Proteins (Protein) 
	#ATH1_bacs_con = TAIR10 Whole Genome (BAC clones) (DNA) 
	#At_upstream_500 = TAIR10 Loci Upstream Sequences -- 500 bp (DNA) 
	#At_upstream_1000 = TAIR10 Loci Upstream Sequences -- 1000 bp (DNA) 
	#At_upstream_3000 = TAIR10 Loci Upstream Sequences -- 3000 bp (DNA) 
	#At_downstream_1000 = TAIR10 Loci Downstream Sequences -- 1000 bp (DNA) 
	#At_downstream_3000 = TAIR10 Loci Downstream Sequences -- 3000 bp (DNA) 
	#At_intergenic = TAIR10 Intergenic (DNA) 
	#At_intron = TAIR10 Intron (DNA) 
	#ATH1_5_UTR = TAIR10 5' UTRs (DNA) 
	#ATH1_3_UTR = TAIR10 3' UTRs (DNA) 
	#At_flanks_DNA = A. thaliana Insertion Flanks (DNA) 
	#At_Uniprot_prot = A. thaliana UniProt (Protein) 
	#At_GB_exp_prot = A. thaliana GB derived from mRNA (Protein) 
	#At_GB_refseq_prot = A. thaliana GB refseq/tpa (Protein) 
	#At_GB_all_prot = A. thaliana GB all (Protein) 
	#At_GB_exp_tx_DNA = A. thaliana GB experimental cDNA/EST (DNA) 
	#At_GB_refseq_tx_DNA = A. thaliana GB refseq/tpa cDNA (DNA) 
	#At_GB_genomic_DNA = A. thaliana GB genomic (DNA) 
	#gp_GB_exp_prot = Green plant GB derived from mRNA (Protein) 
	#gp_GB_refseq_prot = Green plant GB refseq/tpa (Protein) 
	#gp_GB_all_prot = Green plant GB all (Protein) 
	#gp_GB_exp_tx_DNA = Green plant GB experimental cDNA/EST (DNA) 
	#gp_GB_refseq_tx_DNA = Green plant GB refseq/tpa cDNA (DNA) 
	#gp_GB_genomic_DNA = Green plant GB genomic (DNA) 
	#END OF PARAMETERS FOR BLAST#
	my $verbose = 1; #default (0, no info);
	
    my $result = $TB->connect($query,$query_seq,$algorithm,$maxscore,$blast_target_set,$verbose);
    
    #you would like to have the output recorded!
    my $output_file = "output_file.txt";
    $TB->write($output_file);
    
    
    #METHOD #2 for multiple sequences with the help of Bio::DB::Fasta 
    use TAIR::Blast;
	use Bio::DB::Fasta; #You may use different but it is good for parsing your fasta file!
	
	my $TB = TAIR::Blast->new();
	
	my $fasta_file = "seqs.fasta";
	
	my $algorithm = "blastn"; #other set : blastp, blastx, tblastx,tblastn
	my $maxscore; #http://www.arabidopsis.org/help/helppages/BLAST_help.jsp#alignments (OPTIONAL #default=250)
	my $blast_target_set;
	my $verbose = 1;
	
	my $stream  = Bio::DB::Fasta->new($fasta_file)->get_PrimarySeq_stream;
	  while (my $query = $stream->next_seq) {
		my $query_seq = $query->seq;
		#OTHER PARAMETERS ARE SAME AS IN METHOD 2
		my $result = $TB->connect($query,$query_seq,$algorithm,$maxscore,$blast_target_set,$verbose);
	  }
	my $output_file = "output_file.txt";
	$TB->write($output_file);
	
	
    ...

=head1 SUBROUTINES/METHODS

=head2 new

Object-oriented master-hash!

=cut

sub new
  {
    my %master_hash;
    return(bless(\%master_hash, __PACKAGE__));
  }
  
=head2 connect

Blast parameters, do not mess with them unless you know what you are doing.

=cut

sub connect
 {
	my $tair_master = shift;
	my $query = shift;
	my $query_seq = shift;
	my $algorithm = shift;
	#my $out_type = shift; #will be used in v2
	my $maxscore = shift;
	my $blast_target_set = shift;
	my $verbose = shift;
		
	use LWP::UserAgent;
	use HTTP::Request::Common;
	
	if (!defined($blast_target_set)) {
		if ($algorithm =~ /^blast[px]$/i){
			$blast_target_set = "ATH1_pep";
		}elsif ($algorithm =~ /^blastn|tblastx|tblastn$/i){
			$blast_target_set = "At_transcripts"; 	
		}
	}
			
	my $URLtoPostTo = "http://www.arabidopsis.org/cgi-bin/Blast/TAIRblast.pl";
	my %Fields = (
   "Algorithm" => $algorithm,
   "default_db" => "At_transcripts",
   "BlastTargetSet" => $blast_target_set,
   "textbox" => "seq",
   "QueryText" => $query_seq,
	"QueryFilter" => "T",
	"Matrix" => "Blosum62",
	"MaxScore" => "3",
	"Expectation" => "10",
	"MaxAlignments" => $maxscore,
	"NucleicMismatch" => "-3",
	"GappedAlignment" => "T",
	"NucleicMatch" => "1",
	"OpenPenalty" => "0",
	"ExtensionThreshold" => "0",
	"ExtendPenalty" => "0",
	"WordSize" => "0",
	"QueryGeneticCode" => "1",
	"Comment" => "optional, will be added to output for your use",
	"ReplyTo" => "",
	"ReplyVia" => "BROWSER",
	#"ReplyFormat" => $out_type, #will be used in v2
	"ReplyFormat" => "TABULAR",
	"PageType" => "JavaScr"
	);

	my $BrowserName = "TAIR";

	my $Browser = new LWP::UserAgent;

	if($BrowserName) { $Browser->agent($BrowserName); }

	my $Page = $Browser->request(POST $URLtoPostTo,\%Fields);

	if ($Page->is_success) {	
	my $output_raw = $Page->content;
	$output_raw =~ s/TAIR Blast Job Pending//g;
	my @output_array = split("\n", $output_raw); #[0]query/[1]hit/[10]e-value/[11]Score(bits)
	my %output_hsps;
	my @output_hsps;
	my $i = 1;	
		foreach (@output_array){
			my @details = split("\t",$_);
				if (defined($details[1])) {
					my %details = (
					"hsp" => $i,
					"hits" => $details[1],
					"e-value" => $details[10],
					"score" => $details[11]
					);
					shift(@details);
					push (@output_hsps, \%details);	
					$i++;
				}
		}	
	$tair_master->{'output'}->{$query}=\@output_hsps;
	if ($verbose == 1){
		print $query, " is blasted!\n";	
	}
 	}
	else { print $Page->message; }
	
 }

=head2 write

This subroutine handles output of the program, it writes the blast results into a specified file name.

=cut

sub write 
	{
		my $tair_master = shift;
		my $output_file = shift;
		my ($output_file_hn);
		my $output_hash = $tair_master->{'output'};
		open($output_file_hn, ">", $output_file) or die("Cannot open $output_file for writing: $!");
		print $output_file_hn "query\thits\tscore\te-value\n";	
			foreach (keys %$output_hash){		
			my $next_query = $_; 
				foreach(@{$output_hash->{$_}}){
				print $output_file_hn $next_query,"\t",
						$_->{'hits'},"\t",
						$_->{'score'},"\t",
						$_->{'e-value'},"\n";				
				}
			}
	}

=head1 AUTHOR

Haktan Suren, << <hsuren at vt.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-tair-blast at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=TAIR-Blast>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc TAIR::Blast


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=TAIR-Blast>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/TAIR-Blast>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/TAIR-Blast>

=item * Search CPAN

L<http://search.cpan.org/dist/TAIR-Blast/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Haktan Suren.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of TAIR::Blast
