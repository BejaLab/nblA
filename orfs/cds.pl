#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my ($file) = @ARGV;

my $in = Bio::SeqIO->new(-file => $file, -format => 'embl') or die "$!\n";

my ($proteome) = $file =~ /(UP\d+)/;

while (my $chr = $in->next_seq) {
	my $species = $chr->species;
	my ($taxid) = $species->ncbi_taxid;
	my $classification = join ';', reverse $species->classification;
	$classification =~ tr/ /_/;
	my $len = length $chr->seq;
	my %cds = ();
	for my $feature ($chr->get_SeqFeatures) {
		next if $feature->primary_tag ne "CDS";
		my ($translation) = tag($feature, 'translation');
		$translation = $feature->seq->translate(-codontable_id => 2) if not $translation;

		my ($locus_tag)  = tag($feature, 'locus_tag');
		my ($protein_id) = tag($feature, 'protein_id');
		my ($product)    = tag($feature, 'product');
		$product =~ tr/ /_/;

		my $xref = join ';', tag($feature, 'db_xref');
		my ($db, $kb) = ($xref =~ /UniProtKB\/(Swiss-Prot|TrEMBL):(\w+)/);

		my $id = $kb;
		$id = $protein_id if not $id;
		$id = $locus_tag  if not $id;
		my $start = $feature->strand > 0 ? $feature->start : $feature->end;
		my $end   = $feature->strand > 0 ? $feature->end   : $feature->start;

		my $desc = sprintf "[%s - %s] dna=%s taxid=%s proteome=%s classification=%s", $start, $end, $chr->id, $taxid, $proteome, $classification;
		$desc .= sprintf " xref=%s", $xref if $xref;
		$desc .= sprintf " product=%s", $product if $product;

		printf ">%s_%s %s\n%s\n", $chr->id, $id, $desc, $translation;
	}
}

sub tag {
	my $feature = shift;
	my $tag = shift;
	return ("") if not $feature->has_tag($tag);
	return $feature->get_tag_values($tag);
}
