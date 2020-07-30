#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Tree::Interval::Fast;

my ($f1, $f2, $proteome, $min, $max) = @ARGV;

my $cds = Bio::SeqIO->new(-file => $f1, -format => 'fasta') or die "$!\n";
my $orf = Bio::SeqIO->new(-file => $f2, -format => 'fasta') or die "$!\n";
my $out = Bio::SeqIO->new(-fh   => \*STDOUT, -format => 'fasta') or die "$!\n";

my %intervals;
my $taxid = 0;
my $classification = "";

my %ids;

while (my $seq = $cds->next_seq) {
	my ($dna, $id) = split '_', $seq->id;

	while ($ids{$id}) {
		$id = $id . "X";
	}
	$ids{$id} = 1;

	my ($start, $end) = ($seq->description =~ /\[(\d+) - (\d+)\]/);
	my ($xref) = ($seq->description =~ /xref=([^ ]*)/);
	my ($product) = ($seq->description =~ /product=([^ ]+)/);

	($taxid) = ($seq->description =~ /taxid=(\d+)/) if not $taxid;
	($classification) = ($seq->description =~ /classification=([^ ]+)/) if not $classification;

	$intervals{$dna} = Tree::Interval::Fast->new() if not $intervals{$dna};
	my $strand = $start < $end;
	my $low  = $strand ? $start : $end;
	my $high = $strand ? $end : $start;

	my $interval = Tree::Interval::Fast::Interval->new($low, $high, $strand);
	$intervals{$dna}->insert($interval);
	next if $max and length $seq->seq > $max;

	$seq->id($id);
	my $desc = sprintf "range=%s:%d-%d taxid=%s proteome=%s classification=%s", $dna, $start, $end, $taxid, $proteome, $classification;
	$desc .= sprintf " xref=%s", $xref if $xref;
	$desc .= sprintf " product=%s", $product if $product;
	$seq->description($desc);
	$out->write_seq($seq);

}

SEQ: while (my $seq = $orf->next_seq) {
	my ($dna) = split '_', $seq->id;
	# ($dna) = split '.', $dna;
	my ($start, $end) = ($seq->description =~ /\[(\d+) - (\d+)\]/);

	my $len = length $seq->seq;
	next if $max and $len > $max;
	next if $len < $min;
	my $strand = $start < $end;
	my $low  = $strand ? $start : $end;
	my $high = $strand ? $end : $start;

	if ($intervals{$dna}) {
		my $ivals = $intervals{$dna}->findall($low, $high);
		if ($ivals) {
			for my $interval (@{$ivals}) {
				next if $interval->data != $strand;
				next SEQ if $interval->low % 3 == $low % 3;
			}
		}
	}

	my $desc = sprintf "range=%s:%d-%d taxid=%s proteome=%s classification=%s", $dna, $start, $end, $taxid, $proteome, $classification;

	$seq->id(sprintf "%s_%s_%s", $dna, $start, $end);
	$seq->description($desc);
	$out->write_seq($seq);
}
