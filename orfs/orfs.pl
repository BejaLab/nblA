#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-fh => \*STDIN, -format => 'embl') or die "$!\n";

my $min = 50;
my $max = 108;
my %starts = ( ATG => 1, GTG => 2, TTG => 3 );

my %code = (
	TTT => 'F',  TCT => 'S',  TAT => 'Y',  TGT => 'C',
	TTC => 'F',  TCC => 'S',  TAC => 'Y',  TGC => 'C',
	TTA => 'L',  TCA => 'S',  TAA => '*',  TGA => '*',
	TTG => 'L',  TCG => 'S',  TAG => '*',  TGG => 'W',
	CTT => 'L',  CCT => 'P',  CAT => 'H',  CGT => 'R',
	CTC => 'L',  CCC => 'P',  CAC => 'H',  CGC => 'R',
	CTA => 'L',  CCA => 'P',  CAA => 'Q',  CGA => 'R',
	CTG => 'L',  CCG => 'P',  CAG => 'Q',  CGG => 'R',
	ATT => 'I',  ACT => 'T',  AAT => 'N',  AGT => 'S',
	ATC => 'I',  ACC => 'T',  AAC => 'N',  AGC => 'S',
	ATA => 'I',  ACA => 'T',  AAA => 'K',  AGA => 'R',
	ATG => 'M',  ACG => 'T',  AAG => 'K',  AGG => 'R',
	GTT => 'V',  GCT => 'A',  GAT => 'D',  GGT => 'G',
	GTC => 'V',  GCC => 'A',  GAC => 'D',  GGC => 'G',
	GTA => 'V',  GCA => 'A',  GAA => 'E',  GGA => 'G',
	GTG => 'V',  GCG => 'A',  GAG => 'E',  GGG => 'G',
);

my @frames = (1, 2, 3, -1, -2, -3);
my %ipr = (IPR036904 => 1, IPR007574 => 1);

while (my $chr = $in->next_seq) {
	my $len = length $chr->seq;
	my %cds = ();
	for my $feature ($chr->get_SeqFeatures) {
		next if $feature->primary_tag ne "CDS";
		my $kb;
		my @refs = $feature->get_tag_values('db_xref');
		my $nbla = 0;
		my $other = 0;
		for my $db_xref (@refs) {
			my ($id) = $db_xref =~ /:(\w+)/;
			$kb = $id if $db_xref =~ /UniProtKB/;
			if ($db_xref =~ /InterPro/) {
				$nbla  = 1 if $ipr{$id};
				$other = 1 if not $ipr{$id}; 
			}
		}
		die "No uniprot id\n" if not $kb;
		my ($translation) = $feature->get_tag_values('translation');
		my $offset = $feature->strand > 0 ? $feature->start - 1 : $len - $feature->end;
		my $frame = $offset % 3 + 1;
		%{$cds{$frame}{$offset}} = (
			kb     => $kb,
			start  => $feature->start,
			end    => $feature->end,
			strand => $feature->strand,
			nbla   => $nbla,
			other  => $other,
			translation => $translation
		);
	}

	for my $frame (@frames) {
		my $seq = $chr->seq;
		if ($frame < 0) {
			$seq =~ tr/ATGC/TACG/;
			$seq = reverse $seq;
		}
		$seq = substr $seq, abs($frame) - 1;
		my @allCodons = ($seq =~ /.{3}/g);
		my @stops  = grep { $allCodons[$_] =~ /TAG|TGA|TAA/ } 0 .. $#allCodons;
		my $lastStop = -1;
		for my $stop (@stops) {
			my @codons = reverse @allCodons[$lastStop + 1 .. $stop - 1];
			$lastStop = $stop;

			my $long = 0;
			my @pept = ();

			my $pos = 0;
			for my $codon (@codons) {
				my $aa = $code{$codon} ? $code{$codon} : "X";
				push @pept, $code{$codon};
				$long = $pos if $starts{$codon} and $pos >= $min;
				last if $pos == $max;
				$pos++;
			}
			$long = $pos-1 if not $long and $pos <= $max and $pos >= $min;
			if ($long) {
				@pept = reverse(@pept[0 .. $long]);
				printf ">%s frame=%s codons=[%d-%d]\n%s\n", $chr->id, $frame, $stop-$long, $stop-1, join('', @pept);
			}
		}
	}
}
