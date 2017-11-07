# Minor-Variant-Phasing-Tool
Minor-Variant-Phasing-Tool for Pacbio data

### Input
* Bam file: Bam file should be generated by Pacbio blasr program with CCS reads.  
_The CCS sequences for alignment: 10 full pass, barcode quality of >= 45, predicted read accuracy >= 0.99._
* Reference file: Fasta format.
* SNPs file: csv records, with first column records _segment name_, second for _position_, third for _Ref nucleotide_ and the fourth for _Alt nucleotide_. Like below:

	|CHROM|POS|REF|ALT| 
	|:--|:--|:--|:--|
	|Brisbane_PB1_codon_mutant|4|G|A|
	|Brisbane_M_codon_mutant_|92|G|C|
	|brisbane_NS_codon_mutant|100|A|G|

### Output
* Phasing graph in jpeg format.

------
Copyright (c) 2017 Haogao Gu. All rights reserved.