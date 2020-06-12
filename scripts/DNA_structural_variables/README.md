# DNA structural variables

Tools for prediction of DNA structural properties used in the following publications:
* [Zrimec & Lapanje 2018: DNA structure at the plasmid origin-of-transfer indicates its potential transfer range](https://www.nature.com/articles/s41598-018-20157-y)
* [Tosato et al. 2017: Bridge-Induced Translocation between NUP145 and TOP2 Yeast Genes Models the Genetic Fusion between the Human Orthologs Associated With Acute Myeloid Leukemia](https://www.frontiersin.org/articles/10.3389/fonc.2017.00231/full)

## Description
Conformational and physicochemical structural properties of the DNA double helix are important for cell homeostasis as well as bacterial pathogenesis, since they are responsible for regulation of key genetic processes by assisting protein-DNA recognition and binding. The processes include DNA transcription, replication and horizontal transfer of mobile genetic elements that carry virulence and antimicrobial resistance genes. Common to the initiation of these processes is the occurrence of DNA melting bubbles, which form as a consequence of intrinsically low dsDNA stability and extrinsic duplex destabilization in neighboring DNA induced by proteins, topological DNA features and thermal fluctuations (thermally induced duplex destabilization, TIDD, Zrimec & Lapanje 2015). According to our previous results oriT conformational and physicochemical properties were more informative than nucleotide sequences for descriminating group of DNA substrates (see figure below from Zrimec & Lapanje 2018).

<img src="https://github.com/JanZrimec/DNA_structural_variables/blob/master/Figure_3.jpg" width="480">

Here we gathered besides the 6 structural properties that were shown to be able to distinguish between different DNA substrates additional DNA structure models shown to be informative for analysis of DNA-protein interactions. These models were based on nearest neighbor dinucleotide (54) and trinucleotide (4) models and included physicochemical and conformational properties and properties attributed to DNA-protein interactions (see List_structural_variables.csv).

## Usage
```out = get_structures_par(seqs,window)```

where:
* seqs .. fasta file
* window .. size of sliding window
