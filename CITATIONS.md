# bbglab/deepUMIcaller: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

- [Bwa](https://doi.org/10.48550/arXiv.1303.3997)

  > Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv. 2013 May 26. doi: 10.48550/arXiv.1303.3997

- [bwa-mem3](https://github.com/fg-labs/bwa-mem3) (the aligner used by this pipeline, derived from [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2))

  > Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS. 2019. doi: 10.1109/IPDPS.2019.00041

- [bedtools](https://pubmed.ncbi.nlm.nih.gov/20167835/)

  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Feb 8. PubMed PMID: 20167835; PubMed Central PMCID: PMC2832824.

- [Ensembl VEP](https://pubmed.ncbi.nlm.nih.gov/27288722/)

  > McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. doi: 10.1186/s13059-016-0974-4. PubMed PMID: 27288722; PubMed Central PMCID: PMC4893825.

- [fastp](https://pubmed.ncbi.nlm.nih.gov/30423086/)

  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890. doi: 10.1093/bioinformatics/bty560. PubMed PMID: 30423086; PubMed Central PMCID: PMC6129281.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- [FGUMI](https://github.com/fulcrumgenomics/fgumi)

  > Fulcrum Genomics. fgumi: Fulcrum Genomics UMI Tools (Rust implementation of fgbio). Computer software. Vers. 0.7.0. https://github.com/fulcrumgenomics/fgumi. doi: 10.5281/zenodo.18702466

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/)

  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

- [Picard](https://github.com/broadinstitute/picard)

  > Picard Toolkit. Broad Institute, GitHub Repository. https://github.com/broadinstitute/picard

- [Qualimap](https://doi.org/10.1093/bioinformatics/btv566)

  > Okonechnikov K, Conesa A, García-Alcalde F. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics. 2016 Jan 15;32(2):292-4. doi: 10.1093/bioinformatics/btv566. Epub 2015 Oct 1.

- [SeqKit](https://pubmed.ncbi.nlm.nih.gov/27706313/)

  > Shen W, Le S, Li Y, Hu F, Li Q. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Processing. PLoS One. 2016;11(10):e0163962. doi: 10.1371/journal.pone.0163962. PubMed PMID: 27706313.

- [SigProfilerMatrixGenerator](https://doi.org/10.1186/s12864-019-6041-2)

  > Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics. 2019;20(1):685. doi: 10.1186/s12864-019-6041-2.

- [VarDict](https://doi.org/10.1093/nar/gkw227)

  > Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, Dry JR. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016;44(11):e108. doi: 10.1093/nar/gkw227.

## Reference pipelines

deepUMIcaller was not built from scratch: parts of the workflow were taken directly from, or inspired by, the following open-source pipelines. We thank their authors for making their work publicly available.

- [nf-core/fastquorum](https://github.com/nf-core/fastquorum)

  > The 2022 version of this pipeline was taken as a reference backbone to build deepUMIcaller, starting from [nf-core/fastquorum](https://github.com/nf-core/fastquorum) at commit 09a6ae27ce917f2a4b15d2c5396acb562f9047aa. It was originally written by [Nils Homer](https://github.com/nh13) ([Fulcrum Genomics](https://www.fulcrumgenomics.com/)) and implements the [fgbio Best Practices FASTQ to Consensus Pipeline](https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md). If you use nf-core/fastquorum, please cite its Zenodo DOI: 10.5281/zenodo.11267672 (and doi: 10.5281/zenodo.10456900 for [fgbio](https://github.com/fulcrumgenomics/fgbio)).

- [cancerit/NanoSeq](https://github.com/cancerit/NanoSeq)

  > We used the ideas of filtering reads by unique mappability and counting the non-duplex reads to get another measurement of the VAF. The NanoSeq protocol and analysis are described in: Abascal F, Harvey LMR, Mitchell E, et al. Somatic mutation landscapes at single-molecule resolution. Nature. 2021;595(7865):154-159. doi: 10.1038/s41586-021-03477-4.

- [Kennedy-Lab-UW/Duplex-Seq-Pipeline](https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline)

  > We used the code for QCing the mutation calls with the mutations per position plots. The pipeline implements Duplex Sequencing, originally described in: Schmitt MW, Kennedy SR, Fox EJ, Loeb LA. Detection of ultra-rare mutations by next-generation sequencing. Proc Natl Acad Sci U S A. 2012;109(36):14508-13. doi: 10.1073/pnas.1208715109. The step-by-step protocol is described in: Kennedy SR, Schmitt MW, Fox EJ, Kohrn BF, Salk JJ, Ahn EH, Prindle MJ, Kuong KJ, Shen JC, Risques RA, Loeb LA. Detecting ultralow-frequency mutations by Duplex Sequencing. Nat Protoc. 2014 Nov;9(11):2586-606. doi: 10.1038/nprot.2014.170.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
