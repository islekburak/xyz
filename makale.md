- [Evolutionary Analysis reveals Unique Features of Frizzled 4 receptor](#evolutionary-analysis-reveals-unique-features-of-frizzled-4-receptor)
- [ABSTRACT](#abstract)
- [BACKGROUND](#background)
  * [The Frizzled Receptors](#the-frizzled-receptors)
    + [Frizzleds in Cell Signaling](#frizzleds-in-cell-signaling)
    + [Extracellular Frizzled Binding Ligands and the Role of FZD4 in Familial Exudative Vitreoretinopathy](#extracellular-frizzled-binding-ligands-and-the-role-of-fzd4-in-familial-exudative-vitreoretinopathy)
    + [Significance of this study](#significance-of-this-study)
- [METHODS](#methods)
    + [Obtaining Proteomes, Performing BLAST and Clustering](#obtaining-proteomes--performing-blast-and-clustering)
    + [Alignments and Phylogenetic Analysis](#alignments-and-phylogenetic-analysis)
    + [Sequence Retrieval, Mining and Removal](#sequence-retrieval--mining-and-removal)
    + [Subfamily Specific Residues Search](#subfamily-specific-residues-search)
- [RESULTS](#results)
  * [Phylogenetic analysis](#phylogenetic-analysis)
  * [Subfamily Specific Residues of the members of FZD4, FZD9, and FZD10 clade](#subfamily-specific-residues-of-the-members-of-fzd4--fzd9--and-fzd10-clade)
  * [Extracellular Ligands](#extracellular-ligands)
    + [Norrin](#norrin)
    + [WNTs](#wnts)
  * [Structural Comparison of FZD4, FZD9 and FZD10](#structural-comparison-of-fzd4--fzd9-and-fzd10)
    + [TM6 and ECL3](#tm6-and-ecl3)
- [DISCUSSION](#discussion)
- [CONCLUSION](#conclusion)
- [BIBLIOGRAPHY](#bibliography)
- [SUPPLEMENTARY DATA](#supplementary-data)
- [FUNDING SOURCES](#funding-sources)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


Evolutionary Analysis reveals Unique Features of Frizzled 4 receptor
=====
# ABSTRACT
GPCRs are reported as highly related with few physiological pathway
components -such as hormones and neurotransmitters-, regulation of
social graces, regulation of the inflammatory responses and responsible
for vision, olfaction and taste mechanisms.<sup>**1**</sup><sup>**,**</sup><sup>**2**</sup> GPCRs' helices have
back and forth loops that are connected by flexible linkers at the EC
and IC ends, therefore they have defined as 7-transmembrane (7TM)
proteins for their 7 helices structures and localization through the
cell membrane. Now there are various structures of 7TM proteins that
have been solved, approximately 850 of 7TM proteins are defined and they
classified into five subgroups by GRAFS classification roughly.<sup>**3**</sup>
They classified with the features while they have different states,
active to inactive forms in a complex with a G-proteins. In this study,
our purpose is to understand the evolutionary history, functions,
interactions, and significance of GPCRs using different genomic tools
*in silico*.
Accordingly, we have focused on the F class of these GPCRs. Hence we performed BLAST using human sequences for Frizzled/SMO family (11 were obtained from Uniprot-KB database) via high-performance computing. The results of BLAST are evaluated and the protein sequences are clustered by sequence similarities to create a phylogenetic tree. After getting orhtolog information for human Frizzled receptors, we intended to create a
domain pattern for differential diagnosis of Frizzled GPCRs. Thus it
helps to determine functional equivalence to the sets of orthologs of
GPCRs. Then we will perform phylogenetic profiling of the categorized
subfamilies. With the predictions using computer-based genomic tools and
algorithms, we will be able to build an MSA database classified by
subfamilies in the future.

**Key Words:** GPCRs; Frizzled; Norrin, FEVR

The present work was supported by the Research Fund of European
Molecular Biology Organisation (EMBO). Project No. <mark>xxxxx</mark>

# BACKGROUND

GPCRs are reported as highly related with few physiological pathway components -such as hormones and neurotransmitters-, regulation of social graces, regulation of the inflammatory responses and responsible for vision, olfaction and taste mechanisms.<sup>**1,2**</sup> GPCRs' helices have back and forth loops that are connected by flexible linkers at the EC and IC ends, therefore they have defined as 7-transmembrane (7TM) proteins for their 7 helices structures and localization through the cell membrane. Now there are various structures of 7TM proteins that have been solved, approximately 850 of 7TM proteins are defined and they classified into five subgroups by GRAFS classification (glutamate, rhodopsin, secretin, adhesion and Frizzled/SMO receptor family) roughly.<sup>**3**</sup> They classified with the features while they have different states, active to inactive forms in a complex with a G-proteins.

## The Frizzled Receptors
The Frizzled/SMO family of GPCRs is characterized by the long cysteine-rich domain in N-terminus. They are also known with the Fz, Frizzled and SMO_Human domains (UniProt accession: PF01392, PF01534 and ,Q99835) This receptor family is consisting of 10 members: FZD1, FZD2, FZD3, FZD4, FZD5, FZD6, FZD7, FZD8, FZD9, FZD10. The Smoothened receptor is distantly related with this family; they have similar structures as having N-terminus CRD and highly conserved seven transmembrane domains but they are functionally differentiated. While FZDs are the key components of WNT signaling pathway, Smoothened receptors can not bind WNTs as extracellular ligands and are included in Hedgehog signaling pathway.

### Frizzleds in Cell Signaling
Regarding to the existence of β-catenin, Frizzleds can be considered to include in the following two pathways: Canonical and non-canonical WNT signaling. In canonical WNT pathway, signaling can be occured via stimulation by agonist (WNT) or not. Likewise, non-canonical WNT signaling consists of two different pathways: Planar Cell Polarity (PCP), and FZD/Ca<sup>2+</sup> signalings.

| ![](/png/pathway.png)|
|:--:|
| **Figure 1: Schematic Representation of Canonical WNT Signaling** |
The core components of WNT/β-catenin pathway are Frizzled receptors and low-density lipoprotein receptor-related proteins 5 or 6 (LRP5/6) coreceptors. Upon the agonist (WNT) stimulation from the outer membrane, Dishevelled (DVL) protein is activated via signaling through Frizzled and LRP5/6. Dsh inhibits and heads to disassembly of the β-catenin destruction complex that is composed of Axin, adenomatosis polysis coli (APC), casein kinase 1α (CK1α), and glycogen synthase kinase 3β (GSK3-β). Thus β-catenin can not be phosphorylated and targeted to proteasome for degradation. Consequently, unphosphorylated β-catenins translocate into nucleus to bind T-cell factor/Lymphoid enhancer factor (TCF/LEF) transcription factors and induce the signaling cascade via transcription trigger **(Figure 1, left)**.
In the absence of WNT, the β-catenin destruction complex that is located in the cytosol can phosphorylate β-catenins by both GSK3-β and CK1α. Phosphorylated β-catenins are targeted by proteasomes for proteolysis, therefore they are not be able to migrate into nucleus to bind transcription factors and start the transcription **(Figure 1, right)**.<sup>**4,5,6**</sup>
### Extracellular Frizzled Binding Ligands and the Role of FZD4 in Familial Exudative Vitreoretinopathy
There are some evidences that soluble Frizzled related proteins (sFRPs-which act as a negative regulator in the interaction between WNTs and FZDs), R-spondin and Norrin can also bind to Frizzleds whereas it is known fact that the primary agonists (extracellular ligands) of Frizzled receptors are the WNT proteins <sup>**7**</sup>. Norrin ligands specifically are able to bind FZD4 receptors through the whole Frizzled family members. On the N-terminus of Frizzled receptors, there is a region called cysteine rich domain (CRD) that consists of approximately 120-125 amino acids. Activation of downstream signaling upon the stimulation of Frizzled CRD domains by WNTs, maintains the regulation of cell polarity via PCP pathway and modulates the proliferation, cell fate and cell development during embryogenesis <sup>**8**</sup>. Previous studies demonstrated that Norrins structurally mimic WNT ligands for binding to FZD4-CRD <sup>**9**</sup>. 
Stimulation of FZD4-CRD with Norrin ligands, leads to recruitment DVL to FZD4 and thus, disassembly of β-catenin destruction complex can inhibit β-catenins likewise in the canonical WNT signaling. Since the endothelial cells are located in the form of a single cell layer in the surface of blood vessels, the WNT/β-catenin signaling in this tissue has a crucial role in the development of blood vessels and especially in retinal vascularization. Thus, possible defects in the members of WNT/β-catenin pathway that triggered with the binding of Norrin to FZD4,  may result in multiple ocular malformations such as retinal hypovascularization, ophthalmic diseases and a symptom of blindness <sup>**10**</sup>. Eventually, the disruption in Norrin induced canonical WNT signaling, leads to accumulation of β-catenins that can not be collected by the nucleus, in turn the retinal vascularization stops.
Moreover, a member of the erythroblast transformation-specific (ETS) family that is named as ETS-related gene (ERG) helps FZD4 to trigger canonical WNT signalling and to be transcribed by controlling stability of β-catenins in the cytosol. Thus, transcribed ERG and β-catenins lead to interact with vascular endothelial cadherin (VEC) to promote vessel stability, vascular integrity and vessel growth <sup>**11**</sup>.
Familial exudative vitreoretinopathy (FEVR) is an hereditary anomaly due to failure in retinal development and characterized by abnormal peripheral retinal vascularization. Several studies demonstrated that some mutations in NDP gene, FZD4, LRP5 and another coreceptor tetraspanin (TSPAN12) cause FEVR disease via blocking to activation of WNT/β-catenin pathway <sup>**12**</sup>.

### Significance of this study
In this study, we intented to understand the evolutionary history, functions, interactions, and significance of the Frizzled receptor family using different genomic tools in silico. Since FZD4 is the most finely studied, firstly we have focused on the FZD4, FZD9 and FZD10 subclade of these 7TM proteins. Stimulated by recent insights into the activation mechanisms of Class F receptors from structural and functional analysis of Frizzleds, we aim to summarize what we know about molecular details of ligand binding, agonist-driven conformational changes and FZD4 receptor activation. A better understanding of receptor activation mechanisms will allow us to engage in structure- and mechanism-driven drug discovery with the potential to develop more isoform-selective and potentially pathway-selective drugs for human therapy.
Because of their unique divergence profile, GPCRs have not been successfully dissected into subfamily specific phylogenetic trees by the automated tools. Therefore, identifying the set of proteins conserving the ancestral function, and distinguishing them from diverged homologs are essential for accurate evaluation of the conservation status of amino acids.
In the light of the knowledge of their structural behavior, we will be able to consider new approaches to discover new computer-aided drugs to cure disorders or pathological defects in the molecular homeostasis resulting from the changes in the activation mechanisms of Frizzled receptors.



# METHODS
### Obtaining Proteomes, Performing BLAST and Clustering

10 different human frizzled protein sequences (Q9UP38, Q14332, Q9NPG1, Q9ULV1, Q13467, O60353, O75084, Q9H461, O00144, Q9ULW2), 1 smoothened receptor sequence (Q99835) were downloaded with the complete eukaryotic proteomes (n=1202) from Uniprot-KB database<sup>**13**</sup> using FTP and stored in HPC by FASTA format. 11 proteins that were retrieved from Uniprot-KB were used as a query against the eukaryotic lineage subject and BLAST was performed. To see the spectrum of proteins through eukaryotic lineages clearly, the number of maximum target sequences was extended from 50 to 5000. For each Frizzled protein, all target sequences are selected up to a range of third human proteins in BLAST results.
Selected human targets for each Frizzled protein, were evaluated in subgroups by their complete sequence similarities rate using Cluster Omega<sup>**14**</sup>.

### Alignments and Phylogenetic Analysis

The sequences of selected subgroups were aligned using MAFFT version 7 <sup>**15**</sup> with the E-INS-i parameter to assign conserved motifs and carrying multiple domains. Trimal was used for removal of poorly aligned regions and then a phylogenetic tree was generated for each subgroup of Frizzled proteins by using Raxml with maximum likelihood estimation.

### Sequence Retrieval, Mining and Removal

Since only the sets of orthologous genes are expected to reflect the underlying species evolution, we trimmed the phylogenetic tree to get them without paralogous clades by using ete3 toolkit in Python.

### Subfamily Specific Residues Search

A web-server called Spial (specificity in alignments) <sup>**16**</sup> was used (at 90% consensus and 95% specifity treshold) to pair multiple sequence alignments of clades reciprocally that were obtained before. Thus we could highlight positions which are conserved in both alignments and positions which are specific to either alignment. After getting specific residues, the gaps were removed from alignments to a reference sequence and mapped the highlighted positions to residues on 3D structures of human Frizzled proteins. 3D structures of Human Frizzled proteins within CRD domains were generated using Phyre2 web portal <sup>**17**</sup>.


# RESULTS

## Phylogenetic analysis
The Maximum likelihood phylogenetic analysis and sequence homology clustering reveal two supergroup of Frizzled proteins in the presence of SMO outgroup. FZD3, FZD6 cluster and another cluster of FZDs that contains 3 different subclades: FZD1, FZD2, FZD7; FZD4, FZD9, FZD10 and FZD5, FZD8. (**Figure 2**).

| ![](/png/tree.png)|
|:--:|
| **Figure 2: Phylogenetic analysis and homology clustering of frizzled receptors gene family**|
| *This analysis involved 1202 amino acid sequences of 10 frizzled receptor genes in whole eukaryotic lineage. The sequences inferred by Maximum likelihood method were selected by CD-HIT algorithm (with the clustering threshold of 90% identity) as representatives. The outgroup SMO was indicated by red branches and nodes, while FZD4,FZD9,FZD10 clade is shown by blue. Colored labels are implied human sequences (**a**). Sequence homology clustering of only human Frizzled receptors gene family reveals the same clusters. FZD4,FZD9,FZD10 clade is clearly distinguishable in blue frame. The legend indicates the similarity range: "0" means completely similar whereas "1" means completely different.  (**b**).*|

The tree obtained with phylogenetic analysis is advocated that the eukaryotic lineage Frizzled receptors are varied by the following duplication events: The first major duplication event is separating the ancestral genes of FZD3 and FZD6 supergroup from ancestor genes of other FZDs (**Figure 2a**). The second duplication event which has occurred in the other supergroup that comprises the remaining members of the FZDs is dividing the ancestral genes of clade FZD1, FZD2, and FZD7 from the ancestral genes of other subclades. The ensuing distribution is conducted by the third duplication event which results in the separation of the ancestor genes of FZD5, FZD8 subclade from the ancestor genes of FZD4, FZD9, and FZD10 subclade. Throughout the evolutionary process, FZD4 is located as an outlier in the FZD4, FZD9, and FZD10 subclade while FZD9 and FZD10 have partaken of the same ancestor genes with the fourth duplication event.  Sequence homology clustering through human Frizzled receptors has also reinforced these distributions among supergroups and subclades (**Figure 2b**).

## The Subfamily Specific Residues
SPIAL (Specifity in alignments) analysis has revealed the possible subfamily specific residues among the members of Human FZD4, FZD9 and FZD10 clade. We identified a region with 10 specific residues on the Human FZD4 extracellular cysteine-rich domain (CRD; aa 40-161) that are clustered prominently different from the other members of this clade (**Figure 3**).

| ![](/png/figure2.png)|
|:--:|
| **Figure 3: Subfamily Specific Residues of the members of FZD4, FZD9, and FZD10 clade** |
|*Human FZD4, FZD9, and FZD10 receptors and the specific residues through their domains within the secondary structures were shown with lollipop graphs. While lollipops are representing the specific residues for each receptor, CRD domains and the secondary structures were indicated by different colors (**a**). On the right-hand side, the receptors were shown in 3D structures with their CRDs (**b**). The same colors were used in 3D visuals to imply the sub-structures that were shown on the left-hand side, whereas the linker region -the sequence located between the CRD and the TM1- is shown with purple. The gray-colored beads were used to signify the specific residues.*|

*Genel olarak CRD ve TM birleşme bölgelsinde bir kümelenme var hepsinde.

*FZD4^te TM6-ECL3 bölgesi için spesifik resid yokken diğerlerinde var. Özellikle 9'da daha fazla.

*FZD9'un intracellular bölgesinde bir kümlenme mvcut.

*FZD4'te TM üzerindekiler daha çok helix dış bölgesindeyken (hydrophilic head group?), FZD10'da sanki domainler arasında bağlantı kurabilecek bölgelerde (hydrophobic core group).

## Extracellular Ligands
### Norrin




| ![](/png/FZD4_norrin_dimer.png)|
|:--:|
| **Figure 4:** 3D visual of Norrin in Complex with Human FZD4-CRD |
| *FZD4-CRD is shown with brushed yellow whereas Norrin is shown with blue. The subfamily specific residues that were located on Norrin-FZD4 CRD binding region are represented as beads. LYS109 and MET157 were reported as Norrin-FZD4 binding sites. Disulphide bridges were shown in yellow licorice.*|

### WNTs

| ![](/png/FZD-WNT_interactions.png)|
|:--:|
| **Figure 5:** Frizzled receptors and WNT interactions|
|*ABCDEFGH*|


| ![]()|
|:--:|
| **Figure 6:** *WNT8A-PDB* |
|*ABCDEFGH*|

## Structural Comparison of FZD4, FZD9 and FZD10

### TM6 and ECL3


| ![](/png/structural_comparison.png)|
|:--:|
| **Figure 7:** Structural Comparison of Human FZD4, FZD9 and FZD10|
|*ABCDEFGH*|

# DISCUSSION

# CONCLUSION

Subfamily specific rediues of FZD4 that we found are mostly clustered near the disulphide bridges (Figure 3). On the other hand, we found that N144 is specific for FZD4 which occurs as N-acetylation site. It is known that disulfide bonds and N-acetylation sites stabilize the folded structures of proteins. Thus, possible mutations on these sites may affect correct protein folding and moreover may be responsible from the dimerization type of FZD4. K109 and M157 were reported  <sup>**18**</sup> as binding sites of Norrin and FZD4 before (Figure 4), and some mutations (M105V and M157V) have been associated with Familial Exudative Vitreoretinopathy (FEVR) <sup>**19**</sup>.  Hence, possible variations in these residues may induce incorrect folding of FZD4 CRD and, consequently, defects in ligand binding that may have crucial role in FEVR diesase pathway. Further functional studies of FZD4 are required to assess whether a change in these residues plays any functional role in FEVR.

# BIBLIOGRAPHY

**1)** [Rosenbaum, D.M., S.G. Rasmussen, and B.K. Kobilka, *The structure and function of G-protein-coupled receptors.* Nature, 2009.
**459**(7245): p. 356-63.](https://www.nature.com/articles/nature08144)

**2)** [Wu, F., et al., *Structure and Function of
Peptide-Binding G Protein-Coupled Receptors.* J Mol Biol, 2017.
**429**(17): p. 2726-2745.](https://linkinghub.elsevier.com/retrieve/pii/S0022-2836(17)30323-6)

**3)** [Wootten, D., et al., *Mechanisms of
signalling and biased agonism in G protein-coupled receptors.* Nat Rev
Mol Cell Biol, 2018. **19**(10): p. 638-653.](https://www.nature.com/articles/s41580-018-0049-3)

**4)** [Vallée A, Lecarpentier Y. Alzheimer Disease: Crosstalk between the Canonical Wnt/Beta-Catenin Pathway and PPARs Alpha and Gamma. Front Neurosci. 2016;10:459. Published 2016 Oct 19. doi:10.3389/fnins.2016.00459](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069291/)

**5)** [Dolatshad NF, Hellen N, Jabbour RJ, Harding SE, Földes G. G-protein Coupled Receptor Signaling in Pluripotent Stem Cell-derived Cardiovascular Cells: Implications for Disease Modeling. Front Cell Dev Biol. 2015;3:76. Published 2015 Dec 9. doi:10.3389/fcell.2015.00076](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673467/)

**6)** [Kohn AD, Moon RT. Wnt and calcium signaling: beta-catenin-independent pathways. Cell Calcium. 2005 Sep-Oct;38(3-4):439-46. doi: 10.1016/j.ceca.2005.06.022. PMID: 16099039.](https://pubmed.ncbi.nlm.nih.gov/16099039/)

**7)** [Schulte G, Bryja V. The Frizzled family of unconventional G-protein-coupled receptors. Trends Pharmacol Sci. 2007 Oct;28(10):518-25. doi: 10.1016/j.tips.2007.09.001. Epub 2007 Sep 19. PMID: 17884187.](https://pubmed.ncbi.nlm.nih.gov/17884187/)

**8)** [Milhem RM, Ali BR. Disorders of FZ-CRD; insights towards FZ-CRD folding and therapeutic landscape. Mol Med. 2019 Dec 31;26(1):4. doi: 10.1186/s10020-019-0129-7. PMID: 31892318; PMCID: PMC6938638.](https://pubmed.ncbi.nlm.nih.gov/31892318/)

**9)** [Chang TH, Hsieh FL, Zebisch M, Harlos K, Elegheert J, Jones EY. Structure and functional properties of Norrin mimic Wnt for signalling with Frizzled4, Lrp5/6, and proteoglycan. Elife. 2015;4:e06554. Published 2015 Jul 9. doi:10.7554/eLife.06554](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4497409/)

**10)** [Fujimura N. WNT/β-Catenin Signaling in Vertebrate Eye Development. Front Cell Dev Biol. 2016;4:138. Published 2016 Nov 30. doi:10.3389/fcell.2016.00138](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5127792/)

**11)** [Birdsey GM, Shah AV, Dufton N, Reynolds LE, Osuna Almagro L, Yang Y, Aspalter IM, Khan ST, Mason JC, Dejana E, Göttgens B, Hodivala-Dilke K, Gerhardt H, Adams RH, Randi AM. The endothelial transcription factor ERG promotes vascular stability and growth through Wnt/β-catenin signaling. Dev Cell. 2015 Jan 12;32(1):82-96. doi: 10.1016/j.devcel.2014.11.016. PMID: 25584796; PMCID: PMC4292982.](https://pubmed.ncbi.nlm.nih.gov/25584796/)

**12)** [Xiao H, Tong Y, Zhu Y, Peng M. Familial Exudative Vitreoretinopathy-Related Disease-Causing Genes and Norrin/β-Catenin Signal Pathway: Structure, Function, and Mutation Spectrums. J Ophthalmol. 2019 Nov 16;2019:5782536. doi: 10.1155/2019/5782536. PMID: 31827910; PMCID: PMC6885210.](https://pubmed.ncbi.nlm.nih.gov/31827910/)

**13)** [UniProt Consortium. UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Res. 2021 Jan 8;49(D1):D480-D489. doi: 10.1093/nar/gkaa1100. PMID: 33237286; PMCID: PMC7778908.](https://academic.oup.com/nar/article/49/D1/D480/6006196)

**14)** [Sievers F, Higgins DG. Clustal Omega for making accurate alignments of many protein sequences. Protein Sci. 2018 Jan;27(1):135-145. doi: 10.1002/pro.3290. Epub 2017 Oct 30. PMID: 28884485; PMCID: PMC5734385.](https://pubmed.ncbi.nlm.nih.gov/28884485/)

**15)** [Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010. Epub 2013 Jan 16. PMID: 23329690; PMCID: PMC3603318.](https://pubmed.ncbi.nlm.nih.gov/23329690/)

**16)** [Wuster A, Venkatakrishnan AJ, Schertler GF, Babu MM. Spial: analysis of subtype-specific features in multiple sequence alignments of proteins. Bioinformatics. 2010 Nov 15;26(22):2906-7. doi: 10.1093/bioinformatics/btq552. Epub 2010 Sep 29. PMID: 20880955; PMCID: PMC2971580.](https://pubmed.ncbi.nlm.nih.gov/20880955/)

**17)** [The Phyre2 web portal for protein modeling, prediction and analysis Kelley LA et al. Nature Protocols 10, 845-858 (2015)](https://www.nature.com/articles/nprot.2015.053)

**18)** [Shen G, Ke J, Wang Z, Cheng Z, Gu X, Wei Y, Melcher K, Xu HE, Xu W. Structural basis of the Norrin-Frizzled 4 interaction. Cell Res. 2015 Sep;25(9):1078-81. doi: 10.1038/cr.2015.92. Epub 2015 Jul 31. PMID: 26227961; PMCID: PMC4559814.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4559814/)

**19)** [Bang I, Kim HR, Beaven AH, Kim J, Ko SB, Lee GR, Kan W, Lee H, Im W, Seok C, Chung KY, Choi HJ. Biophysical and functional characterization of Norrin signaling through Frizzled4. Proc Natl Acad Sci U S A. 2018 Aug 28;115(35):8787-8792. doi: 10.1073/pnas.1805901115. Epub 2018 Aug 13. Erratum in: Proc Natl Acad Sci U S A. 2018 Nov 6;115(45):E10807. PMID: 30104375; PMCID: PMC6126767.](https://pubmed.ncbi.nlm.nih.gov/30104375/)

# SUPPLEMENTARY DATA

| ![](/png/FZDs_CRD_human_MSA.png)|
|:--:|
| **Figure S1: Multiple sequence alignments of Human Frizzled Receptors** |
|*The clade of FZD4, FZD9, and FZD10 is highlighted with green throughout the whole Frizzled receptor family members in humans. The specific residues are shown as bold letters, whereas the vertical residue numbers correspond to FZD4 specific residues. Cysteine residues are colored yellow while pairs of cysteines are labeled per a bridge by low case characters in the residue numeration line (disulfide bridges). N-glycosylation sites are marked with red.*|


| ![](/png/FZDs_linker_human_MSA.png)|
|:--:|
| **Figure S2:** linker region |
|*ABCDEFGH*|


| ![](/png/FZDs_human_MSA_TM6-ECL3-TM7.png)|
|:--:|
| **Figure S3:** TM6+EC3 sequence |
|*The specific residues on the FZD4,FZD9,and FZD10 clade are shown as bold letters. Cysteine residues are colored yellow while N-glycosylation site is marked with red. Transmembrane domains and Extracellular loop were shown with different colors.*|


||
|:--:|
| **Figure S4:** oligomerization |
|*ABCDEFGH*|


| ![](/png/fzd_domains.png)|
|:--:|
| **Table S1:The residue ranges of domains and regions for the members of FZD4, FZD9, and FZD10 clade**|
|*ABCDEFGH*|

# FUNDING SOURCES
The present work was supported by the Research Fund of European Molecular Biology Organisation (EMBO). Project No. xxxxx


The work was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation – 427840891; KO 5463/1-1) and grants from Karolinska Institutet, the Swedish Research Council (2017–04676; 2019–01190), the Swedish Cancer Society (CAN2017/561), the Novo Nordisk Foundation (NNF17OC0026940, NNF19OC0056122, NNF20OC0063168), Wenner-Gren Foundations (UPD2018-0064), Emil and Wera Cornells Stiftelse.

The funding sources were not involved in the study design; in the collection, analysis and interpretation of data; in the writing of the article; and in the decision to submit the article for publication.
