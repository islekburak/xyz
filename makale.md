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
  * [The Subfamily Specific Residues](#the-subfamily-specific-residues)
  * [Extracellular Ligands](#extracellular-ligands)
  * [Structural Comparison of FZD4, FZD9 and FZD10](#structural-comparison-of-fzd4--fzd9-and-fzd10)
    + [TM6 and ECL3](#tm6-and-ecl3)
- [DISCUSSION](#discussion)
- [CONCLUSION](#conclusion)
- [BIBLIOGRAPHY](#bibliography)
- [SUPPLEMENTARY DATA](#supplementary-data)
- [FUNDING SOURCES](#funding-sources)

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

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/pathway.png)|
|:--:|
| **Figure 1: Schematic Representation of Canonical WNT Signaling** |
The core components of WNT/β-catenin pathway are Frizzled receptors and low-density lipoprotein receptor-related proteins 5 or 6 (LRP5/6) coreceptors. Upon the agonist (WNT) stimulation from the outer membrane, Dishevelled (DVL) protein is activated via signaling through Frizzled and LRP5/6. Dsh inhibits and heads to disassembly of the β-catenin destruction complex that is composed of Axin, adenomatosis polysis coli (APC), casein kinase 1α (CK1α), and glycogen synthase kinase 3β (GSK3-β). Thus β-catenin can not be phosphorylated and targeted to proteasome for degradation. Consequently, unphosphorylated β-catenins translocate into nucleus to bind T-cell factor/Lymphoid enhancer factor (TCF/LEF) transcription factors and induce the signaling cascade via transcription trigger **(Figure 1, left)**.
In the absence of WNT, the β-catenin destruction complex that is located in the cytosol can phosphorylate β-catenins by both GSK3-β and CK1α. Phosphorylated β-catenins are targeted by proteasomes for proteolysis, therefore they are not be able to migrate into nucleus to bind transcription factors and start the transcription **(Figure 1, right)**.<sup>**4,5,6**</sup>
### Extracellular Frizzled Binding Ligands and the Role of FZD4 in Familial Exudative Vitreoretinopathy
There are some evidences that soluble Frizzled related proteins (sFRPs-which act as a negative regulator in the interaction between WNTs and FZDs), R-spondin and Norrin can also bind to Frizzleds whereas it is known fact that the primary agonists (extracellular ligands) of Frizzled receptors are the WNT proteins <sup>**7**</sup>. Norrin ligands specifically are able to bind FZD4 receptors through the whole Frizzled family members. On the N-terminus of Frizzled receptors, there is a region called cysteine rich domain (CRD) that consists of approximately 120-125 amino acids. Activation of downstream signaling upon the stimulation of Frizzled CRD domains by WNTs, maintains the regulation of cell polarity via PCP pathway and modulates the proliferation, cell fate and cell development during embryogenesis <sup>**8**</sup>. Previous studies demonstrated that Norrins structurally mimic WNT ligands for binding to FZD4/CRD <sup>**9**</sup>. 
Stimulation of FZD4/CRD with Norrin ligands, leads to recruitment DVL to FZD4 and thus, disassembly of β-catenin destruction complex can inhibit β-catenins likewise in the canonical WNT signaling. Since the endothelial cells are located in the form of a single cell layer in the surface of blood vessels, the WNT/β-catenin signaling in this tissue has a crucial role in the development of blood vessels and especially in retinal vascularization. Thus, possible defects in the members of WNT/β-catenin pathway that triggered with the binding of Norrin to FZD4,  may result in multiple ocular malformations such as retinal hypovascularization, ophthalmic diseases and a symptom of blindness <sup>**10**</sup>. Eventually, the disruption in Norrin induced canonical WNT signaling, leads to accumulation of β-catenins that can not be collected by the nucleus, in turn the retinal vascularization stops.
Moreover, a member of the erythroblast transformation-specific (ETS) family that is named as ETS-related gene (ERG) helps FZD4 to trigger canonical WNT signalling and to be transcribed by controlling stability of β-catenins in the cytosol. Thus, transcribed ERG and β-catenins lead to interact with vascular endothelial cadherin (VEC) to promote vessel stability, vascular integrity and vessel growth <sup>**11**</sup>.
Familial exudative vitreoretinopathy (FEVR) is an hereditary anomaly due to failure in retinal development and characterized by abnormal peripheral retinal vascularization. Several studies demonstrated that some mutations in NDP gene, FZD4, LRP5 and another coreceptor tetraspanin (TSPAN12) cause FEVR disease via blocking to activation of WNT/β-catenin pathway <sup>**12**</sup>.

### Significance of this study
In this study, we intented to understand the evolutionary history, functions, interactions, and significance of the Frizzled receptor family using different genomic tools in silico. Since FZD4 biological function is the most finely studied among Frizzleds, firstly we have focused on the FZD4/9/10 subclade of these 7TM proteins. Stimulated by recent insights into the activation mechanisms of Class F receptors from structural and functional analysis of Frizzleds, we aim to summarize what we know about molecular details of ligand binding, agonist-driven conformational changes and FZD4 receptor activation. A better understanding of receptor activation mechanisms will allow us to engage in structure- and mechanism-driven drug discovery with the potential to develop more isoform-selective and potentially pathway-selective drugs for human therapy.
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
The Maximum likelihood phylogenetic analysis and sequence homology clustering reveal two supergroup of Frizzled proteins in the presence of SMO outgroup. FZD3/6 cluster and another cluster of FZDs that contains 3 different subclades: FZD1/2/7, FZD4/9/10 and FZD5/8. (**Figure 2**).

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/tree.png)|
|:--:|
| **Figure 2: Phylogenetic analysis and homology clustering of frizzled receptors gene family**|
| *This analysis involved 1202 amino acid sequences of 10 frizzled receptor genes in whole eukaryotic lineage. The sequences inferred by Maximum likelihood method were selected by CD-HIT algorithm (with the clustering threshold of 90% identity) as representatives. The outgroup SMO was indicated by red branches and nodes, while FZD4/9/10 clade is shown by blue. Colored labels are implied human sequences (**a**). Sequence homology clustering of only human Frizzled receptors gene family reveals the same clusters. FZD4/9/10 clade is clearly distinguishable in blue frame. The legend indicates the similarity range: "0" means completely similar whereas "1" means completely different.  (**b**).*|

The evolutionary tree obtained with phylogenetic analysis advocates that the eukaryotic lineage-Frizzled receptors diverged by the following duplication events: The first major duplication event is separating the ancestral genes of FZD3/6 supergroup from ancestor genes of other FZDs (**Figure 2a**). The second duplication event which has occurred in the other supergroup that comprises the remaining members of the FZDs is dividing the ancestral genes of clade FZD1, FZD2, and FZD7 from the ancestral genes of other subclades. The ensuing distribution is conducted by the third duplication event which results in the separation of the ancestor genes of FZD5, FZD8 subclade from the ancestor genes of FZD4/9/10 subclade. Throughout the evolutionary process, FZD4 is located as an outlier in the FZD4/9/10 subclade while FZD9 and FZD10 have partaken of the same ancestor genes with the fourth duplication event.  Sequence homology clustering through human Frizzled receptors has also reinforced these distributions among supergroups and subclades (**Figure 2b**).

## The Subfamily Specific Residues
SPIAL (Specifity in alignments) analysis has revealed the possible subfamily specific residues among the members of Human FZD4/9/10 clade. We identified a region with 10 specific residues on the Human FZD4 extracellular cysteine-rich domain (CRD; aa 40-161) that are clustered prominently different from the other members of this clade (**Figure 3**). These residues are mostly clustered near the disulphide bridges (**Supplemental Figure S1**). On the other hand, we observed that N144 is specific for FZD4 which occurs as N-acetylation site in the squence. It is known that disulfide bonds and N-acetylation sites stabilize the folded structures of proteins. Thus, possible mutations on these sites may affect correct protein folding and moreover may be responsible from the dimerization type of FZD4.

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/figure2.png)|
|:--:|
| **Figure 3: Subfamily Specific Residues of the members of FZD4/9/10 clade** |
|*Human FZD4, FZD9, and FZD10 receptors and the specific residues through their domains within the secondary structures were shown with lollipop graphs. While lollipops are representing the specific residues for each receptor, CRD domains and the secondary structures were indicated by different colors (**a**). On the right-hand side, the receptors were shown in 3D structures with their CRDs. The same colors were used in 3D visuals to imply the sub-structures that were shown on the left-hand side, whereas the linker region -the sequence located between the CRD and the TM1- is shown with purple. The gray-colored beads were used to signify the specific residues(**b**).*|

2D visuals obtained from amino acid sequences has been unveiled the differences between the secondary structures of these receptors. FZD9 and FZD10 have longer TM6 and ECL3 compared to FZD4 (**Figure 3a**). In the clade of FZD4/9/10, we noticed that the specific residues were substantially conglomerated on the coupling site of the CRD region and transmembrane domains (**Figure 3b**). It is clearly seen that there is an abundance of specific residues especially in the intracellular site of FZD9. Although the other members have some specific residues in the region of TM6-ECL3-TM7, notably in FZD9, our analysis has not been revealed any specific residues in FZD4.
The marked residues were located on the exterior surface of transmembrane helices in the FZD4 receptor. FZD9 and FZD10 have outnumbered specific residues on the transmembrane domain, and they mostly observed within inter-helical site (**Figure 3b**).

## Extracellular Ligands

As a functionally distinct member of Frizzled receptor family, the Smoothened receptor does not interact with WNTs by t its CRD region<sup>**A**</sup>. In the previous studies, the other 10 members of Frizzled family have been shown as the essential receptors of WNT proteins <sup>**7,B**</sup>. Wnt proteins were divided into two different classes depending on the pathways in which they are involved. β-catenin dependent pathway (canonical) includes Wnt-1, Wnt-2, Wnt-3, Wnt-8a, Wnt-8b, Wnt-10a and Wnt-10b while Wnt-4, Wnt-5a, Wnt-5b, Wnt-6, Wnt-7a, Wnt-7b and Wnt-11 are the members of β-catenin independent pathway (non-canonical)<sup>**C**</sup>. 

In order to investigate possible individuality and specificity of WNT-FZD combinations, we scanned the reviews and our findings (**Figure 4**) demonstrated that the FZD1 receptor can interact with Wnt-1, Wnt-2, Wnt-3, Wnt-3a, Wnt-5a, Wnt-5b, Wnt-7a and Wnt-7b <sup>**D1,D2**</sup>. The presented data in IUPHAR Database suggests that FZD2 has been proven to bind to Wnt-2, Wnt-3, Wnt-3a, Wnt-5a, Wnt-6, Wnt-7a and Wnt-8b while FZD3 can bind to Wnt-2, Wnt-3a and Wnt-5a <sup>**D3**</sup>. There is a shred of evidence for the direct interaction of FZD4 with Wnt-5a <sup>**D4,D5**</sup>, Wnt-3a <sup>**D6**</sup>, Wnt-2 <sup>**D7**</sup>, Wnt-7b <sup>**D2**</sup>, Wnt-9b <sup>**D8**</sup> and the FZD4 specific ligand Norrin <sup>**D9**</sup>. The immunoprecipitation studies have shown that Wnt-3a, Wnt-5a, Wnt-7a, Wnt-7b and Wnt-9b are able to bind to the FZD5 receptor. Moreover, FZD6 can interact with Wnt-3a, Wnt-4, Wnt-5a, Wnt-5b, and Wnt-7a whereas FZD7 is able to bind to Wnt-3, Wnt-3a, Wnt-5a and Wnt-7a <sup>**D3,D8**</sup>. Besides, we realized that the previous studies have shown possible interactions between the FZD8 receptor and Wnt-2 <sup>**D10**</sup>, Wnt-3a <sup>**D11,D12**</sup>, Wnt-7b <sup>**D8**</sup>, Wnt-9b <sup>**D13**</sup> and Xenopus Wnt-8 (XWNT-8) <sup>**D14**</sup>. It has been previously reported that FZD9 can bind to Wnt-7a and activates the JNK pathway <sup>**D15**</sup> while Wnt-2 induces the TCF transcription led by FZD9 <sup>**D16**</sup>. On the other hand, FZD10 is known to interact with Wnt-7b <sup>**D17**</sup> and Wnt-7a <sup>**D18**</sup>, Wnt-9b <sup>**D8**</sup>.
According to reports published Wnt2 can be an important ligand for FZD5 receptor <sup>**D7**</sup> whereas some information indicates FZD4 and FZD5 can interact with Wnt-10b, especially in MCF-7 adhesion cells <sup>**D19**</sup>. 


| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/FZD-WNT_interactions.png)|
|:--:|
| **Figure 4: The overview of Frizzled receptors-WNT interactions**|
|*The possible combinations between WNTs and Frizzled receptors were examined by several studies. This schematic representation provides WNTs and FZDs that were distributed as a reference to their guide trees which were positioned on the top and the left. The outlier member SMO has been dislocated from the guide tree on the left side. The WNTs (orange) which have been reported into the interaction with FZDs were labeled as **reviewed** whereas the others (blue) were labeled as **unreviewed**.*|

In summary, the data suggests that FZD4 is able to involved in both canonical and non-canonical WNT signaling compared to other members of FZD4/9/10 clade. We could not detect any specific or individual combination with Wnt proteins except Norrin-FZD4 specific binding, however we found a clue for only the combination of WNT-2b and FZD4 <sup>**D20**</sup> interestingly. Also the data retrieved from Human Protein Atlas shows that Wnt-2b is highly expressed in the eye <sup>**D21**</sup>. 

FZD4 is known for the ability of binding both WNTs and Norrin. Yet, the ligand selectivity of FZD4 remains poorly understood. Due to there is still a lack of structural information for WNTs, we structurally aligned FZD8/CRD in complex with XWNT-8 (PDB ID:4F0A) onto the region of FZD4/CRD to explore the residues that would be able to responsible for the FZD4 ligand selection (**Figure 5**). XWNT-8 exhibits two-point attachment on FZD8/CRD through its thumb and index finger . We demonstrated these regions on FZD4/CRD monomer obtained from FZD4/CRD homodimer structure (PDB ID: 5UWG) (**Figure 5a**). On the other hand, upon stimulation by WNT ligands, it has been shown that FZD4/CRD can dimerize. Surprisingly we detected a specific residue (L132) which is located on the FZD4/CRD dimeric interface (**Figure 5b**). Aside from these findings, we also focused on the possible binding site of human FZD4 and FZD4 specific ligand-Norrin. The human FZD4/CRDs combined with Norrin homodimer was obtained from the complex structure (PDB ID: 5CL1) and the specific residues were signified on 3D structure (**Figure 5c**).

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/overview.png)|
|:--:|
|**Figure 5: The Specific Residues on the possible interaction site of WNT, Norrin and FZD4/CRD**|
| *CRD was shown with brushed yellow, Norrin was shown with blue, XWNT-8 was shown with ice blue. The marked residues were shown in boxes indicate the specific resiudes on FZD4/CRD (green beads for Norrin binding region, and red beads for WNT binding region respectively) which are located on the possibble interaction site of FZD/CRD and WNT ligands (**a**). FZD4/CRD dimer structure stimulated by WNT ligand (PDB ID:5UWG) was shown with the specific residues that we found on FZD4/CRD (**b**) Norrin homodimer-FZD4/CRD combination (PDB ID: 5CL1) was inferred with available marked residues on the hypothetical FZD4/CRD interaction site (**c**). Gray beads represent the other specific residues we found on FZD4/CRD.*|

Since FZD4 is the only protein that can bind Norrin ligands in Frizzled receptor family and we found a cluster of marked residues on CRD region, we also examined the interaction of human Norrin and human FZD4 to see whether it is related to dimerization type or not (**Figure 6**).



Even though there is still a lack of knowledge about the trigger mechanism of dimerization, it is previously shown that FZD4 receptors are able to construct a dimer structure within CRD domains to activate cell signaling and biological functions <sup>**9**</sup>. Hence, we examined the specific residues that we found on FZD4/CRD region on dimeric 3D structure of human FZD4 CRD in complex with human Norrin ligands (PDB ID:5BQC). We recognized that 4 of 10 marked residues were positioned on the interaction site of human Norrin and human FZD4/CRD. In the consideration of LYS109 and MET157 were reported <sup>**18**</sup> as binding sites of Norrin and FZD4/CRD before, we supposed that HIS156 and GLY161 also play a selective role in Norrin binding for human FZD4. Besides, some mutations (M105V and M157V) on FZD4/CRD have been associated with Familial Exudative Vitreoretinopathy (FEVR) <sup>**19**</sup> in the previously performed disease-related studies. These pieces of evidence suggest that these residues could be individual for FZD4 to selecting to bind Norrin ligands.

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/FZD4_norrin_dimer.png)|
|:--:|
|**Figure 6: 3D visual of Norrin in Complex with Human FZD4/CRD**|
| *FZD4/CRD was shown with brushed yellow whereas Norrin was shown with blue. The schematic representation of Norrin bound FZD4/CRD (**a**) was depicted in 3D visual on the right-hand side. The subfamily-specific residues located on Norrin-FZD4/CRD binding region were represented as beads. LYS109, HIS156, MET157, and GLY161 were colored in green and the other marked residues on FZD4/CRD were represented as gray beads. Colorless part of 3D Norrin-FZD4/CRD represents the other biological unit of dimer structure(**b**) .*|

Eventually we distributed the specific residues that we found on FZD4/CRD into the following two groups:
In the first group which is located on the site of the WNT thumb, we signified the specific residues V101, L132, and N144. In the second group, K109, H156, and M157 were classified in the Norrin interaction site which corresponds to the WNT index finger site. Additionally, G161 has been observed on the interaction site for both WNT and Norrin ligands. As a result that the second group can interact with both of these extracellular ligands, the residues of the first group seem possibly important for selecting WNT ligands for FZD4.

It would appear that these two different groups seem to be characteristic for ligand bias to activate FZD4. They also may be deterministic for CRD positioning in dimerization mechanism. Norrin homodimers locate two distant FZD4/CRDs in contrary directions (**Figure 6**), whereas FZD4/CRDs construct homodimer structure in WNT-bound form.


## Structural Comparison of FZD4, FZD9 and FZD10

### TM6 and ECL3

In the xSMO cyclopamine-bound CRD structure, the E and F rings of cyclopamine, which acts as an agonist when bound to the CRD, clash with helix VI of the cholesterol-bound multi-domain SMO structure, indicating that helix VI and ECL3 should shift outward even further when cyclopamine is bound to the CRD in a multi-domain SMO (Fig. 3e). (https://www.nature.com/articles/ncomms15383)



| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/structural_comparison.png)|
|:--:|
| **Figure 7:** Structural Comparison of Human FZD4, FZD9 and FZD10|
|*This figure clearly shows that FZD4 has a truncated TM6 domain compared to FZD9 and FZD10. Disulfide bridges were shown with yellow, the linker region colored in magenta, TM6 domains represented with blue, and extracellular loop 3 demonstrated in orange. 3 extra disulfide bonds which FZD4 has were labeled by residue numbers.*|

The members of subclade FZD4/9/10  were modeled with CRDs and compared through their structural differences. Data retrieved from Uniprot shows that FZD9 and FZD10 have 5 disulfide bridges within their CRDs whereas FZD4 has 3 more bridges (Cys181-Cys200, Cys204-Cys282, Cys302-Cys377) (**Figure 7**). 


# DISCUSSION

# CONCLUSION

Subfamily specific rediues of FZD4 that we found are mostly clustered near the disulphide bridges (Figure 3). On the other hand, we found that N144 is specific for FZD4 which occurs as N-acetylation site. It is known that disulfide bonds and N-acetylation sites stabilize the folded structures of proteins. Thus, possible mutations on these sites may affect correct protein folding and moreover may be responsible from the dimerization type of FZD4. K109 and M157 were reported  <sup>**18**</sup> as binding sites of Norrin and FZD4 before (Figure 4), and some mutations (M105V and M157V) have been associated with Familial Exudative Vitreoretinopathy (FEVR) <sup>**19**</sup>.  Hence, possible variations in these residues may induce incorrect folding of FZD4/CRD, and consequently, defects in ligand binding that may have crucial role in FEVR diesase pathway. Further functional studies of FZD4 are required to assess whether a change in these residues plays any functional role in FEVR.

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

**A)** [Povelones M, Nusse R. The role of the cysteine-rich domain of Frizzled in Wingless-Armadillo signaling. EMBO J. 2005 Oct 5;24(19):3493-503. doi: 10.1038/sj.emboj.7600817. Epub 2005 Sep 15. PMID: 16163385; PMCID: PMC1276175.](https://pubmed.ncbi.nlm.nih.gov/16163385/)

**B)**[Bhanot P, Brink M, Samos CH, Hsieh JC, Wang Y, Macke JP, Andrew D, Nathans J, Nusse R. A new member of the frizzled family from Drosophila functions as a Wingless receptor. Nature. 1996 Jul 18;382(6588):225-30. doi: 10.1038/382225a0. PMID: 8717036.](https://pubmed.ncbi.nlm.nih.gov/8717036/)

**C)**[Ackers I, Malgor R. Interrelationship of canonical and non-canonical Wnt signalling pathways in chronic metabolic diseases. Diab Vasc Dis Res. 2018;15(1):3-13. doi:10.1177/1479164117738442](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5752873/)

**D1)**[1](https://pubmed.ncbi.nlm.nih.gov/10557084/)

**D2)** [1](https://pubmed.ncbi.nlm.nih.gov/15923619/)

**D3)**[1](https://pubmed.ncbi.nlm.nih.gov/24032637/)

**D4)**[1](https://pubmed.ncbi.nlm.nih.gov/12958364/)

**D5)**[1](https://pubmed.ncbi.nlm.nih.gov/16602827/)

**D6)**[1](https://pubmed.ncbi.nlm.nih.gov/18772438/)

**D7)**[1](https://pubmed.ncbi.nlm.nih.gov/18302287/)

**D8)**[1](https://pubmed.ncbi.nlm.nih.gov/28733458/)

**D9)**[1](https://pubmed.ncbi.nlm.nih.gov/15035989/)

**D10)**[1](https://pubmed.ncbi.nlm.nih.gov/23815780/)

**D11)**[1](https://pubmed.ncbi.nlm.nih.gov/17576136/)

**D12)**[1](https://pubmed.ncbi.nlm.nih.gov/24885675/)

**D13)**[1](https://pubmed.ncbi.nlm.nih.gov/18509025/)

**D14)**[1](https://pubmed.ncbi.nlm.nih.gov/22653731/)

**D15)**[1](https://pubmed.ncbi.nlm.nih.gov/15705594/)

**D16)**[1](https://pubmed.ncbi.nlm.nih.gov/12138115/)

**D17)**[1](https://pubmed.ncbi.nlm.nih.gov/11786918/)

**D18)**[1](https://pubmed.ncbi.nlm.nih.gov/18567805/)

**D19)**[1](https://pubmed.ncbi.nlm.nih.gov/27853307/)

**D20)**[1](https://pubmed.ncbi.nlm.nih.gov/31359032/)

**D21)**[1](https://www.proteinatlas.org/humanproteome/tissue)

# SUPPLEMENTARY DATA

| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/FZDs_CRD_human_MSA.png)|
|:--:|
| **Figure S1: Multiple sequence alignments of Human Frizzled CRDs** |
|*The clade of FZD4/9/10 is highlighted with green throughout the whole Frizzled receptor family members in humans. The specific residues are shown as bold letters, whereas the vertical residue numbers correspond to FZD4 specific residues. Cysteine residues are colored yellow while pairs of cysteines are labeled per a bridge by low case characters in the residue numeration line (disulfide bridges). N-glycosylation sites are marked with red.*|


| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/FZDs_linker_human_MSA.png)|
|:--:|
| **Figure S2: Multiple sequence alignments of Human Frizzled linker regions** |
|*The linker region positioned between CRD and TM1 for all Frizzled receptors was sequenced and aligned. N-glycosylation sites were indicated with red, while the specific residues marked as bold and vertical amino acid position. Disulfide bond C181-C200 was inferred with lowercase "k". To get a clear aspect, twenty five residues of the human FZD8 are disregarded and showed as "<25>".*|


| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/FZDs_human_MSA_TM6-ECL3-TM7.png)|
|:--:|
| **Figure S3: Multiple sequence alignments of Human Frizzled TM6-ECL3-TM7 regions**|
|*The specific residues on the FZD4/9/10 clade are shown as bold letters. Cysteine residues are colored yellow while N-glycosylation site is marked with red. Transmembrane domains and Extracellular loop were shown with different colors.*|


| ![](/home/islekburak/Desktop/toplantı/makale/src/yeni/png/fzd_domains.png)|
|:--:|
| **Table S1:The Residues of the Human FZD4/9/10 Clade Receptors**|
|*The data is retrieved from GpcrDB and Uniprot respectively*|

# FUNDING SOURCES
The present work was supported by the Research Fund of European Molecular Biology Organisation (EMBO). and grants from The Scientific and Technological Research Council of Turkey(TUBITAK) Project No. xxxxx
The funding sources were not involved in the study design; in the collection, analysis and interpretation of data; in the writing of the article; and in the decision to submit the article for publication.