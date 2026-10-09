# Citation
If you use this resource, please cite:

Stalder, L., Spescha, A., Maurhofer, M., Croll, D. (2026). High-resolution microbial network analysis identifies biocontrol assemblages in the wheat phyllosphere. *Microbiome*.

# Content
This repository contains the code for the critical analysis steps of the study: ASV inference, phyloseq object creation, network inference with rank selection and robustness assessment, network statistics and the pathogen suppressor classification. Code that only produces descriptive figures is not included.

00 - Helper functions for the network analyses, sourced by the 03 scripts

01 - Scripts and files for the dada2 pipeline to create ASV tables

02 - Scripts and files for phyloseq object creation

03 - SpiecEasi network calculation with the latent graphical model (sparse and low-rank, SLR) for the bacterial-fungal, Pseudomonas-fungal and Pseudomonas-Z. tritici networks (one script each): rank selection by EBIC, robustness analysis (relative Hamming distance, Kurtz et al. 2019) and export of the selected edges with partial correlations. The model fits take several hours per rank. The Pseudomonas-fungal script also runs the sensitivity analysis without the P. congelans ASVs S44 and S45 (set `exclude_asvs`).

04 - Network statistics (edge counts, positive/negative and within/cross-kingdom shares, hub degrees) and modularity tests on the edge tables written by 03; with `--ranks`, the network-level sensitivity to the SLR rank on the saved models of all ranks

05 - Pathogen suppressor, facilitator and stabilizer classification (Figure 5A/B) on the selected Pseudomonas-fungal model, and on the model without S44/S45. Uses `05_Meta_info_..._trophy.txt`, the phytopathogen annotation of all network taxa (Supplementary Table 8)

06 - Extension of the λ path beyond the StARS-selected boundary for the selected Pseudomonas-Z. tritici model (Supplementary Table 19)

Scripts are numbered in the order in which they are run.

## Software
R 4.3.3, SpiecEasi 1.1.3, pulsar 0.3.11, igraph 2.0.3, phyloseq 1.42.0, Matrix, dplyr.
SpiecEasi 1.1.3 does not check whether its singular value decomposition succeeded (`softSVT2` in `ADMM.cpp`). On some LAPACK builds this causes the error "Mat::operator(): index out of bounds" for high ranks; replacing the call by a fallback to the standard SVD routine (`svd_econ(..., "std")` if the default fails) solves this without changing results when the default routine converges.

# Summary
*Background*: Plant-associated microbiota comprise diverse microbial species that coexist and interact, influencing host community structure and plant health. However, our understanding of these interactions in field conditions and at strain resolution remains limited. Laboratory findings often fail to translate due to insufficient interaction network knowledge. Here, we apply pangenome-informed, taxon-specific long-read amplicons to resolve a cross-kingdom co-occurrence network within the wheat phyllosphere microbiota sampled at five timepoints across a growing season. 

*Results*: We performed in-depth monitoring of strains from the hub genus *Pseudomonas*, revealing strain-specificic cross-kingdom associations of *Pseudomonas*. Through negative co-occurrence modelling, we identified an assembly of eleven biocontrol taxa with the potential to suppress eight fungal pathogens under field conditions. Additional taxa with positive co-occurrence associations with this assembly were identified as putative persistence-enhancing stabilizer strains. We tested predicted interactions of *Pseudomonas* with the major fungal pathogen *Zymoseptoria tritici* in co-inoculation experiments with *Pseudomonas* isolates from the same field. Consistent with a species-level prediction, one of four *P. poae* isolates inhibited the pathogen *in vitro* more strongly than all other field isolates on one *Z. tritici* strain and showed the strongest trend towards reduced disease in planta. 

*Conclusions*: Our study demonstrates that high-resolution network inference can map microbial co-occurrence networks and identify candidate biocontrol assemblies at strain resolution across a growing season under field conditions. Our approach supports the design of sustainable biocontrol solutions.


*Keywords*: Network inference, co-occurrence, amplicon sequencing, microbiota, microbiome, wheat, phyllosphere, *Pseudomonas*, *Zymoseptoria tritici*
