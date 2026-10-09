# Citation
If you use this resource, please cite:

Stalder, L., Spescha, A., Maurhofer, M., Croll, D. (2026). High-resolution microbial network analysis identifies biocontrol assemblages in the wheat phyllosphere.

# Content
This repository contains the code for the critical analysis steps of the study: ASV inference, phyloseq object creation, network inference with rank selection and robustness assessment, network statistics and the pathogen suppressor classification. Code that only produces descriptive figures is not included.

01 - Scripts and files for the dada2 pipeline to create ASV tables

02 - Scripts and files for phyloseq object creation

00 - Helper functions for the network analyses (sourced by the 03 scripts)

03 - SpiecEasi network calculation with the latent graphical model (sparse and low-rank, SLR) for the bacterial-fungal, Pseudomonas-fungal and Pseudomonas-Z. tritici networks: rank selection by EBIC, robustness analysis (relative Hamming distance, Kurtz et al. 2019), export of the selected edges with partial correlations. The Pseudomonas-fungal script also runs the sensitivity analysis without the P. congelans ASVs S44 and S45 (set `exclude_asvs`).

04 - Network statistics (edge counts, positive/negative and within/cross-kingdom shares, hub degrees), modularity tests, and network-level sensitivity to the SLR rank

05 - Pathogen suppressor, facilitator and stabilizer classification and Figure 5A/B. Uses `05_Meta_info_..._trophy.txt`, the phytopathogen annotation of all network taxa (Supplementary Table 8)

09 - Extension of the λ path beyond the StARS-selected boundary for the Pseudomonas-Z. tritici network (Supplementary Table 19)

## Order of execution
1. `03_SPIEC_EASI_network_*.Rmd` (one per network; the model fits take several hours per rank)
2. `04_network_statistics.R` on the edge tables written by 03 (and with `--ranks` on the saved models of all ranks)
3. `05_pathogen_suppressor_assessment.R` on the selected Pseudomonas-fungal model (and on the model without S44/S45)
4. `09` on the selected Pseudomonas-Z. tritici model (Supplementary Table 19)

## Software
R 4.3.3, SpiecEasi 1.1.3, pulsar 0.3.11, igraph 2.0.3, phyloseq 1.42.0, Matrix, dplyr.
SpiecEasi 1.1.3 does not check whether its singular value decomposition succeeded (`softSVT2` in `ADMM.cpp`). On some LAPACK builds this causes the error "Mat::operator(): index out of bounds" for high ranks; replacing the call by a fallback to the standard SVD routine (`svd_econ(..., "std")` if the default fails) solves this without changing results when the default routine converges.

# Summary
- Plant-associated microbiota comprise diverse microbial species that coexist and interact, influencing host community structure and plant health. However, our understanding of these interactions in field conditions and at strain resolution remains limited. Laboratory findings often fail to translate due to insufficient interaction network knowledge. Here, we apply pangenome-informed, taxon-specific long-read amplicons to resolve a cross-kingdom co-occurrence network within the wheat phyllosphere microbiota.
- We performed in-depth monitoring of strains from the hub genus Pseudomonas, revealing a high degree of strain-specificity of Pseudomonas interactions both within and across kingdoms.
- Through negative co-occurrence modelling, we identified a candidate assembly of eleven biocontrol taxa negatively associated with eight fungal pathogen genera. Additional taxa with positive co-occurrence associations with these candidates were identified as putative persistence-enhancing stabilizer strains. We validated the strain-specific interactions of Pseudomonas with the major fungal pathogen Zymoseptoria tritici using co-inoculation experiments with genotypes retrieved from the same field. Consistent with our prediction, a P. poae isolate was the most antagonistic towards the pathogen in vitro and showed the strongest trend towards reduced disease in planta.
- Our study demonstrates that high-resolution network inference can map microbial co-occurrence networks and predict strain-specific interaction patterns of biocontrol candidates that are persistently present under field conditions. Our approach supports the design of sustainable biocontrol assemblies.

Keywords: Network inference, co-occurrence, amplicon sequencing, microbiota, microbiome, wheat, phyllosphere, Pseudomonas, Zymoseptoria tritici
