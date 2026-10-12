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

06 - Extension of the λ path beyond the StARS-selected boundary for the selected Pseudomonas-Z. tritici model (Supplementary Table 14)

Scripts are numbered in the order in which they are run.

## Software
R 4.2.2, SpiecEasi 1.1.3, pulsar 0.3.11, igraph 2.0.3, phyloseq 1.42.0, Matrix, dplyr.
SpiecEasi 1.1.3 does not check whether its singular value decomposition succeeded (`softSVT2` in `ADMM.cpp`). On some LAPACK builds this causes the error "Mat::operator(): index out of bounds" for high ranks; replacing the call by a fallback to the standard SVD routine (`svd_econ(..., "std")` if the default fails) solves this without changing results when the default routine converges.

# Summary
*Background*: Background: Plant-associated microbiota comprise diverse microbial species that shape community structure and host health. Yet biocontrol agents selected in the laboratory frequently fail in the field, in part because the interaction networks they enter are known only at coarse taxonomic resolution and from single timepoints. Here we resolve a cross-kingdom co-occurrence network in the wheat phyllosphere at strain level across a growing season, combining pangenome-informed, taxon-specific long-read amplicons for *Pseudomonas* and *Zymoseptoria tritici* with full-length 16S rRNA gene and ITS markers across 478 leaves from five timepoints, eight cultivars and three canopy heights in a naturally infected field.

*Results*: At genus level, co-occurrence was dominated by negative associations (81.2% of 165 edges). At strain resolution the pattern inverted: positive associations accounted for 88.4% of the *Pseudomonas*–fungal network and 96.8% of 917 Pseudomonas–Z. tritici associations. Neither interaction partners nor interaction strength were phylogenetically conserved, consistent with strain- rather than clade-level specificity. Negative co-occurrence modeling identified an assemblage of eleven candidate biocontrol taxa, comprising five *Pseudomonas* ASVs and six fungal genera, negatively associated with eight phytopathogenic fungal genera. A further nine Pseudomonas ASVs were positively associated with seven of these suppressors and may act as stabilizers that enhance their persistence. Testing these predictions with isolates from the same field (238 of ~3,000 species-typed), co-inoculation on leaf discs increased *Pseudomonas* biomass in 15 of 24 isolate–strain combinations, consistent with the predicted positive dominance. As predicted for *P. poae*, isolate Poa1 inhibited *Z. tritici* strain 1A5 in vitro more strongly than all 17 other field isolates (adjusted p ≤ 0.001) and was, with the second *P. poae* isolate, the only field isolate associated with reduced disease severity in planta.

*Conclusions*: Strain-resolved network inference recovers structure invisible at genus level, with potential antagonism concentrated in a small, strain-specific subset of the community and not predictable from phylogeny. By identifying both candidate antagonists and the taxa that may sustain them, this approach prioritizes consortia assembled from field co-occurrence rather than from antagonist screens alone, yielding a ranked set of candidates for targeted experimental testing.


*Keywords*: Network inference, co-occurrence, amplicon sequencing, microbiota, microbiome, wheat, phyllosphere, *Pseudomonas*, *Zymoseptoria tritici*
