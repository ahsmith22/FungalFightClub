# Fungal Fight Club

Data and analysis to accompany Smith et al., "Fungal Fight Club: Phylogeny and growth rate predict competitive outcomes among ectomycorrhizal fungi."

Code is split into 4 markdown files that can be run independently: 

- FFC_Model.Rmd: Code for calculating the growth rate based on modeling

- FFC_Plots.Rmd: Code for generating all plots and linear models

- FFC_Stats.Rmd: Code for all statistical analyses

- FFC_eoc.Rmd: Code for calculating the effect of competition (eoc) metric

### Data

- FFC Treatments - Controls SVS.csv: Raw data for self versus self control growth rate

- FFC Treatments - Controls.csv: Raw data for single control growth rate

- FFC Treatments - Pairwise.csv: Raw data for competition plate growth rate

- FFC Treatments - Reaction Types.csv: Raw data for reaction type scoring (used for index of antagonism metric)

- ffc_n_pair.csv: number of viable plates for each competitive treatment

- final_data_eoc.csv: Raw growth rate data, modeled growth rate calculation, and calculated eoc values

- finaldata1.csv: modeled growth rate calculation

### Plots

- RAxML_Tree.nex: Phylogenetic tress derived by RAxML analysis

- eocgenplot_r.png: Plot of eoc values across species

- eocplot_r.png/svg: Figure 3. Plot of eoc values for each competitive treatment

- infile.nex.con.tre: Phylogenetic tree file derived using ape package

- kiteplot.png/svg: Figure 2. Plot summarizing growth rates of single and svs controls across species

- mrBayesPost.nex: Phylogenetic tress derived by mrBayes analysis

- phyloplot.png/svg: Phylogenetic tree plot derived using ape package

- plot_EOCr_distbay.png/svg: Figure 5a. Plot of model results comparing eoc and phylogenetic distance

- plot_EOCr_rdiff.png/svg: Figure 5b. Plot of model results comparing eoc and growth rate distance

- plot_dist_r.png: Plot of model results comparing growth rate distance and phylogenetic distance

- sinplot_r.png: Summary of single control growth rates by species

- svsplot_r.png: Summary of svs control growth rates by species

### Tables

- Tuk_eocphdiff_table.htm: Supplemental Table 6. Significant TukeyHSD results comparing effects of pH on the EoC metric on growth rate

- Tuk_final.csv: Results of TukeyHSD testing on all competitions treatments

- Tuk_geoc_fin.csv: Results of TukeyHSD testing across species

- Tuk_geoc_table.htm: Results of TukeyHSD testing across species in table form

- Tuk_kite_fin.csv: Supplemental Table 5. Significant TukeyHSD results comparing single and svs control growth rates

- Tuk_sin_fin.csv: Results of TukeyHSD testing on single controls

- Tuk_sincont_table.htm: Supplemental Table 3. Significant TukeyHSD results comparing effects of pH on single control growth rate

- Tuk_svscont_table.htm: Supplemental Table 4. Significant TukeyHSD results comparing effects of pH on SvS control growth rate

- Tuk_svsdiff_table.htm: Supplemental Table 7. TukeyHSD results comparing effects of interspecific competition on growth rate compared to growth rate in intraspecific competition

- lmeocr.htm: Supplemental Table 8. Linear model results showing the significance of fungal identity, competitor identity, and pH in performance in competition 

- mod_EOCrdiff_distbay.htm: Table 1. ANOVA table comparing the two linear mixed effects models

- mod_rdist_phylodiff.htm: Supplemental Table 9. Linear model results showing the correlation between growth rate distance and phylogenetic distance

- sig_eoc_fin.csv: Results of TukeyHSD testing on all competitive treatments between pH conditions

- sig_svs_fin.csv: Results of TukeyHSD testing on svs controls
