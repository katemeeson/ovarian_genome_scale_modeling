# ovarian_genome_scale_modeling
## Repository to accompany transcriptomics integration method developed during the third year of my PhD

###   Constraint-based modelling performed on the 'Human-GEM-annotated.xml' 2020 version of the Human1 genome-scale model, developed by Robinson et al, 2020 (PMID: 32209698) [(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7331181/)]

###   Transcriptomics data was used to constrain generic GEM to six cell line-specific ovarian cancer models (CCLE, DepMap Public 22Q2, 'CCLE_expression.csv') [(https://depmap.org/portal/download/all/)]

###   Six cell line-specific models were constrained: three high-grade (CAOV3, COV318, OAW28), and three low-grade (59M, HEYA8, OV56)

###   Constraint-based modeling used FBA and own integration code to create reaction constraint dictionaries, which were integrated using COBRApy and MEWpy
