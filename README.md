# personalized-GE-behaviour
Personalized cellular gene expression behaviour

## Flow of work
- METHODOLOGY
  - Catalogue expression datasets
    - microarray
      - affy
        - platform
      - illumina
        - platform
    - RNAseq
      - platform
  - Use the processed gene expression data
  - Identify responsive and responsivity genes
  - Quality check if any variability exist between replicates (if the study have biological replicates)
  - Perform delta co-expression analysis or WGCNA analysis or other network approaches
  - Create coposite network using both repsonsivity and responsive genes.
- ANALYSIS
  - Identify key plays in the network that can be associated with any
  phenotypes (viral titer or others)
  - Find TFs and co-TFs in the responsive and responsivity genes in the network.
  - Use RNAseq data to look for any long non-coding miRNA that can regulate downstream genes.


## LATER
- Download raw dataset
- Perform appropriate normalization method
- Remove batch effect and retain gender, age and other informations
- Step (3) to 
