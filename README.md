What is opy_Targets?
===============
Opentargets.org uses GraphQL API to explore it's content via coding.
This ensemble of functions aim is to make it easy to use the API via our
beloved coding language - python <3 



Python dependencies
-----------------

(requires Python 3)

* pandas
* requests
* json


Usage
=====

Example: <br>

 Lets say you want to get SNPS of gene you think is related to a disease


    import opy_targets as opy
    
    ENDOMETRIOSIS_EFO_ID = 'EFO_0001065'
    GENE_ESNG_ID = ENSG00000196208
    
    opy.get_snp_data(GENE_ESNG_ID,ENDOMETRIOSIS_ID).head(4)
    
output: <br>


|    | variantRsId   | variantId      | studyId      |   studySampleSize | publicationFirstAuthor   | label              |   chr |   location | disease_ID   | gene_related    |
|---:|:--------------|:---------------|:-------------|------------------:|:-------------------------|:-------------------|------:|-----------:|:-------------|:----------------|
|  0 | rs11674184    | 2_11581409_T_G | GCST004549   |            208903 | Sapkota Y                | intron_variant     |     2 |   11581409 | EFO_0001065  | ENSG00000196208 |
|  1 | rs11674184    | 2_11581409_T_G | GCST004549_2 |            208903 | Sapkota Y                | intron_variant     |     2 |   11581409 | EFO_0001065  | ENSG00000196208 |
|  2 | rs13394619    | 2_11587381_G_A | GCST001720   |             13997 | Nyholt DR                | intron_variant     |     2 |   11587381 | EFO_0001065  | ENSG00000196208 |
|  3 | rs77294520    | 2_11520829_G_C | GCST004549_2 |            208903 | Sapkota Y                | intergenic_variant |     2 |   11520829 | EFO_0001065  | ENSG00000196208 |


Most of the chances you don't know the genes related to the disease. NO PROBLEM! <br>

Lets say you want to know which SNP associate with a disease you study


    import opy_targets as opy
    
    #
    
    ENDOMETRIOSIS_EFO_ID = 'EFO_0001065'
    
    opy.get_SNP_df(ENDOMETRIOSIS_ID).head()
    
    
output:

|    | variantRsId   | variantId            | studyId                      |   studySampleSize | publicationFirstAuthor   | label              |   chr |   location | disease_ID   | gene_related    |
|---:|:--------------|:---------------------|:-----------------------------|------------------:|:-------------------------|:-------------------|------:|-----------:|:-------------|:----------------|
|  0 | rs11674184    | 2_11581409_T_G       | GCST004549                   |            208903 | Sapkota Y                | intron_variant     |     2 |   11581409 | EFO_0001065  | ENSG00000196208 |
|  1 | rs13394619    | 2_11587381_G_A       | GCST001720                   |             13997 | Nyholt DR                | intron_variant     |     2 |   11587381 | EFO_0001065  | ENSG00000196208 |
|  2 | rs77294520    | 2_11520829_G_C       | GCST004549_2                 |            208903 | Sapkota Y                | intergenic_variant |     2 |   11520829 | EFO_0001065  | ENSG00000196208 |
|  3 | rs58502716    | 2_11592421_G_GAATCAC | FINNGEN_R5_N14_ENDOMETRIOSIS |             77257 | FINNGEN_R5               | intron_variant     |     2 |   11592421 | EFO_0001065  | ENSG00000196208 |



You can change parameters such as the amount of genes retrived for each disease (25 default) and the sorting method of them (genetic association). check the functions docstring. 

Enjoy ;)

