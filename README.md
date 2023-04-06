# ProtPredict

Proteomic profiles reflect the functional readout of the physiological state of an organism, and
thus, increased understanding of what controls and defines the accumulated abundances of proteins is of
high scientific interest. Saccharomyces cerevisiae is a well studied model organism, as such there is a large
amount of structured knowledge on yeast systems biology contained in databases such as the Saccharomyces
Genome Database, and in highly curated genome-scale metabolic models like Yeast8. These are the product
of decades worth of experiments on multiple different modalities, these are abundant in information, and
adhere to semantically meaningful ontologies.

![alt text](https://github.com/DanielBrunnsaker/ProtPredict/blob/main/Schematic.png?raw=true)

By representing this prior knowledge in a richly expressive Datalog database (/knowledgeBase) we generated data
descriptors using relational learning that, when combined with standard approaches, allows us to accurately
predict protein abundances in S. cerevisiae in an explainable manner, connecting them to functional
annotations of the genotype and phenotypical observations, such as α-amino acid concentrations and
deviations in chronological lifespan. We further showcase this methodology on the proteins His4 and Ilv2,
successfully connecting qualitative biological concepts to quantified abundances. We argue that explainable
models improve the understanding of the implications of knowledge-based systems biology representations,
and enables them to be rationally improved

This repository contains the code need to reproduce the results of the study.

Additional data will need to be downloaded at Zenodo (link incoming).

# Frequent Pattern Mining

Frequent patterns are mined from Datalog database created from the Saccharomyces Genome Database (SGD), BioGRID, 
and Yeast8[1,2,3] and then used as features in a protein abundance predictions.  These patterns are mined using Aleph in Prolog using the WARMR algorithm, and by using sample meta-data (deletant strains) from a dataset by Messner et al. as positive examples [4,5].

In order to generate the features, you SWi-prolog needs to be installed, and you need to follow the following commands:
In the folder of the relevant dataset (feature_generation/proteomics or feature_generation/proteomics_noAA):

```
$ swipl
[aleph_orig].
?- read_all(proteomics).
?- induce_features.
```

In order to save the features as a .txt document, do the following:

**To save the positive examples in order:**
```
?- saveQueries('name.txt'). 
?- show(pos).
?- stopQueriesSaving().
```

**Saving the binary features:**
```
?- saveQueries('name_features.txt'). 
?- show(train_pos).
?- stopQueriesSaving().
```

**Saving the logic programs connected to each feature:**
```
?- saveQueries('name_expl.txt'). 
?- show(features).
?- stopQueriesSaving().
```

# Model training & analysis

Models are trained using XGBoost[cite] on both standard propositional data (protein abundances and metabolite concentrations) 
and relational features (binary features).

See the following notebooks:

**scripts/notebooks/ProteinsFromAA.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using metabolite (amino acid) concentrations.

**scripts/notebooks/ProteinsFromILP.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using propositionalized logic programs. Additonally, visualize feature importances using SHAP and Gain.

**scripts/notebooks/ProteinsFromAA.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using a combination of propositionalized logic programs and metabolite concentrations. Additionally, visualize comparative performance between featuresets, global feature importances and the predictive capacity of amino acids, and beeswarm-plots for SHAP-values for specific proteins.


# References

[1] J. Michael Cherry, et al. Saccharomyces Genome Database: the genomics resource of budding yeast. Nucleic Acids Research, 40(Database issue):D700–705, January 2012

[2] Rose Oughtred, et al. The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. 
Protein Science: A Publication of the Protein Society, 30(1):187–200, January 2021

[3] Hongzhong Lu, et al. A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism. Nature Communications, 10(1):3586, August 2019. Number: 1

[4] Ashwin Srinivasan. The Aleph Manual.

[5] Ross D. King, et al. Warmr: a data mining tool for chemical data. Journal of Computer-Aided Molecular Design, 15(2):173– 181, February 2001.



