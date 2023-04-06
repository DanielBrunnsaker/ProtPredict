# ProtPredict

Proteomic profiles reflect the functional readout of the physiological state of an organism, and
thus, increased understanding of what controls and defines the accumulated abundances of proteins is of
high scientific interest. Saccharomyces cerevisiae is a well studied model organism, as such there is a large
amount of structured knowledge on yeast systems biology contained in databases such as the Saccharomyces
Genome Database, and in highly curated genome-scale metabolic models like Yeast8. These are the product
of decades worth of experiments on multiple different modalities, these are abundant in information, and
adhere to semantically meaningful ontologies.

![alt text](https://github.com/DanielBrunnsaker/ProtPredict/blob/main/Schematic.png?raw=true)

By representing this prior knowledge in a richly expressive Datalog database we generated data
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
and the metabolic atlas[cite] and then used as features in a protein abundance predictions.  These patterns are mined using Aleph in Prolog[cite], and by
using sample meta-data (deletant strains) from a dataset by Messner et al. as positive examples.

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

?- saveQueries('name.txt'). 

?- show(pos).

?- stopQueriesSaving().


**Saving the binary features:**

?- saveQueries('namec_features.txt'). 

?- show(train_pos).

?- stopQueriesSaving().


**Saving the logic programs connected to each feature:**

?- saveQueries('name_expl.txt'). 

?- show(features).

?- stopQueriesSaving().


# Model training & analysis

Models are trained using XGBoost[cite] on both standard propositional data (protein abundances and metabolite concentrations) 
and relational features (binary features).

See the following notebooks:

**scripts/notebooks/ProteinsFromAA.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using metabolite (amino acid) concentrations.

**scripts/notebooks/ProteinsFromILP.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using propositionalized logic programs.

**scripts/notebooks/ProteinsFromAA.ipynb** - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using a combination of propositionalized logic programs and metabolite concentrations.






