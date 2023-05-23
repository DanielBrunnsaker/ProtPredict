# ProtPredict

Proteomic profiles reflect the functional readout of the physiological state of an organism.
An increased understanding of what controls and defines protein abundances is of high scientific interest.
Saccharomyces cerevisiae is a well studied model organism, and there is a large amount of structured
knowledge on yeast systems biology in databases such as the Saccharomyces Genome Database, and highly
curated genome-scale metabolic models like Yeast8. These data-sets, the result of decades of experiments,
are abundant in information, and adhere to semantically meaningful ontologies.

![alt text](https://github.com/DanielBrunnsaker/ProtPredict/blob/main/graphical_abstract.png?raw=true)

By representing this knowledge in an expressive Datalog database we generated data descriptors
using relational learning that, when combined with supervised machine learning, enables us to predict protein
abundances in an explainable manner. We learnt predictive relationships between protein abundances,
function and phenotype; such as α-amino acid accumulations and deviations in chronological lifespan.
We further demonstrate the power of this new methodology on the proteins His4 and Ilv2, connecting
qualitative biological concepts to quantified abundances. Interpretable models enable us to better understand
knowledge-based representations in systems biology.

This repository contains the code need to reproduce the results of the study.

Additional data will need to be downloaded at Zenodo (link incoming).

# Install Python dependencies

Add requirements.txt etc.

# Install SWI-Prolog

Follow download and install instructions [here](https://www.swi-prolog.org/download/stable). It can also be installed using package managers such as apt, snap, and brew.

# Frequent Pattern Mining

Frequent patterns are mined from Datalog database created from the Saccharomyces Genome Database (SGD), BioGRID, 
and Yeast8[1,2,3] and then used as features in a protein abundance predictions.  These patterns are mined using Aleph in Prolog using the WARMR algorithm, and by using sample meta-data (deletant strains) from a dataset by Messner et al. as positive examples [4,5,6].

A pattern in this context can take the following form, represented as a logic program:

```
Gene(A) :=
  RegulatedBy(A, B, Transcription factor),
  nullPhenotype(B, Abnormal chronological lifespan),
  InvolvedIn(A, One − carbon metabolic process)
```

This pattern can be interpreted as: genes (A) which are involved
in the one-carbon metabolic process, and which are regulated by
a transcription factor (B) whose deletion causes the cell to have
an abnormal chronological lifespan.

## Feature generation in aleph

In order to generate the features, SWI-prolog (tested with v7.6.3) needs to be installed on your system, and you need to run the following commands in the folder of the relevant dataset (`feature_generation/proteomics` for the relational features only analysis or `feature_generation/proteomics_noAA` for use in combination with metabolite concentration values):

```
$ swipl
[aleph_orig].
?- read_all(proteomics).
?- induce_features.
```

In order to save the features as .txt documents, run the following commands in sequence:

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

<Write the shell command as well>

Alternatively, the data/feature-sets accompanied with explanatory logic programs can be downloaded at [LINK].

# Model training & analysis

Models are trained using XGBoost on both standard propositional data (protein abundances and metabolite concentrations) and relational features (binary features/logic programs).

See the following notebooks:

`scripts/notebooks/ProteinsFromAA.ipynb` - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using metabolite (amino acid) concentrations.

`scripts/notebooks/ProteinsFromILP.ipynb` - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using propositionalized logic programs. Additonally, visualize feature importances using SHAP and Gain.

`scripts/notebooks/ProteinsFromAAILP.ipynb` - Train and evaluate models used to predict protein abundances in *S. cerevisiae* using a combination of propositionalized logic programs and metabolite concentrations. Additionally, visualize comparative performance between featuresets, global feature importances and the predictive capacity of amino acids, and beeswarm-plots for SHAP-values for specific proteins.

# References

[1] J. Michael Cherry, et al. Saccharomyces Genome Database: the genomics resource of budding yeast. Nucleic Acids Research, 40(Database issue):D700–705, January 2012

[2] Rose Oughtred, et al. The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. 
Protein Science: A Publication of the Protein Society, 30(1):187–200, January 2021

[3] Hongzhong Lu, et al. A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism. Nature Communications, 10(1):3586, August 2019. Number: 1

[4] Ashwin Srinivasan. The Aleph Manual.

[5] Ross D. King, et al. Warmr: a data mining tool for chemical data. Journal of Computer-Aided Molecular Design, 15(2):173– 181, February 2001.

[6] Christoph B. Messner, et al. The Proteomic Landscape of Genome-Wide Genetic Perturbations, May 2022. Pages: 2022.05.17.492318
