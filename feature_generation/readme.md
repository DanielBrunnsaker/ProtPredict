# Feature generation/frequent pattern mining

Code and files needed to perform relational learning and frequent pattern mining over the content of the Datalog-database.

`proteomics/*` - files needed to recreate the feature-set used for Figures 2,3 and 6

`proteomics_noAA/*` - files needed to recreate the feature-set used for Figures 4, 5 and 6.


## File explanatons

`/aleph_orig.pl` - Source file for Aleph version 5 (needed to initialize aleph in SWI-Prolog).


`/proteomics.b` - Background file containing all of the information pertaining the pattern search. Explanations for the terms can be seen as comments inside the files themselves.


`/proteomics.f` - File containing all of the positive examples (deletant strains) used in the search. These are the examples that the frequent pattern miner tries to entail.


`/main.pl` - File that is called in `proteomics.b`. This file contains paths to all of the used relations (i.e. links to `knowledgeBase/Datalog/*.dl`). These are called to establish all existing relations in the search-space.


`saveOutput.pl` - Helper-function used to assist in saving the generated features to .txt files. This function is loaded in `proteomics.b`.


