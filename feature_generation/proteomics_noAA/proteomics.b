%search constraints
:- set(evalfn, coverage).
:- set(i,20).
:- set(search,ar).
:- set(pos_fraction, 0.00025).
:- set(max_features, 2049). %2048 features, 2049 because it starts from one
:- set(clauselength,10).

%mode declarations
:- modeh(*, networked(+gene)).
:- modeb(*, gene_metabolite(+gene, #metabolite)).
:- modeb(*, gene_pathway(+gene, #pathway)).
:- modeb(*, gene_chromosome(+gene, #chromosome)).
:- modeb(*, gene_disease(+gene, #disease)).
:- modeb(1, enzymatic_interaction(+gene,-gene)).
:- modeb(1, regulates(+gene, -gene, #type)). %arity 3
:- modeb(1, regulatedBy(+gene, -gene, #type)). %arity 3
:- modeb(*, gene_nullphenotype(+gene, #phenotype)).
:- modeb(*, gene_nullphenotypechemical(+gene, #phenotype,#chemical)). %arity 3
:- modeb(*, gene_hasproteindomain(+gene, #domain)).
:- modeb(*, involved_in(+gene,#function)).
:- modeb(*, enables(+gene, #function)).
:- modeb(*, part_of(+gene, #function)).
:- modeb(*, contributes_to(+gene, #function)).
:- modeb(*, colocalizes_with(+gene, #function)).
:- modeb(*, acts_upstream_of_or_within(+gene, #function)).
:- modeb(*, not_involved_with(+gene, #function)).
:- modeb(*, acts_upstream_of(+gene, #function)).
:- modeb(*, acts_upstream_of_negative_effect(+gene, #function)).
:- modeb(*, acts_upstream_of_or_within_positive_effect(+gene, #function)).
:- modeb(*, acts_upstream_of_positive_effect(+gene, #function)).
:- modeb(*, is_active_in(+gene, #function)).
:- modeb(*, located_in(+gene, #function)).

%determinations
:- determination(networked/1, gene_metabolite/2).
:- determination(networked/1, enzymatic_interaction/2).
:- determination(networked/1, regulates/3).
:- determination(networked/1, regulatedBy/3).
:- determination(networked/1, gene_pathway/2).
:- determination(networked/1, gene_nullphenotype/2).
:- determination(networked/1, gene_nullphenotypechemical/3).
:- determination(networked/1, gene_disease/2).
:- determination(networked/1, gene_chromosome/2).
:- determination(networked/1, gene_hasproteindomain/2).
:- determination(networked/1, located_in/2).
:- determination(networked/1, involved_in/2).
:- determination(networked/1, enables/2).
:- determination(networked/1, part_of/2).
:- determination(networked/1, contributes_to/2).
:- determination(networked/1, colocalizes_with/2).
:- determination(networked/1, acts_upstream_of_or_within/2).
:- determination(networked/1, not_involved_with/2).
:- determination(networked/1, acts_upstream_of/2).
:- determination(networked/1, acts_upstream_of_negative_effect/2).
:- determination(networked/1, acts_upstream_of_or_within_positive_effect/2).
:- determination(networked/1, acts_upstream_of_positive_effect/2).
:- determination(networked/1, is_active_in/2).

% show examples as boolean vectors, included to provide an output of the "induce features" mode, as per the aleph manual

:- set(portray_examples,true).
aleph_portray(train_pos):-
        setting(train_pos,File),
        show_features(File,positive).
aleph_portray(train_neg):-
        setting(train_neg,File),
        show_features(File,negative).
show_features(File,Class):-
        open(File,read,Stream),
        repeat,
        read(Stream,Example),
        (Example = end_of_file -> close(Stream);
                write_features(Example,Class),
                fail).
write_features(Example,_):-
        feature(_,(Example:- Body)),
        (Body -> write(1), write(' '); write(0), write(' ')),
        fail.
write_features(_,Class):-
	writeq(Class), nl.

:- ['main.pl']. %load kb
:- ['../saveOutput.pl']. % function to allow for saving the output