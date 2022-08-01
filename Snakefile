# pipeline for host-parasite STRING

# this pipeline performs homology correction of cooccurrence scores, 
# calculates direct functional textmining scores from co-occurrence and NLP scores, 
# performs orthology transfer of functional and physical textmining scores 
# and aggregates all results 

# TODO
# load modules
# config.file with directories: no input file directories in rules
# make sure files are zipped in the end
# write everything to snakemake results folder
# C2 profile https://github.com/Snakemake-Profiles/pbs-torque


import os
from subprocess import check_output

configfile: "config.yaml"
print("config: ", config)

workdir: config['workdir']
print("PWD: ", os.getcwd())

# I need the pipeline's path to execute the scripts
print("pipeline basedir: ", workflow.basedir)


# find number of parasite-host pairs
p_h_pairs_file = config["parasite_host_pairs_path"]
n_p_h_pairs = int(check_output(["wc", "-l", p_h_pairs_file]).split()[0])
# create the wildcards to iterate over splits
PARASITE_HOST_SPLITS = list(range(int(n_p_h_pairs / 50) + 1))
# pad with 0
PARASITE_HOST_SPLITS = [str(i).rjust(3, "0") for i in PARASITE_HOST_SPLITS]
print("PARASITE_HOST_SPLITS", PARASITE_HOST_SPLITS)

# this is for security: command snakemake in this directory will only run this rule
rule first: 
    shell:
        "echo 'first rule'"

###################################
####### Textmining channel ########
###################################

# sort cooccurrence scores by taxid for more efficient homology correction
rule do_sort_cooccurrence:
    input:
        # Cooccurrence scores of host-parasite PPIs (A-B and B-A)
        cooccur_scores = config["cooccur_scores_path"]
    output:
        # Cooccurrence scores of host-parasite PPIs, sorted by taxids
        cooccur_sorted = "all_host_parasite_database_cooccur_interact_global_only_parasites_sorted.tsv"
    run: 
        # Sort input file by taxids
        shell("{workflow.basedir}/scripts/sort_cooccurrence_scores_taxids.sh {input.cooccur_scores} {output.cooccur_sorted} 2> do_sort_cooccurrence.err")


# Calculate homology score for all PPIs with direct co-occurrence score and apply homology correction
rule do_homology_correction:
    input:
        # Cooccurrence scores of host-parasite PPIs (A-B and B-A), sorted by taxids
        cooccur_sorted = rules.do_sort_cooccurrence.output.cooccur_sorted,
        # Self normalized bitscores from SIMAP all vs. all searches, in both directions
        simap_bitscores = config["simap_data_path"] + "selfnorm-bitscores/"

    output: 
        # Co-occurrence scores with additional column homology score and corrected co-occurrence score
        cooccur_corrected = "all_host_parasite_database_cooccur_interact_global_only_parasites_corrected.tsv"

    run: 
        # Run homology_correction_cooccurrence.py
        shell("{workflow.basedir}/scripts/homology_correction_cooccurrence.py {input.simap_bitscores} {input.cooccur_sorted} {output.cooccur_corrected}")


# double nlp scores so that they are A-B and B-A, just like the cooccurrence scores
rule do_double_nlp:
    input:
        nlp_scores_single = config["nlp_scores_path"]

    output:
        nlp_scores = "host_parasite_physical_textmining_direct.tsv"

    shell:
        # repeat all lines with taxids and proteins switched
        """cat <(cat {input.nlp_scores_single}) <(awk -F"\\t" '{{OFS="\\t";}} {{print $1,$4,$5,$2,$3,$6}}' {input.nlp_scores_single}) > {output.nlp_scores}"""

# Raise corrected direct cooccurrence scores to direct NLP Complex Portal 
# calibrated scores to obtain direct functional textmining score
rule do_raise_cooccurrence_to_nlp:
    input:
        # Homology corrected co-occurrence scores
        cooccur_corrected = "all_host_parasite_database_cooccur_interact_global_only_parasites_corrected.tsv",
        # NLP scores both ways around
        nlp_scores = rules.do_double_nlp.output.nlp_scores

    output:
        # direct functional textmining score
        # -> Corrected co-occurrence score with additional column with raised score
        cooccur_corrected_raised = "all_host_parasite_database_cooccur_interact_global_only_parasites_corrected_raised.tsv"

    run:
        shell("{workflow.basedir}/scripts/raise_cooccur_to_nlp.py {input.cooccur_corrected} {input.nlp_scores} {output.cooccur_corrected_raised}")


rule do_select_column_textmining:
    input:
        # Corrected co-occurrence score with additional column with raised score
        # "textmining", "taxid1", "protein1", "taxid2", "protein2", "score", "homology_score", "score_corrected", "score_raised"
        cooccur_corrected_raised = rules.do_raise_cooccurrence_to_nlp.output.cooccur_corrected_raised

    output:
        # direct functional textmining score
        functional_textmining_direct = "host_parasite_functional_textmining_direct.tsv"

    run:
        shell("cut -f1,2,3,4,5,9 {input.cooccur_corrected_raised} > {output.functional_textmining_direct}")


# Splitting the parasite-host pairs should be a rule since the input file might change.
rule do_split_parasite_host_pairs:
    input:
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"]

    output:
        parasite_host_split = expand("transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split{split}.tsv", split = PARASITE_HOST_SPLITS)

    run:
        # shuffle the pairs and split into files of size max 50
        shell("""cat {input.parasite_host_pairs} | shuf > transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host_shuffled.tsv""")
        shell("""split -l50 transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host_shuffled.tsv transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split --additional-suffix=.tsv -da 3""")


# Orthology transfer of direct functional textmining scores to obtain transferred functional textmining scores
rule do_transfer_functional_textmining:
    input:
        # Direct functional textmining score
        functional_textmining_direct = rules.do_select_column_textmining.output.functional_textmining_direct,
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"],
        # Species pairs between which to transfer
        transfer_pairs = config["transfer_pairs_path"],
        # Chunks of parasite-host pairs for parallelization
        parasite_host_split = "transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split{split}.tsv",
        # simap rbh
        simap_rbh = config["simap_data_path"] + "rbh/"

    output:
        # Columns: "textmining", "taxid1", "prot1", "taxid2", "prot2", "score"
        # functional_textmining_transferred = expand("host_parasite_functional_textmining_transferred_split{split}.tsv.gz", split = PARASITE_HOST_SPLITS)
        functional_textmining_transferred = "host_parasite_functional_textmining_transferred_split{split}.tsv.gz"

    run: 
        shell("{workflow.basedir}/scripts/orthology_transfer_host_parasite_interactions_generalized.py {input.parasite_host_pairs} {input.transfer_pairs} {input.functional_textmining_direct} {input.simap_rbh} {input.parasite_host_split} 0.041 {output.functional_textmining_transferred}")


# Aggregate transferred functional textmining scores
rule do_aggregate_functional_textmining:
    input: 
        # Columns: "textmining", "taxid1", "prot1", "taxid2", "prot2", "score"
        functional_textmining_transferred = expand("host_parasite_functional_textmining_transferred_split{split}.tsv.gz", split = PARASITE_HOST_SPLITS)
        # for testing
        # functional_textmining_transferred = expand("host_parasite_functional_textmining_transferred_split{split}.tsv.gz", split = ["000"])
    output:
        functional_textmining_transferred_aggregated = "host_parasite_functional_textmining_transferred_aggregated.tsv"
    run:
        # combine all chunks from previous process
        shell("zcat host_parasite_functional_textmining_transferred_split*.tsv.gz > host_parasite_functional_textmining_transferred.tsv")
        # sort scores by taxid and protein
        shell("sort -S30G --parallel=20 -k2,5 host_parasite_functional_textmining_transferred.tsv > host_parasite_functional_textmining_transferred_sorted.tsv")
        # aggregate scores from transfer so that we have one score per PPI
        shell("{workflow.basedir}/scripts/orthology_score_aggregation_host_parasite_interactions.py host_parasite_functional_textmining_transferred_sorted.tsv {output.functional_textmining_transferred_aggregated} 0.041")
        # clean up
        shell("rm host_parasite_functional_textmining_transferred_sorted.tsv")


# Orthology transfer of direct NLP score (= direct physical textmining score) scores to obtain transferred physical textmining scores
rule do_transfer_physical_textmining:
    input:
        # Direct physical textmining score
        nlp_scores = rules.do_double_nlp.output.nlp_scores,
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"],
        # Species pairs between which to transfer
        transfer_pairs = config["transfer_pairs_path"],
        # Chunks of parasite-host pairs for parallelization
        parasite_host_split = "transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split{split}.tsv",
        # simap rbh
        simap_rbh = config["simap_data_path"] + "rbh/"

    output:
        # Columns: "textmining", "taxid1", "prot1", "taxid2", "prot2", "score"
        physical_textmining_transferred = "host_parasite_physical_textmining_transferred_split{split}.tsv.gz"

    run:
        shell("{workflow.basedir}/scripts/orthology_transfer_host_parasite_interactions_generalized.py {input.parasite_host_pairs} {input.transfer_pairs} {input.nlp_scores} {input.simap_rbh} {input.parasite_host_split} 0.006 {output.physical_textmining_transferred}")


# Aggregate transferred physical textmining scores
rule do_aggregate_physical_textmining:
    input: 
        # Columns: "textmining", "taxid1", "prot1", "taxid2", "prot2", "score"
        physical_textmining_transferred = expand("host_parasite_physical_textmining_transferred_split{split}.tsv.gz", split = PARASITE_HOST_SPLITS)
    output:
        physical_textmining_transferred_aggregated = "host_parasite_physical_textmining_transferred_aggregated.tsv"
    run:
        # combine all chunks from previous process
        shell("zcat host_parasite_physical_textmining_transferred_split*.tsv.gz > host_parasite_physical_textmining_transferred.tsv")
        # sort scores by taxid and protein
        shell("sort -S30G --parallel=20 -k2,5 host_parasite_physical_textmining_transferred.tsv > host_parasite_physical_textmining_transferred_sorted.tsv")
        # aggregate scores from transfer so that we have one score per PPI
        shell("{workflow.basedir}/scripts/orthology_score_aggregation_host_parasite_interactions.py host_parasite_physical_textmining_transferred_sorted.tsv {output.physical_textmining_transferred_aggregated} 0.006")
        # clean up
        shell("rm host_parasite_physical_textmining_transferred_sorted.tsv")


# Combine direct and transferred functional textmining scores to obtain final functional score from textmining
rule do_combine_functional_textmining:
    input: 
        functional_textmining_transferred_aggregated = rules.do_aggregate_functional_textmining.output.functional_textmining_transferred_aggregated,
        functional_textmining_direct = rules.do_select_column_textmining.output.functional_textmining_direct
    output:
        functional_textmining = "host_parasite_functional_textmining.tsv"
    run:
        shell("{workflow.basedir}/scripts/combine_direct_transferred.py {input.functional_textmining_direct} {input.functional_textmining_transferred_aggregated} 0.041 {output.functional_textmining}")

# Combine direct and transferred physical textmining scores to obtain final physical score from textmining
rule do_combine_physical_textmining:
    input: 
        physical_textmining_transferred_aggregated = rules.do_aggregate_physical_textmining.output.physical_textmining_transferred_aggregated,
        nlp_scores = rules.do_double_nlp.output.nlp_scores
    output:
        physical_textmining = "host_parasite_physical_textmining.tsv"
    run:
        shell("{workflow.basedir}/scripts/combine_direct_transferred.py {input.nlp_scores} {input.physical_textmining_transferred_aggregated} 0.006 {output.physical_textmining}")

###################################
####### Experiments channel #######
###################################

# Filter functional experiments data for host-parasite interactions
rule do_filter_functional_experiments:
    input:
        functional_experiments_unfiltered = config["experiments_functional_path"],
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"]
    output:
        functional_experiments_direct = "host_parasite_functional_experiments_direct.tsv"

    run:
        shell("""
        awk -F"\\t" 'NR==FNR{{c[$1$2]++;next}};c[$2$4] > 0 || c[$4$2] > 0' {input.parasite_host_pairs} {input.functional_experiments_unfiltered} > {output.functional_experiments_direct}
        """)


# Filter physical experiments data for host-parasite interactions
rule do_filter_physical_experiments:
    input:
        physical_experiments_unfiltered = config["experiments_physical_path"],
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"]
    output:
        physical_experiments_direct = "host_parasite_physical_experiments_direct.tsv"

    run:
        shell("""
        awk -F"\\t" 'NR==FNR{{c[$1$2]++;next}};c[$2$4] > 0 || c[$4$2] > 0' {input.parasite_host_pairs} {input.physical_experiments_unfiltered} > {output.physical_experiments_direct}
        """)


# Orthology transfer of direct functional experiments scores to obtain transferred functional experiments scores
rule do_transfer_functional_experiments:
    input:
        # Direct functional experiments score
        functional_experiments_direct = rules.do_filter_functional_experiments.output.functional_experiments_direct,
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"],
        # Species pairs between which to transfer
        transfer_pairs = config["transfer_pairs_path"],
        # Chunks of parasite-host pairs for parallelization
        parasite_host_split = "transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split{split}.tsv",
        # simap rbh
        simap_rbh = config["simap_data_path"] + "rbh/"

    output:
        # Columns: "experiments", "taxid1", "prot1", "taxid2", "prot2", "score"
        functional_experiments_transferred = "host_parasite_functional_experiments_transferred_split{split}.tsv.gz"

    run: 
        shell("{workflow.basedir}/scripts/orthology_transfer_host_parasite_interactions_generalized.py {input.parasite_host_pairs} {input.transfer_pairs} {input.functional_experiments_direct} {input.simap_rbh} {input.parasite_host_split} 0.041 {output.functional_experiments_transferred}")


# Aggregate transferred functional experiments scores
rule do_aggregate_functional_experiments:
    input: 
        # Columns: "experiments", "taxid1", "prot1", "taxid2", "prot2", "score"
        functional_experiments_transferred = expand("host_parasite_functional_experiments_transferred_split{split}.tsv.gz", split = PARASITE_HOST_SPLITS)
    output:
        functional_experiments_transferred_aggregated = "host_parasite_functional_experiments_transferred_aggregated.tsv"
    run:
        # combine all chunks from previous process
        shell("zcat host_parasite_functional_experiments_transferred_split*.tsv.gz > host_parasite_functional_experiments_transferred.tsv")
        # sort scores by taxid and protein
        shell("sort -S30G --parallel=20 -k2,5 host_parasite_functional_experiments_transferred.tsv > host_parasite_functional_experiments_transferred_sorted.tsv")
        # aggregate scores from transfer so that we have one score per PPI
        shell("{workflow.basedir}/scripts/orthology_score_aggregation_host_parasite_interactions.py host_parasite_functional_experiments_transferred_sorted.tsv {output.functional_experiments_transferred_aggregated} 0.041")
        # clean up
        shell("rm host_parasite_functional_experiments_transferred_sorted.tsv")


# Orthology transfer of direct physical experiments scores to obtain transferred physical experiments scores
rule do_transfer_physical_experiments:
    input:
        # Direct physical experiments score
        physical_experiments_direct = rules.do_filter_physical_experiments.output.physical_experiments_direct,
        # Host parasite taxid pairs 
        parasite_host_pairs = config["parasite_host_pairs_path"],
        # Species pairs between which to transfer
        transfer_pairs = config["transfer_pairs_path"],
        # Chunks of parasite-host pairs for parallelization
        parasite_host_split = "transfer_species_pairs_parasite_host/transfer_species_pairs_parasite_host.split{split}.tsv",
        # simap rbh
        simap_rbh = config["simap_data_path"] + "rbh/"

    output:
        # Columns: "experiments", "taxid1", "prot1", "taxid2", "prot2", "score"
        physical_experiments_transferred = "host_parasite_physical_experiments_transferred_split{split}.tsv.gz"

    run: 
        shell("{workflow.basedir}/scripts/orthology_transfer_host_parasite_interactions_generalized.py {input.parasite_host_pairs} {input.transfer_pairs} {input.physical_experiments_direct} {input.simap_rbh} {input.parasite_host_split} 0.006 {output.physical_experiments_transferred}")


# Aggregate transferred physical experiments scores
rule do_aggregate_physical_experiments:
    input: 
        # Columns: "experiments", "taxid1", "prot1", "taxid2", "prot2", "score"
        physical_experiments_transferred = expand("host_parasite_physical_experiments_transferred_split{split}.tsv.gz", split = PARASITE_HOST_SPLITS)
    output:
        physical_experiments_transferred_aggregated = "host_parasite_physical_experiments_transferred_aggregated.tsv"
    run:
        # combine all chunks from previous process
        shell("zcat host_parasite_physical_experiments_transferred_split*.tsv.gz > host_parasite_physical_experiments_transferred.tsv")
        # sort scores by taxid and protein
        shell("sort -S30G --parallel=20 -k2,5 host_parasite_physical_experiments_transferred.tsv > host_parasite_physical_experiments_transferred_sorted.tsv")
        # aggregate scores from transfer so that we have one score per PPI
        shell("{workflow.basedir}/scripts/orthology_score_aggregation_host_parasite_interactions.py host_parasite_physical_experiments_transferred_sorted.tsv {output.physical_experiments_transferred_aggregated} 0.006")
        # clean up
        shell("rm host_parasite_physical_experiments_transferred_sorted.tsv")

# Combine direct and transferred physical experiments scores to obtain final physical score from experiments
rule do_combine_functional_experiments:
    input: 
        functional_experiments_transferred_aggregated = rules.do_aggregate_functional_experiments.output.functional_experiments_transferred_aggregated,
        functional_experiments_direct = rules.do_filter_functional_experiments.output.functional_experiments_direct
    output:
        functional_experiments = "host_parasite_functional_experiments.tsv"
    run:
        shell("{workflow.basedir}/scripts/combine_direct_transferred.py {input.functional_experiments_direct} {input.functional_experiments_transferred_aggregated} 0.041 {output.functional_experiments}")

# Combine direct and transferred physical experiments scores to obtain final physical score from experiments
rule do_combine_physical_experiments:
    input: 
        physical_experiments_transferred_aggregated = rules.do_aggregate_physical_experiments.output.physical_experiments_transferred_aggregated,
        physical_experiments_direct = rules.do_filter_physical_experiments.output.physical_experiments_direct
    output:
        physical_experiments = "host_parasite_physical_experiments.tsv"
    run:
        shell("{workflow.basedir}/scripts/combine_direct_transferred.py {input.physical_experiments_direct} {input.physical_experiments_transferred_aggregated} 0.006 {output.physical_experiments}")


# Combine functional experiments and textmining scores to obtain final functional score
rule do_combine_functional_channels:
    input: 
        functional_experiments = rules.do_combine_functional_experiments.output.functional_experiments,
        functional_textmining = rules.do_combine_functional_textmining.output.functional_textmining
    output:
        functional = "host_parasite_functional.tsv"
    run:
        shell("cat {input.functional_experiments} {input.functional_textmining} > host_parasite_functional_concat.tsv")
        shell("{workflow.basedir}/scripts/combine_channels.py host_parasite_functional_concat.tsv 0.041 {output.functional}")
        shell("rm host_parasite_functional_concat.tsv")


# Combine physical experiments and textmining scores to obtain final physical score
rule do_combine_physical_channels:
    input: 
        physical_experiments = rules.do_combine_physical_experiments.output.physical_experiments,
        physical_textmining = rules.do_combine_physical_textmining.output.physical_textmining
    output:
        physical = "host_parasite_physical.tsv"
    run:
        shell("cat {input.physical_experiments} {input.physical_textmining} > host_parasite_physical_concat.tsv")
        shell("{workflow.basedir}/scripts/combine_channels.py host_parasite_physical_concat.tsv 0.006 {output.physical}")
        shell("rm host_parasite_physical_concat.tsv")

# rule to run all rules
rule all:
    input: 
        functional = rules.do_combine_functional_channels.output.functional,
        physical = rules.do_combine_physical_channels.output.physical
