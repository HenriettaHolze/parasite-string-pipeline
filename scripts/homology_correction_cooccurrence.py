#!/usr/bin/env python
import math
import pandas as pd
from datetime import datetime
import sys

##############
### DESCRIPTION:
# Calculate homology score and corrected co-occurrence score for all PPIs with cooccurrence score.
# Formula is from Damian, based on self-normalized bitscore from SIMAP. 
# 
### RUN INFO:
# Is called from snakemake
# 
### INPUT: 
# - sorted co-occurrence scores
# - simap-processed/selfnorm-bitscores/ (self normalized bitscores from SIMAP all vs. all searches, in both directions)
# 
### OUTPUT:
# - host_parasite_func_textmining_direct_transferred_homology.tsv (same as input but with homology score as additional column)
##############

# For the snakemake pipeline the input and output filenames are parsed by the pipeline

# self normalized bit scores
# dir_bitscores = "simap-processed/selfnorm-bitscores/"
dir_bitscores = sys.argv[1]
# sort file by taxids so you only have to read in one SIMAP file at a time
# input_file = "host_parasite_func_textmining_direct_transferred_sorted.tsv"
input_file = sys.argv[2]
# output_file = open("host_parasite_func_textmining_direct_transferred_homology.tsv", "w")
output_file = open(sys.argv[3], "w")

# This formula is from Damian and has been calibrated for STRING (intra-species interactions)
def calc_homology_score(x):
    return 0.537868 + (0.999993 - 0.537868) / (1 + math.exp(2.32771 * (-1.58792 - math.log(x))))

# homology correction on score w/o prior
prior = 0.041
def compute_prior_away(score, prior):

    # there should be no scores below the prior in the input file
    if score < prior: score = prior
    score_no_prior = (score - prior) / (1 - prior)

    return score_no_prior


taxid1_old, taxid2_old = None, None

# reading the file line by line avoids loading it into buffer
for line in open(input_file):

    cooccur, taxid1, prot1, taxid2, prot2, score = line.strip().split("\t")

    # new PPI
    if (taxid1, taxid2) != (taxid1_old, taxid2_old):
        # check how long reading the file takes by printing time before and after

        print(taxid1, taxid2)
        # read in simap file
        simap_result = pd.read_csv(dir_bitscores + "/{}.{}_selfnorm_bitscore.tsv.gz".format(taxid1, taxid2), 
                                   sep="\t", 
                                   header=None, 
                                   names=["qseqid", "sseqid", "normbitscore", "percident", "percsim", "qstart", "qend", "sstart", "send"], 
                                   usecols = ["qseqid", "sseqid", "normbitscore"], 
                                   compression="gzip")

        # set qseqid and sseqid as index for faster look up
        simap_result = simap_result.set_index(["qseqid", "sseqid"])
    
    try:
        # look up using index
        # this will throw an error if the index does not exist. Thus try except statement
        norm_bitscore = simap_result.loc[(".".join([taxid1, prot1]), ".".join([taxid2, prot2])), "normbitscore"]
        # norm_bitscore = simap_result.loc[(simap_result["qseqid"] == ".".join(taxid1, prot1)) & (simap_result["sseqid"] == ".".join(taxid2, prot2)), "normbitscore"]

        # many PPIs will not have a hit and will not be corrected at all --> homology score of 0

        homology_score = calc_homology_score(norm_bitscore)

        score_no_prior = compute_prior_away(score = float(score), prior = prior)

        score_corrected = score_no_prior * (1 - homology_score)

        score_corrected = score_corrected * (1.0 - prior) + prior


    except KeyError:
        homology_score = 0
        score_corrected = score

    except Exception as e:
        print("Error")
        print(e)
        break

    finally:
        output_file.write('\t'.join([cooccur, taxid1, prot1, taxid2, prot2, score, str(homology_score), str(score_corrected)]) + "\n")

        taxid1_old, taxid2_old = taxid1, taxid2

output_file.close()

