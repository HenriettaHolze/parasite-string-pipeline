#!/usr/bin/python3

######### DESCRIPTION #########
# Combine STRING scores from multiple channels for each protein pair
# Adapted from combine_subscores.py (https://string-db.org/download/combine_subscores.py)
# The memory scales with the input file size, run on C2. 
# 
######### INPUT ###############
# A concatenated file with inter-species PPIs from multiple channels with colums [channel, taxid1, 
# protein1, taxid2, protein2, score], the prior and the location of the output file. 
# 
######### OUTPUT ##############
# A file with aggregated STRING scores with columns [taxid1, protein1, taxid2, protein2, score]

import sys

def compute_prior_away(score, prior):

    if score < prior: score = prior
    score_no_prior = (score - prior) / (1 - prior)

    return score_no_prior

def collect_score_products(input_file, prior):
    '''For each protein pair in scores file, pre-process and multiply all individual scores. 
    Return a dictionary with (taxid1, protein1, taxid2, protein2), product of processed individual scores.'''
    score_products_one_minus = {}

    with open(input_file, 'r') as i:
        for line in i:
            (_, taxid1, protein1, taxid2, protein2, score) = line.strip('\n').split('\t')

            # remove prior from individual score, do 1 - score conversion
            score_no_prior = compute_prior_away(float(score), prior)
            score_one_minus = 1 - score_no_prior

            # multiply all individual scores of each protein pair
            if (taxid1, protein1, taxid2, protein2) not in score_products_one_minus:
                score_products_one_minus[(taxid1, protein1, taxid2, protein2)] = score_one_minus
            else:
                score_products_one_minus[(taxid1, protein1, taxid2, protein2)] *= score_one_minus

    return score_products_one_minus

def calc_combined_scores(score_products_one_minus, prior, output_file):

    with open(output_file, "w") as o:
        '''Calculate the final aggregated score from pre-processed individual scores and print to stdout.'''
        for taxid_strings, score_product_one_minus in score_products_one_minus.items():
            # do 1-score conversion
            score_combined = 1 - score_product_one_minus
            # scale down
            score_combined *= (1 - prior)
            # add prior back once
            score_combined += prior
            # round
            score_combined_rounded = round(score_combined, 3)
            o.write('\t'.join(taxid_strings) + "\t" + str(score_combined_rounded) + "\n")


if __name__ == '__main__':

    input_file = sys.argv[1]
    prior = float(sys.argv[2])
    output_file = sys.argv[3]

    score_products_one_minus = collect_score_products(input_file, prior)
    calc_combined_scores(score_products_one_minus, prior, output_file)
