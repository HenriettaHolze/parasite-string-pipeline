#!/usr/bin/env python
import sys

# combine direct and transferred scores
# iterate over file with transferred score and check if a direct score exists
# combine scores if necessary and print to output file
# The memory scales with the file size of the direct evidence and runtime with the 
# size of the transferred evidence. 

# I don't write the PPIs with only direct evidence...

def make_direct_dict(input_direct):
    # read in big file with cooccurrence scores as dict
    direct_dict = {}
    with open(input_direct, "r") as f:
        for line in f:
            _, taxid1, prot1, taxid2, prot2, direct_score = line.strip().split("\t")
            direct_dict["..".join([taxid1, prot1, taxid2, prot2])] = float(direct_score)
    return direct_dict


def compute_prior_away(score, prior):
    if score < prior: score = prior
    score_no_prior = (score - prior) / (1 - prior)
    return score_no_prior


def combine_scores(direct_score, transferred_score, prior):
    both_prior_corrected = 1 - (1 - compute_prior_away(direct_score, prior)) * (1 - compute_prior_away(transferred_score, prior))
    combined_score = both_prior_corrected * (1 - prior) + prior
    return round(combined_score, 3)

def main():
    input_direct = sys.argv[1]
    input_transferred = sys.argv[2]
    PRIOR = float(sys.argv[3])
    output_file = open(sys.argv[4], "w")

    direct_dict = make_direct_dict(input_direct)

    combined_direct_ppis = []

    with open(input_transferred, "r") as f:
        for line in f:
            channel, taxid1, prot1, taxid2, prot2, transferred_score = line.strip().split("\t")
            transferred_score = float(transferred_score)
            dict_key = "..".join([taxid1, prot1, taxid2, prot2])

            if dict_key in direct_dict:
                # calculate combined score
                combined_score = combine_scores(direct_score = direct_dict[dict_key], transferred_score = transferred_score, prior = PRIOR)
                # write PPI with both direct and transferred score
                output_file.write("\t".join([channel, taxid1, prot1, taxid2, prot2, str(combined_score)]) + "\n")
                # keep track of which direct scores were already added
                combined_direct_ppis.append(dict_key)

            else:
                # write PPI with only transferred score
                output_file.write("\t".join([channel, taxid1, prot1, taxid2, prot2, str(transferred_score)]) + "\n")

    # write PPIs with only direct score
    combined_direct_ppis = set(combined_direct_ppis)
    for ppi in direct_dict.keys():
        if ppi not in combined_direct_ppis:
            output_file.write("\t".join([channel] + ppi.split("..") + [str(direct_dict[ppi])]) + "\n")

    output_file.close()

if __name__ == "__main__":
    main()
    
