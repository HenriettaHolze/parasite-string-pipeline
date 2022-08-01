#!/usr/bin/env python
import sys

# functional textmining score is max(cooccurrence_score, theoretical_functional_score)
# with (1-theoretical_functional_score) = (1-physical_score) / (1-0.006) * (1 - 0.041)
# we cannot simply raise the cooccurrence score to the NLP score because of the different priors


def make_cooccur_dict(input_cooccur_corrected):
    # read in big file with cooccurrence scores as dict
    cooccur_dict = {}
    with open(input_cooccur_corrected, "r") as f:
        for line in f:
            _, taxid1, prot1, taxid2, prot2, score, homology_score, score_corrected = line.strip().split("\t")
            cooccur_dict["..".join([taxid1, prot1, taxid2, prot2])] = [float(score), float(homology_score), float(score_corrected), float(score_corrected)]
    return cooccur_dict


def raise_scores(input_nlp_scores_doubled, cooccur_dict):
    # iterate over smaller file with NLP scores
    with open(input_nlp_scores_doubled, "r") as f:
        for line in f:
            _, taxid1, prot1, taxid2, prot2, nlp_score = line.strip().split("\t")
            nlp_score = float(nlp_score)
            dict_key = "..".join([taxid1, prot1, taxid2, prot2])
            # if physical PPI has functional PPI and if score is lower, update score
            if dict_key in cooccur_dict:
                theoretical_functional_score = 1 - (1 - nlp_score) / (1 - 0.006) * (1 - 0.041)
                cooccur_dict[dict_key][3] = max(cooccur_dict[dict_key][3], theoretical_functional_score)
            
            # if physical PPI does not have functional PPI (should not happen)
            if dict_key not in cooccur_dict:
                Warning("Physical PPI that does not have functional score:", dict_key.split(".."))
                cooccur_dict[dict_key] = [0, 0, 0, nlp_score]
    
    return cooccur_dict


def main():
    input_cooccur_corrected = sys.argv[1]
    input_nlp_scores_doubled = sys.argv[2]
    output_file = open(sys.argv[3], "w")

    cooccur_dict = make_cooccur_dict(input_cooccur_corrected)
    cooccur_dict = raise_scores(input_nlp_scores_doubled, cooccur_dict)

    # write dict with raised scores to file
    for key_id in cooccur_dict.keys():
        output_file.write("\t".join(["textminig"] + key_id.split("..") + [str(s) for s in cooccur_dict[key_id]]) + "\n")

    output_file.close()

if __name__ == "__main__":
    main()
