#!/usr/bin/env python
import sys

input_file = sys.argv[1]
output_file = open(sys.argv[2], "w")
PRIOR = float(sys.argv[3])

channel_old, taxid1_old, prot1_old, taxid2_old, prot2_old, score_old = None, None, None, None, None, None
running_prod = None

# reading the file line by line avoids loading it into buffer
for line in open(input_file):

    channel, taxid1, prot1, taxid2, prot2, score = line.split("\t")
    score = float(score)

    # new PPI
    if (taxid1, prot1, taxid2, prot2) != (taxid1_old, prot1_old, taxid2_old, prot2_old):
        # finish previous score and write to file (except for first line)
        if running_prod != None:
            final_score = 1 - (1 - PRIOR) * running_prod
            # write both ways
            output_file.write('\t'.join([channel, taxid1_old, prot1_old, taxid2_old, prot2_old, str(final_score)]) + "\n")
            output_file.write('\t'.join([channel, taxid2_old, prot2_old, taxid1_old, prot1_old, str(final_score)]) + "\n")

        # start new running product
        running_prod = (1 - score) / (1 - PRIOR)

    else:
        running_prod *= (1 - score) / (1 - PRIOR)

    taxid1_old, prot1_old, taxid2_old, prot2_old = taxid1, prot1, taxid2, prot2

# end of file
final_score = 1 - (1 - PRIOR) * running_prod
# write both ways to be consistent with other files
output_file.write('\t'.join([channel, taxid1, prot1, taxid2, prot2, str(final_score)]) + "\n")
output_file.write('\t'.join([channel, taxid2, prot2, taxid1, prot1, str(final_score)]) + "\n")

output_file.close()
