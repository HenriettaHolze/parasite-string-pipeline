#!/usr/bin/env python
import pandas as pd
import sys
from collections import defaultdict

# PSEUDO CODE
# transfer inter-species PPI interactions to orthologs of interacting species
# for each PPI:
#   transfer to other host:
#      for each target host:
#           find ortholog(s) (same score hits)
#           for each ortholog:
#               calculate new score, add interaction
#   transfer to other parasite:
#      for each target parasite:
#           find ortholog(s) (same score hits)
#           for each ortholog:
#               calculate new score, add interaction
#   transfer to other host-parasite pairs:
#       for each parasite-host pair:
#           find ortholog(s) (same score hits)
#           for each ortholog in host:
#               for each ortholog in parasite:
#                   calculate new score, add interaction



def orthology_transfer(parasite_host_pairs_path, transfer_pairs_path, evidence_path, parasite_host_pairs_split_path, simap_rbh_path, PRIOR):

    # transfer inter-species PPIs 
    alpha = 0.2

    # dict and set of host_parasite_pairs for fast look up
    parasite_host_pairs = pd.read_csv(parasite_host_pairs_path, 
                                    sep="\t", 
                                    header=None, 
                                    names=["parasite_taxid", "host_taxid"], 
                                    dtype={"parasite_taxid": str, "host_taxid":str})

    parasite_host_pairs_split = pd.read_csv(parasite_host_pairs_split_path, 
                                    sep="\t", 
                                    header=None, 
                                    names=["parasite_taxid", "host_taxid"], 
                                    dtype={"parasite_taxid": str, "host_taxid":str})

    # dict of transfer_pairs
    transfer_pairs = pd.read_csv(transfer_pairs_path,
                                sep="\t",
                                header=None,
                                names=["class", "taxid1", "taxid2", "level"], 
                                dtype={"taxid1": str, "taxid2":str})

    # load the evidence 
    evidence = pd.read_csv(evidence_path,
                                        sep="\t",
                                        header=None,
                                        names = ["channel", "taxid1", "protein1", "taxid2", "protein2", "score"],
                                        dtype={"taxid1": str, "taxid2":str})


    ############
    # identify for each host-parasite pair to which host-parasite pairs I can transfer (when transferring both host and parasite)

    # dictionary of parasites between which I can transfer
    pp_pairs_dict = defaultdict(list)
    for row in transfer_pairs.loc[transfer_pairs["class"] == "parasite"].iterrows():
        pp_pairs_dict[row[1]["taxid1"]].append(row[1]["taxid2"])

    # list of host pairs between which I can transfer
    # this is the longest list so we want to do the final look up in here
    hh_pairs_set = {
        (row[1]["taxid1"], row[1]["taxid2"])
        for row in transfer_pairs.loc[transfer_pairs["class"] == "host"].iterrows()
    }
    hh_pairs_dict = defaultdict(list)
    for row in transfer_pairs.loc[transfer_pairs["class"] == "host"].iterrows():
        hh_pairs_dict[row[1]["taxid1"]].append(row[1]["taxid2"])


    # dictionary which parasites infect which hosts
    parasite_host_dict = defaultdict(list)
    for row in parasite_host_pairs.iterrows():
        parasite_host_dict[row[1]["parasite_taxid"]].append(row[1]["host_taxid"])

    host_parasite_dict = defaultdict(list)
    for row in parasite_host_pairs.iterrows():
        host_parasite_dict[row[1]["host_taxid"]].append(row[1]["parasite_taxid"])

    target_dict = {}
    for row in parasite_host_pairs.iterrows():
        target_hp = []
        parasite = row[1]["parasite_taxid"]
        host = row[1]["host_taxid"]

        target_ps = pp_pairs_dict[parasite]

        for target_p in target_ps:
            target_hs = parasite_host_dict[target_p]

            for target_h in target_hs:
                if (host, target_h) in hh_pairs_set:
                    target_hp.append((target_h, target_p))

        if len(target_hp) != 0:
            target_dict[host +"_"+ parasite] = target_hp

    # only 559 host-parasite pairs that can have transfer via host and parasite

    transferred_ppis = []

    # iterate over host-parasite pairs
    for row in parasite_host_pairs_split.iterrows():
        parasite = row[1]["parasite_taxid"]
        host = row[1]["host_taxid"]

        print(host, parasite)

        # select evidence from host-parasite pair
        evidence_host_parasite = evidence.loc[(evidence["taxid1"] == host) & (evidence["taxid2"] == parasite)]

        if len(evidence_host_parasite) == 0:
            continue

        # get all hosts and parasite between which we will transfer
        # for transferring either parasite or host
        target_hosts = set(hh_pairs_dict[host]).intersection(set(parasite_host_dict[parasite]))
        target_parasites = set(pp_pairs_dict[parasite]).intersection(set(host_parasite_dict[host]))
        # for transferring on both sides
        target_hp = target_dict.get(host +"_"+ parasite, [])
        target_hs = []
        target_ps = []
        for target_hp_pair in target_hp:
            target_hs.append(target_hp_pair[0])
            target_ps.append(target_hp_pair[1])
        target_ps = set(target_ps)
        target_hs = set(target_hs)

        # load orthology data of relevant organisms in memory
        rbh_dict = dict()

        for target_host in target_hs.union(target_hosts):
            rbh = pd.read_csv(simap_rbh_path + "/{}.{}_rbh.tsv".format(host, target_host), sep="\t", header=None, names=["qseqid", "sseqid", "normbitscore", "percident", "percsim", "qstart", "qend", "sstart", "send"])
            # protein IDs of evidence does not contain taxid
            rbh['qseqid'].replace({r"^[0-9]+\.": ''}, inplace = True, regex = True)
            rbh['sseqid'].replace({r"^[0-9]+\.": ''}, inplace = True, regex = True)

            rbh_dict[host + "_" + target_host] = defaultdict(list)
            for row in rbh.iterrows():
                # possibly multiple orthologs (same score)
                rbh_dict[host + "_" + target_host][row[1]["qseqid"]].append((row[1]["sseqid"], row[1]["normbitscore"]))

        for target_parasite in target_ps.union(target_parasites):
            rbh = pd.read_csv(simap_rbh_path + "/{}.{}_rbh.tsv".format(parasite, target_parasite), sep="\t", header=None, names=["qseqid", "sseqid", "normbitscore", "percident", "percsim", "qstart", "qend", "sstart", "send"])
            # protein IDs of evidence does not contain taxid
            rbh['qseqid'].replace({r"^[0-9]+\.": ''}, inplace = True, regex = True)
            rbh['sseqid'].replace({r"^[0-9]+\.": ''}, inplace = True, regex = True)

            rbh_dict[parasite + "_" + target_parasite] = defaultdict(list)
            for row in rbh.iterrows():
                # possibly multiple orthologs (same score)
                rbh_dict[parasite + "_" + target_parasite][row[1]["qseqid"]].append((row[1]["sseqid"], row[1]["normbitscore"]))

        # change default return value
        for key in rbh_dict.keys():
            rbh_dict[key].default_factory = lambda:False

        # iterate over evidence from host-parasite pair
        for ppi in evidence_host_parasite.iterrows():
            ppi = ppi[1]
            host_prot, parasite_prot, score = (ppi["protein1"], ppi["protein2"], ppi["score"])

            # one-way transfer
            for target_host in target_hosts:
                # check if protein has rbh in other host
                target_proteins = rbh_dict[host + "_" + target_host][host_prot]
                if target_proteins:
                    for target_prot in target_proteins:
                        new_score = (score - PRIOR) * (target_prot[1]) ** alpha + PRIOR
                        transferred_ppis.append((target_host, target_prot[0], parasite, parasite_prot, new_score))

            for target_parasite in target_parasites:
                # potentially multiple rbh's (identical sequences)
                target_proteins = rbh_dict[parasite + "_" + target_parasite][parasite_prot]
                if target_proteins:
                    for target_prot in target_proteins:                       
                        new_score = (score - PRIOR) * (target_prot[1]) ** alpha + PRIOR
                        transferred_ppis.append((host, host_prot, target_parasite, target_prot[0], new_score))

            for target_hp_pair in target_hp:
                target_proteins_host = rbh_dict[host + "_" + target_hp_pair[0]][host_prot]
                target_proteins_parasite = rbh_dict[parasite + "_" + target_hp_pair[1]][parasite_prot]
                if target_proteins_host and target_proteins_parasite:
                    # iterate over rbh's in host and parasite
                    for target_prot_host in target_proteins_host:
                        for target_prot_parasite in target_proteins_parasite:
                            new_score = (score - PRIOR) * (target_prot_host[1] * target_prot_parasite[1]) ** alpha + PRIOR
                            transferred_ppis.append((target_hp_pair[0], target_prot_host[0], target_hp_pair[1], target_prot_parasite[0], new_score))

    transferred_ppis_df = pd.DataFrame(transferred_ppis, columns = ["taxid1", "prot1", "taxid2", "prot2", "score"])
    transferred_ppis_df["channel"] = evidence["channel"][0]
    transferred_ppis_df = transferred_ppis_df[["channel", "taxid1", "prot1", "taxid2", "prot2", "score"]]

    return(transferred_ppis_df)


if __name__ == "__main__":

    parasite_host_pairs_path = sys.argv[1]
    transfer_pairs_path = sys.argv[2]
    evidence_path = sys.argv[3]
    simap_rbh_path = sys.argv[4]
    parasite_host_pairs_split_path = sys.argv[5]
    PRIOR = float(sys.argv[6])
    output_file = sys.argv[7]

    transferred_ppis_df = orthology_transfer(parasite_host_pairs_path, transfer_pairs_path, evidence_path, parasite_host_pairs_split_path, simap_rbh_path, PRIOR)
    
    transferred_ppis_df.to_csv(output_file, sep="\t", compression="gzip", index=False, header=False)
