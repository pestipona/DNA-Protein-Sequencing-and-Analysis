import protein_tests as test
import json

### PHASE 1 ###

def read_file(filename):
    # open the file
    f = open(filename, 'r')
    result = ""
    for line in f:
        line = line.rstrip('\n')
        result += line
    f.close()
    return result

def dna_to_rna(dna, start_index):
    dna = dna[start_index:]         # slice dna based on starting index
    dna = dna.replace('T', 'U')     # replace each T base to U base
    codon = ""                      # create an empty codon
    codon_list = []                 # create an empty codon list

    # loop through entire dna sequence
    for char in dna:
        if (len(codon) < 3):            # combine 3 bases to a codon
            codon += char               # build codon with chars

        if (len(codon) == 3):
            codon_list.append(codon)    # append codon to codon_list
            if (codon == "UAA" or       # check for stop codons
                codon == "UAG" or
                codon == "UGA"):
                return codon_list       # return list with stop codon
            codon = ""  # reset codon string

    return codon_list   # return list with no stop codon

def make_codon_dictionary(filename):
    # open the JSON file
    f = open(filename, 'r')

    # return JSON object as a dictionary
    dict_old = json.load(f)

    # create a new empty dict to store codons to amino acids
    dict_new = dict()

    # iterate through the old dictionary
    for key in dict_old:
        # iterate through each codon in the list
        for codon in dict_old[key]:
            # replace all T's to U's in the Codon
            codon = codon.replace('T', 'U')
            # add codon & amino-acid value pair to new dict
            dict_new[codon] = key

    # return new dictionary
    return dict_new

def generate_protein(codons, codon_dict):
    # create empty protein list to hold amino-acids
    protein = []

    # iterate explicitly through the RNA list
    for i in range(len(codons)):
        # check first codon in the list for special-case
        if (i == 0):
            # if first codon is "AUG"
            if (codons[i] == "AUG"):
                # Start of protein sequence
                protein.append("Start")
            else:
                # if not AUG use amino acid (Met)
                protein.append(codon_dict[codons[i]])
        else:
            # add associated amino acid (based on codon_dict) to new list.
            protein.append(codon_dict[codons[i]])

    # return protein sequence as a list
    return protein

def synthesize_proteins(dna_filename, codon_filename):
    # read the dna file
    dna = read_file(dna_filename)

    # create codon dictionary
    codon_dictionary = make_codon_dictionary(codon_filename)

    # create empty protein list to protein sequence
    proteins = []

    # creat place value holders
    index = 0
    num_codons = 0
    num_bases = len(dna)

    # iterate through the dna string
    while(index < len(dna)):
        # check if codon is ATG
        if dna[index:index+3] == "ATG":
            # create rna list based on dna
            rna = dna_to_rna(dna, index)
            # get number of codons in the rna list
            num_codons = num_codons + len(rna)
            # update new index value
            index = index + len(rna)*3
            # create protein list based on rna and dict
            protein = generate_protein(rna, codon_dictionary)
            # append to proteins list
            proteins.append(protein)
        else:
            index = index + 1

    unused_bases = num_bases - num_codons*3
    syn_proteins = len(proteins)

    #print("There are", num_bases, "total bases,")
    #print(unused_bases, "unused bases,", "and", syn_proteins, "synthesized proteins")

    # return list of protein sequences
    return proteins

### PHASE 2 ###

def common_proteins(protein_list1, protein_list2):
    # create blank list for storing common proteins
    common = []

    # compare both input list1 and list2
    for i in protein_list1:
        # go through each item in list2
        for j in protein_list2:
            # compare the amino acids in list1 & list2
            if (i == j and i not in common):
                # add item to the list
                common.append(i)

    return common

def combine_proteins(protein_list):
    # create blank list for storing combine proteins
    combine = []

    # combine together both lists
    # go through each sublist in protein_list:
    for sub_list in protein_list:
        # go through each codon in the sublist
        for codon in sub_list:
            combine.append(codon)

    return combine

def amino_acid_dictionary(aa_list):
    # create blank dict for storing occurrences
    d = dict()
    # create blank list for keys
    keys = []

    # go through the protein list
    for amino_acid in aa_list:
        if amino_acid not in keys:
            keys.append(amino_acid)

    # get the occurrence for each key
    for key in keys:
        x = aa_list.count(key)
        d[key] = x

    return d

def find_amino_acid_differences(protein_list1, protein_list2, cutoff):
    # create combined proteins in each list
    combined_list1 = combine_proteins(protein_list1)
    combined_list2 = combine_proteins(protein_list2)

    # create dictionary of occurrences for each list
    occur_dict1 = amino_acid_dictionary(combined_list1)
    occur_dict2 = amino_acid_dictionary(combined_list2)

    # create dictionary of frequency for each list
    freq_dict1 = dict()
    freq_dict2 = dict()

    # create list of keys for each list
    keys1 = occur_dict1.keys()
    keys2 = occur_dict2.keys()

    # create list of values & sum for each list
    values1 = occur_dict1.values()
    sum1 = 0
    for value in values1:
        sum1 += value

    values2 = occur_dict2.values()
    sum2 = 0
    for value in values2:
        sum2 += value

    # populate frequency dictionaries
    for key in keys1:
        frequency = occur_dict1[key] / sum1
        freq_dict1[key] = frequency

    for key in keys2:
        frequency = occur_dict2[key] / sum2
        freq_dict2[key] = frequency

    # do not include "Start" and "Stop" amino acids
    freq_dict1.pop("Start")
    freq_dict1.pop("Stop")
    freq_dict2.pop("Start")
    freq_dict2.pop("Stop")

    # get set intersection of both lists
    set1 = set(keys1)
    set2 = set(keys2)
    set3 = set1.intersection(set2)

    # add amino acids that are common in both (intersection):
    sublist = []
    superlist = []
    for amino in set3:
        if amino in freq_dict1 and amino in freq_dict2:
            diff = abs(freq_dict1[amino] - freq_dict2[amino])
            if diff > cutoff:
                sublist.append(amino)
                sublist.append(freq_dict1[amino])
                sublist.append(freq_dict2[amino])
                superlist.append(sublist)
            sublist = []

    # add amino acids that don't occur in list2
    set1_diff = set1.difference(set2)
    sublist = []
    for amino in set1_diff:
        if freq_dict1[amino] > cutoff:
            sublist.append(amino)
            sublist.append(freq_dict1[amino])
            sublist.append(0)
            superlist.append(sublist)

    # add amino acids that don't occur in list1
    set2_diff = set2.difference(set1)
    sublist = []
    for amino in set2_diff:
        if freq_dict2[amino] > cutoff:
            sublist.append(amino)
            sublist.append(0)
            sublist.append(freq_dict2[amino])
            superlist.append(sublist)

    superlist.sort()
    return superlist

def display_commonalities(commonalities):
    # remove "Start" and "Stop" codons
    for proteins in commonalities:
        proteins.remove("Start")
        proteins.remove("Stop")

    # arrange proteins in alphabetical order
    commonalities.sort()

    # print out the results
    for proteins in commonalities:
        if len(proteins) > 1:
            separator = "-"
            x = separator.join(proteins)
            print(x)
        else:
            for protein in proteins:
                print(protein)

    return

def display_differences(differences):

    # print out results
    for proteins in differences:
        amino_acid = proteins[0]
        diff_1 = "{:.2%}".format(proteins[1])
        diff_2 = "{:.2%}".format(proteins[2])
        print(amino_acid, ":", diff_1, "in Seq1,", diff_2, "in Seq2")

    return

### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    test.test_all()
    test.run()