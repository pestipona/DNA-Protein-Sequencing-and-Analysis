import sys

sys.path.append("D:\\Portfolio Projects\\DNA Protein Sequencing and Analysis")
from protein import *


def test_read_file():
    print("Testing read_file()...", end="")
    # Check that the function works on the provided test file
    text1 = read_file("data/test_dna.txt")
    assert (text1 == "ATGGATGGACTCTAACGCAATGCCCTTTTAG")

    # Now check the human DNA file you loaded
    text2 = read_file("data/human_p53.txt")
    assert (text2[:10] == "GATGGGATTG")  # the whole sequence is too long to check here!
    assert (len(text2) == 19149)
    # If the length is not correct, check that you're
    # removing newlines, and that you copied the whole sequence
    print("... done!")


def test_dna_to_rna():
    print("Testing dna_to_rna()...", end="")
    # Test a basic sequence
    dna = "ATGGATGGACTCTAA"
    assert (dna_to_rna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    # Test two mRNA strands in a row, with a random codon in between
    dna = "ATGGATGGACTCTAACTCATGCCCTTTTAG"
    assert (dna_to_rna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    assert (dna_to_rna(dna, 18) == ["AUG", "CCC", "UUU", "UAG"])
    # Test a DNA strand that doesn't end properly
    dna = "CCTATGGACCAT"
    assert (dna_to_rna(dna, 3) == ["AUG", "GAC", "CAU"])
    # Test a DNA strand with random bases in between
    dna = "ATGGATGGACTCTAACGCAATGCCCTTTTAGAAA"
    assert (dna_to_rna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    assert (dna_to_rna(dna, 19) == ["AUG", "CCC", "UUU", "UAG"])
    print("... done!")


def test_make_codon_dictionary():
    print("Testing make_codon_dictionary()...", end="")
    d = make_codon_dictionary("data/codon_table.json")
    assert (d["AAA"] == "Lys")
    assert (d["GGA"] == "Gly")
    assert (d["AUG"] == "Met")
    assert (d["UAA"] == "Stop")
    print("... done!")


def test_generate_protein():
    print("Testing generate_protein()...", end="")
    codon_dict = make_codon_dictionary("data/codon_table.json")
    rna = ["AUG", "GAU", "GGA", "CUC", "UAA"]
    assert (generate_protein(rna, codon_dict) == ["Start", "Asp", "Gly", "Leu", "Stop"])
    rna = ["AUG", "CCC", "UUU", "UAG"]
    assert (generate_protein(rna, codon_dict) == ["Start", "Pro", "Phe", "Stop"])
    rna = ["AUG", "GAC", "CAU"]
    assert (generate_protein(rna, codon_dict) == ["Start", "Asp", "His"])
    # Note: "AUG" only maps to "Start" if it's at the beginning. "AUG" in the middle should become "Met"
    rna = ["AUG", "CGA", "AUG", "GGG", "UGG", "UGA"]
    assert (generate_protein(rna, codon_dict) == ["Start", "Arg", "Met", "Gly", "Trp", "Stop"])
    print("... done!")


def test_synthesize_proteins():
    print("Testing synthesize_proteins()...", end="")
    # First, test on the provided test data
    proteins1 = synthesize_proteins("data/test_dna.txt", "data/codon_table.json")
    # The function should say there are 31 total bases,
    # 4 unused bases, and 2 synthesized proteins
    assert (proteins1 == [['Start', 'Asp', 'Gly', 'Leu', 'Stop'],
                          ['Start', 'Pro', 'Phe', 'Stop']])

    # Now test on the actual data
    proteins2 = synthesize_proteins("data/human_p53.txt", "data/codon_table.json")
    # The function should say there are 19149 total bases,
    # 10560 unused bases, and 119 synthesized proteins
    assert (len(proteins2) == 119)
    assert (proteins2[0] == ['Start', 'Gly', 'Leu', 'Gly', 'Phe', 'Ser', 'Pro',
                             'Pro', 'Met', 'Cys', 'Ser', 'Arg', 'Leu', 'Ala',
                             'Leu', 'Lys', 'Val', 'Leu', 'Ser', 'Phe', 'Ser',
                             'Lys', 'Val', 'Stop'])
    assert (proteins2[1] == ['Start', 'Ser', 'Pro', 'Leu', 'Stop'])
    assert (proteins2[118] == ['Start', 'Met', 'Ile', 'Trp', 'Ile', 'His', 'Gln',
                               'Asp', 'Leu', 'Phe', 'Tyr', 'Ala', 'Gln', 'Gly',
                               'Gln', 'Phe', 'Leu', 'Phe', 'Ser', 'Phe', 'Phe',
                               'Phe', 'Phe', 'Phe', 'Phe', 'Phe', 'Phe', 'Phe',
                               'Glu', 'Thr', 'Gly', 'Ser', 'Arg', 'Phe', 'Val',
                               'Ala', 'Gln', 'Ala', 'Gly', 'Val', 'Glu', 'Trp',
                               'Arg', 'Asp', 'Leu', 'Gly', 'Leu', 'Leu', 'Gln',
                               'Pro', 'Leu', 'Pro', 'Pro', 'Arg', 'Leu', 'Glu',
                               'Gln', 'Ser', 'Cys', 'Leu', 'Ser', 'Leu', 'Arg',
                               'Ser', 'Ser', 'Trp', 'Asp', 'His', 'Arg', 'Phe',
                               'Met', 'Pro', 'Pro', 'Trp', 'Pro', 'Ala', 'Asn',
                               'Phe', 'Cys', 'Met', 'Phe', 'Cys', 'Arg', 'Asp',
                               'Gly', 'Val', 'Ser', 'Gln', 'Cys', 'Cys', 'Pro',
                               'Gly', 'Trp', 'Ser', 'Gln', 'Thr', 'Pro', 'Gly',
                               'Leu', 'Arg', 'Arg', 'Ser', 'Thr', 'Cys', 'Leu',
                               'Ser', 'Leu', 'Pro', 'Glu', 'Cys', 'Trp', 'Asp',
                               'Tyr', 'Asn', 'Cys', 'Glu', 'Pro', 'Pro', 'Arg',
                               'Pro', 'Ala', 'Gly', 'Arg', 'Val', 'Asn', 'Ile',
                               'Phe', 'Tyr', 'Ile', 'Leu', 'Gln', 'Ala', 'His',
                               'Leu', 'His', 'Phe', 'His', 'Pro', 'Thr', 'Leu',
                               'Pro', 'Leu', 'Leu', 'Leu', 'Pro', 'Phe', 'Tyr',
                               'Ile', 'Pro', 'Phe', 'Leu', 'Tyr', 'Arg', 'Ser',
                               'Leu', 'Ile', 'Leu', 'Gln', 'Stop'])
    print("... done!")


def test_common_proteins():
    print("Testing common_proteins()...", end="")
    plist1 = [["Start", "Pro", "Val", "Stop"], ["Start", "Phe", "Stop"],
              ["Start", "Asp", "Glu", "Stop"], ["Start", "His", "Stop"]]
    plist2 = [["Start", "Cys", "Cys", "Tyr", "Stop"], ["Start", "Glu", "Asp", "Stop"],
              ["Start", "His", "Stop"], ["Start", "Stop"], ["Start", "Met", "Leu", "Stop"]]
    plist3 = [["Start", "Asp", "Glu", "Stop"], ["Start", "Phe", "Stop"],
              ["Start", "Asp", "Glu", "Stop"], ["Start", "Lys", "Stop"],
              ["Start", "Asn", "Asn", "Asn", "Asn", "Stop"]]
    assert (common_proteins(plist1, plist2) == [["Start", "His", "Stop"]])
    assert (sorted(common_proteins(plist1, plist3)) == [["Start", "Asp", "Glu", "Stop"],
                                                        ["Start", "Phe", "Stop"]])
    assert (common_proteins(plist2, plist3) == [])
    print("... done!")


def test_combine_proteins():
    print("Testing combine_proteins()...", end="")
    plist1 = [["Start", "Pro", "Val", "Stop"], ["Start", "Phe", "Stop"],
              ["Start", "Asp", "Glu", "Stop"], ["Start", "His", "Stop"]]
    plist2 = [["Start", "Cys", "Cys", "Tyr", "Stop"], ["Start", "Glu", "Asp", "Stop"],
              ["Start", "His", "Stop"], ["Start", "Stop"], ["Start", "Met", "Leu", "Stop"]]
    plist3 = [["Start", "Asp", "Glu", "Stop"], ["Start", "Phe", "Stop"],
              ["Start", "Asp", "Glu", "Stop"], ["Start", "Lys", "Stop"],
              ["Start", "Asn", "Asn", "Asn", "Asn", "Stop"]]
    assert (combine_proteins(plist1) == ["Start", "Pro", "Val", "Stop", "Start",
                                         "Phe", "Stop", "Start", "Asp", "Glu",
                                         "Stop", "Start", "His", "Stop"])
    assert (combine_proteins(plist2) == ["Start", "Cys", "Cys", "Tyr", "Stop",
                                         "Start", "Glu", "Asp", "Stop", "Start",
                                         "His", "Stop", "Start", "Stop", "Start",
                                         "Met", "Leu", "Stop"])
    assert (combine_proteins(plist3) == ["Start", "Asp", "Glu", "Stop", "Start",
                                         "Phe", "Stop", "Start", "Asp", "Glu",
                                         "Stop", "Start", "Lys", "Stop", "Start",
                                         "Asn", "Asn", "Asn", "Asn", "Stop"])
    print("... done!")


def test_amino_acid_dictionary():
    print("Testing amino_acid_dictionary()...", end="")
    aa_list1 = ["Start", "Pro", "Val", "Stop", "Start", "Phe", "Stop", "Start",
                "Asp", "Glu", "Stop", "Start", "His", "Stop"]
    aa_list2 = ["Start", "Cys", "Cys", "Tyr", "Stop", "Start", "Glu", "Asp",
                "Stop", "Start", "His", "Stop", "Start", "Stop", "Start", "Met",
                "Leu", "Stop"]
    aa_list3 = ["Start", "Asp", "Glu", "Stop", "Start", "Phe", "Stop", "Start",
                "Asp", "Glu", "Stop", "Start", "Lys", "Stop", "Start", "Asn",
                "Asn", "Asn", "Asn", "Stop"]
    assert (amino_acid_dictionary(aa_list1) == {"Start": 4, "Pro": 1, "Val": 1,
                                                "Stop": 4, "Phe": 1, "Asp": 1, "Glu": 1, "His": 1})
    assert (amino_acid_dictionary(aa_list2) == {"Start": 5, "Cys": 2, "Tyr": 1,
                                                "Stop": 5, "Glu": 1, "Asp": 1, "His": 1, "Met": 1, "Leu": 1})
    assert (amino_acid_dictionary(aa_list3) == {"Start": 5, "Asp": 2, "Glu": 2,
                                                "Stop": 5, "Phe": 1, "Lys": 1, "Asn": 4})
    print("... done!")


def test_find_amino_acid_differences():
    print("Testing find_amino_acid_differences()...", end="")
    set1 = [['Start', 'Gly', 'Leu', 'Gly', 'Phe', 'Ser', 'Pro', 'Pro', 'Met',
             'Cys', 'Ser', 'Arg', 'Leu', 'Ala', 'Leu', 'Lys', 'Val', 'Leu',
             'Ser', 'Phe', 'Ser', 'Lys', 'Val', 'Stop'],
            ['Start', 'Ser', 'Pro', 'Leu', 'Stop'],
            ['Start', 'Glu', 'Ala', 'Trp', 'Leu', 'Glu', 'Gly', 'Ser', 'Ser', 'Stop'],
            ['Start', 'Met', 'Gly', 'Met', 'Leu', 'Gly', 'Pro', 'Ser', 'Glu',
             'Leu', 'Lys', 'Val', 'Glu', 'Arg', 'Leu', 'Gly', 'Arg', 'Gly',
             'Val', 'Glu', 'Leu', 'Trp', 'Gly', 'Thr', 'Leu', 'Ser', 'Arg',
             'Pro', 'Lys', 'Ala', 'Tyr', 'Phe', 'Phe', 'Ala', 'His', 'Pro',
             'Pro', 'Gly', 'Ala', 'Gly', 'Arg', 'Arg', 'Glu', 'Ser', 'Leu',
             'Lys', 'Stop'],
            ['Start', 'His', 'Lys', 'Ala', 'Leu', 'Arg', 'Ser', 'Glu', 'Thr',
             'Phe', 'Gly', 'Ser', 'Arg', 'Asn', 'Ile', 'Glu', 'Asn', 'Ser', 'Stop']]
    set2 = [['Start', 'Ala', 'Stop'],
            ['Start', 'Phe', 'Ser', 'Ile', 'Asn', 'Ser', 'Thr', 'Leu', 'Ala',
             'Ala', 'Leu', 'Val', 'Cys', 'Arg', 'Thr', 'Ser', 'Pro', 'Pro',
             'Gln', 'Asn', 'Pro', 'Gly', 'Ser', 'Leu', 'Arg', 'Ser', 'Leu',
             'Leu', 'Phe', 'His', 'Ser', 'Leu', 'Ser', 'Ala', 'Ser', 'Pro',
             'Leu', 'Pro', 'Thr', 'Gly', 'Lys', 'Leu', 'Leu', 'Ala', 'Leu',
             'Thr', 'Cys', 'His', 'Gly', 'Asp', 'Cys', 'Pro', 'Ala', 'Leu',
             'Cys', 'Gln', 'Lys', 'Pro', 'Arg', 'Gly', 'Gly', 'Cys', 'Trp',
             'Asp', 'Trp', 'Glu', 'Phe', 'Pro', 'Phe', 'Pro', 'Cys', 'Ala',
             'His', 'Thr', 'Gly', 'Ala', 'Lys', 'Ser', 'Phe', 'Gln', 'Leu',
             'Phe', 'Lys', 'Ser', 'Pro', 'Lys', 'Pro', 'Pro', 'Ser', 'Trp',
             'Leu', 'Gln', 'Leu', 'Ala', 'Ala', 'Gly', 'Leu', 'Trp', 'Arg',
             'Tyr', 'Leu', 'Val', 'Ser', 'Gly', 'Leu', 'Gly', 'Pro', 'Cys',
             'Phe', 'Gln', 'Gly', 'Arg', 'Leu', 'His', 'Ala', 'Arg', 'Leu',
             'Arg', 'Phe', 'Gly', 'Stop'],
            ['Start', 'Ser', 'Pro', 'Leu', 'Stop'],
            ['Start', 'Phe', 'Arg', 'Ala', 'Leu', 'Gly', 'Val', 'Glu', 'Stop'],
            ['Start', 'Leu', 'Val', 'Pro', 'Ala', 'Asp', 'Leu', 'Glu', 'Leu', 'Stop']]
    result1 = find_amino_acid_differences(set1, set2, 0.02)  # 2% difference
    result1.sort()
    assert (len(result1) == 12)
    assert ((result1[0][0] == "Ala") and (0.057 < result1[0][1] < 0.058) and (0.087 < result1[0][2] < 0.088))
    assert ((result1[1][0] == "Arg") and (0.076 < result1[1][1] < 0.077) and (0.054 < result1[1][2] < 0.055))
    assert ((result1[11][0] == "Ser") and (0.123 < result1[11][1] < 0.124) and (0.087 < result1[11][2] < 0.088))

    result2 = find_amino_acid_differences(set1, set2, 0.05)  # 5% difference
    assert (len(result2) == 1)
    assert ((result2[0][0] == "Glu") and (0.076 < result2[0][1] < 0.077) and (0.020 < result2[0][2] < 0.021))

    result3 = find_amino_acid_differences(set1, set2, 0.005)  # 0.5% difference
    assert (len(result3) == 18)
    print("... done!")


def test_display_commonalities():
    human_proteins = synthesize_proteins("data/human_p53.txt", "data/codon_table.json")
    elephant_proteins = synthesize_proteins("data/elephant_p53.txt", "data/codon_table.json")
    commonalities = common_proteins(human_proteins, elephant_proteins)
    display_commonalities(commonalities)
    assert (commonalities)


def test_all():
    test_read_file()
    test_dna_to_rna()
    test_make_codon_dictionary()
    test_generate_protein()
    test_synthesize_proteins()
    test_common_proteins()
    test_combine_proteins()
    test_amino_acid_dictionary()
    test_find_amino_acid_differences()


def run():
    print("\n-----\n")
    print("Now let's try running everything!\n")

    human_proteins = synthesize_proteins("data/human_p53.txt", "data/codon_table.json")
    elephant_proteins = synthesize_proteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = common_proteins(human_proteins, elephant_proteins)
    print("The following proteins occurred in both DNA Sequences:")
    display_commonalities(commonalities)

    differences = find_amino_acid_differences(human_proteins, elephant_proteins, 0.005)
    print("\nThe following amino acids occurred at very different rates in the two DNA sequences:")
    display_differences(differences)