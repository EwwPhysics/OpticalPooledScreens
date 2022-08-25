import tile_analysis
from collections import defaultdict
import os
import pandas as pd

barcode_counts = defaultdict(int)
home = os.path.split(os.getcwd())[:-1]
os.chdir(os.path.join(*home, 'projects', "steph"))

for tile in range(1):
    # change these values based on your data!
    tile_analysis.run(
        well=3,
        tile=tile,
        cycles=11,
        data_path=f'data/10x_Cycle*_Well{3}_Point3_{str(tile).rjust(4, "0")}*.ome.tif',  # change so it will match
        # the paths of your images
        project_name="steph",  # name of subdirectory under `projects` where you put your data
        threshold_dapi=2000,  # adjust if nuclei.tif doesn't detect all nuclei or shows background spots
        threshold_cell=2500,  # adjust if cells are too big or too small
        threshold_reads=600,  # adjust if percentage of reads are in library is unexpected
        nucleus_area=(40, 400),
        DAPI_index=0,
        barcode_counts=barcode_counts,
    )
    print(sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True))
    print(sum(barcode_counts.values()) / (pd.read_csv("process_ipynb/cells.csv").shape[0] * 2) * 100)

os.remove("process_ipynb/cells.csv")

# 500 - 7.54
# 600 - 7.82
# 700 - 7.77
# [('ATAATAGCATT', 26), ('GCCCCGCCGCC', 19), ('CTCAACCACTT', 17), ('GTTGCCAGCTT', 13), ('GAAGAAGGATG', 12), ('AGGTCTGCCCA', 12), ('GACAGACTGGT', 11), ('GATCGCTTTCC', 11), ('TCTAGCTGAGT', 11), ('AAAGCTATGTT', 9), ('CCCACCACCCC', 9), ('CTCTCAGACCT', 8), ('GCCTGCCCTAA', 8), ('GCACTATTAAT', 7), ('ACTACTAGCTC', 7), ('CGAAAGATACC', 6), ('CTCAAACTTAT', 6), ('ATGAGTCAAAG', 6), ('TTGCCTCACCC', 6), ('CTCCCCGACTA', 5), ('TGGACAGTGCC', 5), ('CTGATCAGGAT', 5), ('ACCATCGCATG', 4), ('GCGATGGCCTC', 4), ('GTTGATCGAAA', 4), ('CCCACACCCCC', 4), ('GGAATGGGAAG', 4), ('TCTTTCCAACT', 4), ('TTAGCCGGATG', 4), ('ACTCTCTAATA', 3), ('GTATGGGCTGC', 3), ('CCAACAACGTA', 3), ('CCCTCCCCGCC', 3), ('CCAGCTTTGAT', 3), ('CGCCTGCGCAC', 3), ('AGCCCCTCCTT', 3), ('GTATATCCTCT', 3), ('ACACAGTGTTT', 3), ('GGTGAAGTCGT', 3), ('GGGCATGAGAT', 3), ('AGGACCCCGCC', 2), ('ACTACAGTCCA', 2), ('GGTCGTGGGTC', 2), ('TGTCGACAAAT', 2), ('CACGCCCCCCA', 2), ('ACTACTAGCAC', 2), ('AGGTGGGGCTG', 2), ('AAGTTCTTGAA', 2), ('AGGCCTCGTGC', 2), ('GGTTAACTTAA', 2), ('GAAAGCACGCT', 2), ('GTGCCACAAAC', 2), ('GGGAGGCTGAT', 2), ('GGTCCCGCAGC', 2), ('GGGCTTCTGCT', 2), ('AGAACCGGTCA', 2), ('ATTCGCCCATC', 2), ('GAAATCGCCCC', 2), ('TTCTCTATGTA', 2), ('ACATACGTGAT', 2), ('AGATGCGAGCT', 2), ('GGGCCAAGAAT', 2), ('ACACACACGCT', 2), ('GATCTCAACTC', 2), ('CTGGCCCCCAT', 2), ('GTGAGAGCGTA', 2), ('CCTCTCCCAAT', 2), ('AACCGAGAGGG', 1), ('AGAAAGCCCGC', 1), ('CACCGCACACC', 1), ('CCCCCACGAGA', 1), ('CCCCCAAAGCG', 1), ('ACCCGACAAAT', 1), ('CCCATCTCTAT', 1), ('CCACAACAGAA', 1), ('CAAGGCTATGG', 1), ('AAGAGCTGGCC', 1), ('TAACTACTGTC', 1), ('ACCCCGCGGGC', 1), ('ACGCAGGTTAA', 1), ('GAACTGTGATG', 1), ('TGAGCTTCTGG', 1), ('GGTCCCCAACC', 1), ('AATAAGGCATT', 1), ('GAGGCCCACCC', 1), ('CATAGGAGGAG', 1), ('CCGCCACGCCC', 1), ('ATGGCCCACCC', 1), ('GTTTCGAAACT', 1), ('GTTTTACATCT', 1), ('GTGGCAGGGGA', 1), ('TGCACAGAGAT', 1), ('GACAGAACCAT', 1), ('CTACAATTCCA', 1), ('TCAAGATCACA', 1), ('CCACTATTCCC', 1), ('AACATGGCACT', 1), ('GATGGCTAACC', 1), ('AAAGAAAAAAC', 1), ('CCCACCACACC', 1), ('CCACCAAGTCA', 1), ('ACCCAGATTAG', 1), ('CCCGCCCCATC', 1), ('AAACGCCTTTC', 1), ('CCCACTGTCCT', 1), ('AAATGCAGAGC', 1), ('TGCTAATGCTG', 1), ('GACACACTGAC', 1), ('ATTGTTACACA', 1), ('GAACACGAGAT', 1), ('AACTAACTACT', 1), ('GAGAAGATAGC', 1), ('GCCCCCCCAGG', 1), ('GAGAGGCAGCC', 1), ('ACACCCCGCCC', 1), ('ATGCGCAGCTC', 1), ('GCTGGTTCATC', 1), ('GGCAGTACTTG', 1), ('CCCAGGAGACT', 1), ('ACAGCTTGAGT', 1), ('TTGTGCTCTCC', 1), ('ATGCTGATCTA', 1), ('TTGTGCTGTCC', 1), ('CACCCCCAGGC', 1), ('AACTTGCAACC', 1), ('GAGCAACATTT', 1), ('CTTGGCAGGAG', 1), ('GAATCGAGAGG', 1), ('CATGGCCTACG', 1), ('GTGATAATGTT', 1), ('CCTGGGTCTAC', 1), ('TTCACACATCT', 1), ('GGGCGGGGGAA', 1), ('TCCACCACCAC', 1), ('GTCAAGGAATT', 1), ('TTGGACTGGGC', 1), ('CGGAGAGGGGA', 1), ('GAGTCCATCAA', 1), ('CGGGGACTCAG', 1), ('TTCTCAAACTC', 1), ('ACCCAGGGGCG', 1), ('CAGCAGGAGCG', 1), ('TTCACCGTCCA', 1), ('AAACCTCGCTG', 1), ('GCCCCGCGGAC', 1), ('AGATCAGGTGG', 1), ('TGCTCAGAGAT', 1), ('CAGGGGGCATT', 1), ('CAGCCGCCGCC', 1), ('GCCCCAAGAAT', 1), ('CCCCAAAGGAA', 1), ('ACAGTTAAATT', 1), ('CCCACCAAAGA', 1), ('GCTTCAAGTAT', 1), ('TATGCTAATGC', 1), ('GCGGCCCACAC', 1), ('GAGGCCCCCTC', 1), ('GAACTCGGACA', 1), ('TGTTATCCTGT', 1), ('ATGAGCAGTGC', 1), ('CCCACCAACGG', 1), ('ATGTGGAAAAA', 1), ('TCGCTCGCCAC', 1), ('TGCCGCTTGCC', 1), ('ACTCCAAATAA', 1), ('CCGCCCGCCAG', 1), ('GAGAGGCGGGG', 1), ('CAACCAGTACC', 1), ('GTGCAATATGA', 1), ('GATTTGGTCTC', 1), ('ATGAGCACCAC', 1)]
# [('ATAATAGCATT', 60), ('CTCAACCACTT', 28), ('GAAGAAGGATG', 22), ('CCCACCACCCC', 20), ('GTTGCCAGCTT', 19), ('GCCCCGCCGCC', 19), ('GGGCCAAGAAT', 18), ('TCTAGCTGAGT', 18), ('AAAGCTATGTT', 14), ('TTAGCCGGATG', 14), ('GATCGCTTTCC', 13), ('AACATGGCACT', 13), ('AGGTCTGCCCA', 13), ('GACAGACTGGT', 12), ('CCCACACCCCC', 12), ('TCTTTCCAACT', 11), ('GTGATTCCTTC', 11), ('GGGCATGAGAT', 10), ('CGAAAGATACC', 9), ('AGCTTACAATA', 9), ('ACTACAGTCCA', 8), ('CTCCCCGACTA', 8), ('CTCTCAGACCT', 8), ('ACTACTAGCTC', 8), ('CTCAAACTTAT', 8), ('GCCTGCCCTAA', 8), ('ATGAGTCAAAG', 8), ('GCTAGTACCCA', 8), ('GCACTATTAAT', 7), ('TGGACAGTGCC', 7), ('GGAATGGGAAG', 7), ('TTGCCTCACCC', 7), ('CCAGCTTTGAT', 7), ('CGCCTGCGCAC', 7), ('GCTTCAAGTAT', 7), ('TCACAGTTTTC', 7), ('GTCCAATTTTC', 7), ('GTTGATCGAAA', 6), ('AGAGAATTCTG', 6), ('GACCCCTAGAT', 6), ('ACCATTAGGCC', 6), ('GGGCTTCTGCT', 5), ('CTGATCAGGAT', 5), ('AGCCCCTCCTT', 5), ('AGCCTGACCTT', 5), ('GGTCGTGGGTC', 4), ('TGTCGACAAAT', 4), ('ACCATCGCATG', 4), ('GCGATGGCCTC', 4), ('CCCTCCCCGCC', 4), ('GTATATCCTCT', 4), ('ACACAGTGTTT', 4), ('TGATCTCCCAT', 4), ('CAGCAGAAGTG', 4), ('CAATGAAAGTG', 4), ('CACGCCCCCCA', 3), ('AGGCCTCGTGC', 3), ('ACTCTCTAATA', 3), ('GTATGGGCTGC', 3), ('CCAACAACGTA', 3), ('GTGCCACAAAC', 3), ('TTCTCTATGTA', 3), ('AGATGCGAGCT', 3), ('GATCTCAACTC', 3), ('TTCACCGTCCA', 3), ('GGTGAAGTCGT', 3), ('CAGCCGCCGCC', 3), ('CCTCTCCCAAT', 3), ('CAGAGTAATAT', 3), ('AACCGAGCTCT', 3), ('CTAATCCCTGA', 3), ('TAGTCTACATG', 3), ('CCAATCAGTGC', 3), ('GAACCACCACC', 3), ('AGGACCCCGCC', 2), ('ACTACTAGCAC', 2), ('AGGTGGGGCTG', 2), ('AAGTTCTTGAA', 2), ('CCGCCACGCCC', 2), ('GGTTAACTTAA', 2), ('GAAAGCACGCT', 2), ('GGGAGGCTGAT', 2), ('GGTCCCGCAGC', 2), ('ACACCCCGCCC', 2), ('ATGCGCAGCTC', 2), ('AGAACCGGTCA', 2), ('ATTCGCCCATC', 2), ('GAAATCGCCCC', 2), ('AACTTGCAACC', 2), ('ACATACGTGAT', 2), ('ACACACACGCT', 2), ('TGCTCAGAGAT', 2), ('CTGGCCCCCAT', 2), ('GTGAGAGCGTA', 2), ('ATGAGCACCAC', 2), ('TCGCAAGGAAG', 2), ('AACAACCACCC', 2), ('CCAGCCCACCC', 2), ('CCGACCAAAAA', 2), ('CACTATTAATG', 2), ('GTTCACCAACA', 2), ('TCTAATAGTAT', 2), ('GAGGAGCTTGA', 2), ('CAAGGACCCTA', 2), ('CCCGTTCATCC', 2), ('GAGAAAAAGAA', 2), ('CGTCTAGGAGG', 2), ('TTGACATTGCC', 2), ('GCTTCGAGGCT', 2), ('CGTCTGTTCAG', 2), ('AACCGAGAGGG', 1), ('AGAAAGCCCGC', 1), ('CACCGCACACC', 1), ('CCCCCACGAGA', 1), ('CCCCCAAAGCG', 1), ('ACCCGACAAAT', 1), ('CCCATCTCTAT', 1), ('CCACAACAGAA', 1), ('CAAGGCTATGG', 1), ('AAGAGCTGGCC', 1), ('TAACTACTGTC', 1), ('ACCCCGCGGGC', 1), ('ACGCAGGTTAA', 1), ('GAACTGTGATG', 1), ('TGAGCTTCTGG', 1), ('GGTCCCCAACC', 1), ('AATAAGGCATT', 1), ('GAGGCCCACCC', 1), ('CATAGGAGGAG', 1), ('ATGGCCCACCC', 1), ('GTTTCGAAACT', 1), ('GTTTTACATCT', 1), ('GTGGCAGGGGA', 1), ('TGCACAGAGAT', 1), ('GACAGAACCAT', 1), ('CTACAATTCCA', 1), ('TCAAGATCACA', 1), ('CCACTATTCCC', 1), ('GATGGCTAACC', 1), ('AAAGAAAAAAC', 1), ('CCCACCACACC', 1), ('CCACCAAGTCA', 1), ('ACCCAGATTAG', 1), ('CCCGCCCCATC', 1), ('AAACGCCTTTC', 1), ('CCCACTGTCCT', 1), ('AAATGCAGAGC', 1), ('TGCTAATGCTG', 1), ('GACACACTGAC', 1), ('ATTGTTACACA', 1), ('GAACACGAGAT', 1), ('AACTAACTACT', 1), ('GAGAAGATAGC', 1), ('GCCCCCCCAGG', 1), ('GAGAGGCAGCC', 1), ('GCTGGTTCATC', 1), ('GGCAGTACTTG', 1), ('CCCAGGAGACT', 1), ('ACAGCTTGAGT', 1), ('TTGTGCTCTCC', 1), ('ATGCTGATCTA', 1), ('TTGTGCTGTCC', 1), ('CACCCCCAGGC', 1), ('GAGCAACATTT', 1), ('CTTGGCAGGAG', 1), ('GAATCGAGAGG', 1), ('CATGGCCTACG', 1), ('GTGATAATGTT', 1), ('CCTGGGTCTAC', 1), ('TTCACACATCT', 1), ('GGGCGGGGGAA', 1), ('TCCACCACCAC', 1), ('GTCAAGGAATT', 1), ('TTGGACTGGGC', 1), ('CGGAGAGGGGA', 1), ('GAGTCCATCAA', 1), ('CGGGGACTCAG', 1), ('TTCTCAAACTC', 1), ('ACCCAGGGGCG', 1), ('CAGCAGGAGCG', 1), ('AAACCTCGCTG', 1), ('GCCCCGCGGAC', 1), ('AGATCAGGTGG', 1), ('CAGGGGGCATT', 1), ('GCCCCAAGAAT', 1), ('CCCCAAAGGAA', 1), ('ACAGTTAAATT', 1), ('CCCACCAAAGA', 1), ('TATGCTAATGC', 1), ('GCGGCCCACAC', 1), ('GAGGCCCCCTC', 1), ('GAACTCGGACA', 1), ('TGTTATCCTGT', 1), ('ATGAGCAGTGC', 1), ('CCCACCAACGG', 1), ('ATGTGGAAAAA', 1), ('TCGCTCGCCAC', 1), ('TGCCGCTTGCC', 1), ('ACTCCAAATAA', 1), ('CCGCCCGCCAG', 1), ('GAGAGGCGGGG', 1), ('CAACCAGTACC', 1), ('GTGCAATATGA', 1), ('GATTTGGTCTC', 1), ('GCATCACTAAA', 1), ('AGGACCCGTAT', 1), ('AGAAAAAATGT', 1), ('CAACATGCCGA', 1), ('GCAGGGGCGGG', 1), ('GCAACAGAGCC', 1), ('GTTGGGTGGCA', 1), ('GCACAAAACCA', 1), ('AAACTTCATCT', 1), ('CCGAGCTGAGT', 1), ('AGTTCTAGAAT', 1), ('CTGAGCAGGAC', 1), ('GGGCCGTTAAA', 1), ('ACACCCTTGGT', 1), ('CAGCACGTGAG', 1), ('ATAACAGTGGG', 1), ('TTGGATCATGG', 1), ('AGTTTGTATAG', 1), ('TGGTTAAGTGA', 1), ('AGGATGCTGAA', 1), ('CCCCAACAAAA', 1), ('GGGGGGAATCA', 1), ('CATCCATGAAC', 1), ('TACGACAAGCT', 1), ('CACGAACTCAC', 1), ('AATATGTCGTA', 1), ('AGGCAAAAAAA', 1), ('CCCCCGCCGGG', 1), ('GAGGCACAGAC', 1), ('AGATCCTTGCC', 1), ('CAGTGCGTTGA', 1), ('GGGACAAGAAG', 1), ('TTGAGGGGGAG', 1), ('CCCAAACCCCC', 1), ('GGCATCCAATG', 1), ('GGGCCCAGCAG', 1), ('GGGTGAATAAT', 1), ('TTGGGCCAACC', 1), ('CCCCAACAGGC', 1), ('CCTAGTATCCA', 1), ('CTGGCTCGGGA', 1), ('GTGCGGGGTCG', 1), ('ATAAAAATAAC', 1), ('TTGCTAGGGAA', 1), ('GGCAGGGAGAG', 1), ('GGTACTTGTGC', 1), ('TCACTTAAAAA', 1), ('TCAGTTAAAAA', 1), ('CCCGATCCCAA', 1), ('CCCAAGACGCC', 1), ('CTCTCATAGAC', 1), ('GGGAAAGGGAC', 1), ('CCCACAGGAGG', 1), ('CCAGGGTAGGC', 1), ('GTACTTCATGA', 1), ('ACCTTACAATA', 1), ('AGCAAAAACTA', 1), ('CACACAACGCC', 1), ('CATAACCCTTG', 1), ('AACAGGAGAGA', 1), ('TCTTCTACATG', 1), ('AAGTCACAAAG', 1), ('TGGTCCTCGAT', 1), ('AAGTAAAAACT', 1), ('CGCCACCCCGC', 1), ('TGCTTTTCATG', 1), ('CAGCCACTGAC', 1), ('ATCCCCGCCAG', 1), ('ACACAACAGAC', 1), ('GGGAAACTTAA', 1), ('TCCCCCCACCC', 1), ('CACTCTTCGGT', 1), ('CCCCCACAAAC', 1), ('CCCCAAAAGAG', 1), ('CCCAGTCTCAA', 1), ('GCAGTAGGGGT', 1), ('TCTATTTGAAA', 1), ('AGGAGGGGGAG', 1), ('GCAAAAACCCA', 1), ('CCAACCCCGAC', 1)]