import tile_analysis
from collections import defaultdict
import os
import pandas as pd

barcode_counts = defaultdict(int)
home = os.path.split(os.getcwd())[:-1]
os.chdir(os.path.join(*home, 'projects', "example"))

if os.path.exists("process_ipynb/cells.csv"):
    os.remove("process_ipynb/cells.csv")

for tile in range(1):
    # change these values based on your data!
    tile_analysis.run(
        well=1,
        tile=107,
        cycles=12,
        data_path=f'input/*/*.tif',  # change so it will match
        # the paths of your images. This example would match data/10x_Cycle1_Well3_Point3_0000_ChannelDAPI,G-ISS,
        # T-ISS,A-ISS,C-ISS_Seq0784.ome for example
        project_name="example",  # name of subdirectory under `projects` where you put your data
        barcode_csv_name="barcodes.csv",  # name of csv file with barcode library.
        threshold_dapi=2000,  # adjust if nuclei.tif doesn't detect all nuclei or shows background spots
        threshold_cell=1500,  # adjust if cells are too big or too small
        threshold_reads=500,  # adjust if percentage of reads that are in library is unexpected
        nucleus_area=(40, 400),
        DAPI_index=0,
        align_method="SBS_mean",  # should be "SBS_mean" or "DAPI"
        barcode_counts=barcode_counts,
    )

    print(sum(barcode_counts.values()) / sum(pd.read_csv("process_ipynb/cells.csv")["barcode_count"]) * 100)  # print
    # percentage of reads are found in barcode library

print(sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True))  # print all barcodes which are found in
# the barcode library, in order from most to least common

# 500 - 7.54
# 600 - 7.82
# 700 - 7.77
