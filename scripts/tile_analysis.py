from collections import defaultdict
from glob import glob
import os

from natsort import natsorted
import numpy as np
import pandas as pd

from ops import io, annotate, filenames
from ops.firesnake import Snake


def run(well: int,
        tile: int,
        cycles: int,
        project_name: str = "steph",
        threshold_dapi: int = 2000,
        threshold_cell: int = 2500,
        threshold_reads: int = 50,
        nucleus_area: tuple = (40, 400),
        DAPI_index: int = 0,
        barcode_counts: defaultdict = None,
        ):
    if project_name not in os.getcwd():
        home = os.path.split(os.getcwd())[:-1]
        os.chdir(os.path.join(*home, 'projects', project_name))
        print(os.getcwd())

    wildcards = {"well": well, "tile": tile}
    if barcode_counts is None:
        barcode_counts = defaultdict(int)

    # color of bases
    # lut = "lookup table", used to map one color to another like a filter
    luts = [
        io.GRAY,
        io.GREEN,
        io.RED,
        io.MAGENTA,
        io.CYAN
    ]

    # for formatting tif images when they are saved?
    DISPLAY_RANGES = [
        [500, 15000],
        [100, 10000],
        [100, 20000],
        [100, 8000],
        [100, 6000]
    ]

    barcodes = pd.read_csv('design.csv')  # list of barcodes along with which gene they target
    barcodes["prefix"] = barcodes["barcode"].apply(lambda x: x[:cycles])
    barcode_set = set(barcodes["prefix"])

    # find sbs images and print paths
    search = f'data/10x_Cycle*_Well{wildcards["well"]}_Point3_{str(wildcards["tile"]).rjust(4, "0")}*.ome.tif'
    input_files = natsorted(glob(search))

    # used to format output filenames
    description = {'mag': "10X", "well": wildcards["well"], 'tile': wildcards['tile'],
                   'subdir': f'process_ipynb/tile{wildcards["tile"]}', 'ext': 'tif'}

    data = np.array([io.read_stack(f) for f in input_files])

    data = Snake.align_SBS(data, method="SBS_mean")  # rigid alignment of sequencing cycles
    # save(name(description, tag='aligned'), data, display_ranges=DISPLAY_RANGES, luts=LUTS)

    logged = Snake.transform_log(data, skip_index=0)  # apply Laplacian-of-Gaussian filter from scipy.ndimage.
    # save(name(description, tag='log'), logged, display_ranges=DISPLAY_RANGES, luts=LUTS)

    maxed = Snake.max_filter(logged, 3, remove_index=0)  # apply a maximum filter in a window of `width`.
    # Conventionally operates on Laplacian-of-Gaussian filtered SBS data, dilating sequencing channels to compensate
    # for single-pixel alignment error.
    # save(name(description, tag='maxed'), maxed, display_ranges=DISPLAY_RANGES[1:], luts=LUTS[1:])

    std = Snake.compute_std(logged, remove_index=0)  # use standard deviation over cycles, followed by mean across
    # channels to estimate sequencing read locations.
    # save(name(description, tag='std'), std)

    peaks = Snake.find_peaks(std)  # where are the spots
    # save(name(description, tag='peaks'), peaks)

    # Find nuclei from DAPI (fluorescent stain)
    # change first argument if DAPI staining is only done for a certain cycle, e.g. 11th cycle would be data[10]
    nuclei = Snake.segment_nuclei(data[DAPI_index], threshold_dapi,
                                  area_min=nucleus_area[0], area_max=nucleus_area[1])

    io.save_stack(filenames.name_file(description, tag='nuclei'), nuclei, compress=1)

    cells = Snake.segment_cells(data[DAPI_index], nuclei, threshold_cell)  # Matches cell labels to nuclei labels.
    io.save_stack(filenames.name_file(description, tag='cells'), cells, compress=1)

    # Find the signal intensity from `maxed` at each point in `peaks` above `threshold_peaks`.
    df_bases = Snake.extract_bases(maxed, peaks, cells,
                                   threshold_reads, wildcards=wildcards)
    # df_bases.to_csv(name(description, tag='bases', ext='csv'), index=None)

    df_reads = Snake.call_reads(df_bases, peaks=peaks)  # call reads by compensating for channel cross-talk and
    # calling the base with the highest corrected intensity for each cycle. Q = quality?
    filename = filenames.name_file(description, tag='reads', ext='csv')
    df_reads.to_csv(filename, index=None)

    df_cells = Snake.call_cells(df_reads, df_pool=barcodes)  # gets the two most-common barcode reads for each cell
    # prioritizes barcodes found in gene pool
    df_cells.to_csv(filenames.name_file(description, tag='cells', ext='csv'), index=None)

    if not os.path.exists("process_ipynb/cells.csv") or os.stat("process_ipynb/cells.csv").st_size == 0:
        with open("process_ipynb/cells.csv", "w") as f:
            f.write("cell,tile,well,peak,barcode_0,count_0,barcode_1,count_1,barcode_count,sgRNA_0,gene_symbol_0,"
                    "sgRNA_1,gene_symbol_1\n")
    if wildcards["tile"] not in set(pd.read_csv("process_ipynb/cells.csv")["tile"]):
        df_cells.to_csv("process_ipynb/cells.csv", mode="a", index=False, header=False)

    # last channel annotates base calls
    annotate_luts = luts + [annotate.GRMC, io.GRAY]
    annotate_display_ranges = [(a / 4, b / 4) for a, b in DISPLAY_RANGES] + [(0, 4)]
    annotate_SBS = Snake.annotate_SBS(log=logged, df_reads=df_reads)

    io.save_stack(filenames.name_file(description, tag='annotate_SBS'), annotate_SBS,
                  display_ranges=annotate_display_ranges, luts=annotate_luts, compress=1)

    in_library = 0

    for index, row in df_cells.iterrows():
        b0 = row["barcode_0"]
        b1 = row["barcode_1"]
        if b0 in barcode_set:
            barcode_counts[b0] += int(row["count_0"])
            in_library += int(row["count_0"])
        if b1 in barcode_set:
            barcode_counts[b1] += int(float(row["count_1"]))
            in_library += int(float(row["count_1"]))

    # print(in_library / (df_cells.shape[0] * 2) * 100)
