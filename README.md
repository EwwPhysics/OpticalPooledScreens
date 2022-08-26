## Optical Pooled Screens

Analysis resources for the publication [*Pooled genetic perturbation screens with image-based phenotypes*](https://pubmed.ncbi.nlm.nih.gov/35022620/). This is an updated and expanded version of the repository for the original 2019 publication [*Optical pooled screens in human cells*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6886477/), available [here](https://github.com/feldman4/OpticalPooledScreens_2019).

### Installation (OSX)

Download the repository (e.g., on Github use the green "Clone or download" button, then "Download ZIP").

In Terminal, go to the OpticalPooledScreens project directory and create a Python 3 virtual environment using a command like:

```bash
python3 -m venv venv
```

If the python3 command isn't available, you might need to specify the full path. E.g., if [Miniconda](https://conda.io/miniconda.html) is installed in the home directory:

```bash
~/miniconda3/bin/python -m venv venv
```

This creates a virtual environment called `venv` for project-specific resources. The commands in `install.sh` add required packages to the virtual environment:

```bash
sh install.sh
```

Within this script, the `ops` package is installed with `pip install -e`, so the source code in the `ops/` directory can be modified in place.

Once install is complete, activate the virtual environment from the project directory:

```bash
source venv/bin/activate
```

Additionally, if using the CellPose segmentation method, this must be installed in the virtual environment:
```bash
pip install cellpose[gui]
```

## Running the pipeline on example data

Simply run the `well_analysis` file! The provided example data only contains one tile.

## Using your own data!

1. In the `projects/` directory, create a new subdirectory for your experiment. You should add a csv file with your barcodes to this folder (specifically, with columns `barcode`, `sgRNA`, and `gene_symbol` - see `barcodes.csv` in `projects/example` for an example). Within this new subdirectory, you should add another subdirectory called `data`. Then add your images to this folder. 

2. In the `scripts` directory, open `well_analysis.py`. Here, there are several factors that you will have to change based on your data. 

3. The final argument passed to `os.path.join()` should be changed to whatever your experiment’s directory is called. (`os.chdir(os.path.join(*home, 'projects', "example"))`)
4. The file contains a loop, which iterates over a range to run the analysis on each tile number. This should be changed accordingly based on your data. For example, if you had 100 tiles, the loop declaration should be: `for tile in range(100)`.Note that if your numbering starts from 1, not 0, you will have to use a slightly different range: `for tile in range(1, 101)`
5. The main loop in  well_analysis.py calls a function from tile_analysis.py with many keyword arguments. Here, I will describe what each argument means
   1. `well` - should be the well number from your experiment. This is important for naming files and reading your images.
   2. `tile` - you shouldn’t have to change this, if you followed the instructions in step 4.
   3. `cycles` - the number of imaging cycles in your experiment.
   4. `data_path` - this should reflect what the path to your images looks like, including the file names of the images. The goal is to match the images of the current tile from each cycle. It should start with the name of the subdirectory in your project folder that has all of your data. From here, it can vary.
      1. \* matches any number of any character. It would match “”, “h”, “hello”, “hhh”, etc. This is useful if there’s a value that can vary that you don’t care about, or if you don’t want to write out the whole path.
      2. If your cycles are each in different folders, should start with something like data/cycle*/
      3. See the example in the well_analysis file, this should help you figure out how to change the pattern so that it’ll match your files
      4. Here’s an article about f-strings if you’re confused about the curly braces in the string! 
   5. `project_name` - this should be the name of the subdirectory in `projects/` that you created
   6. `barcode_csv_name` - the name of the csv file which has the barcode library, which you added in step 2
   7. Thresholds - you will have to change these if the program isn’t detecting nuclei/cells/reads correctly. You won’t be able to determine whether they have to be changed until you have tried running the program. Once you run it, it should produce a nuclei.tif file, a cells.tif file, and a cells.csv file (among other files). If background spots are appearing on the nuclei file, for example, you should raise the dapi threshold. 
   8. `nucleus_area` probably shouldn’t have to be changed. Idk
   9. `DAPI_index` will change if you didn’t stain with DAPI on the first cycle. Remember that this is 0-indexed (the first cycle is zero). So if you stained on the 11th cycle, you would change this argument to 10. If you stained on all cycles, the number shouldn’t matter.
   10. `align_method` can be changed to DAPI if you stained with DAPI on all cycles. But it should work either way.
   11. `barcode_counts` shouldn’t have to be changed.
   
After running, all the output files should be in the `process_ipynb/` subdirectory of your experiment’s folder.
