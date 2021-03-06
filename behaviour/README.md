# Rebreathing behaviour analysis pipeline

### About the raw data (videos)
These analyses are based on video data collected in Colombia, Costa Rica, Ecuador, the Dominican Republic, Haiti, Jamaica, and Mexico. Videos were taken of semi-aquatic and non-aquatic lizards (primarily anoles) under natural or experimental dive conditions. Lizards performed a variety of behaviours while submerged.

<br/>

### Analysis of videos
All videos were analyzed using BORIS software (http://www.boris.unito.it/pages/download.html). See the supplement to this paper for a full ethogram.

If you are interested in viewing the data in BORIS, we have provided a .boris file in this repository (rebreathing_final.boris) that contains all trials.

All trials were exported as .csv files with the trial name identifier as the filename. In this repository, they can be found in the "boris_csvs" folder

Since our video collection comprises several hundred files totalling around 700 GB, video files will only be made available upon request, and are not part of this repository.

<br/>

### BORIS file processing

The csv files output by BORIS have a regular structure and were processed using a custom Python script (version 3.7.4). The script is extensively commented.

To re-extract data from the BORIS files, run the "boris_py.py" script. Paths are tailored to the repository, and will need to be adjusted if used separately. The script makes use of file with all the names of trials to be processed ("boris_names.txt"). This file should contain either the absolute paths to the files, or just the names (and the os.chdir command on line 70 should reflect the folder in which the files are stored).
To run the script, use a command like:

```
  python3 boris_py.py
```
The script creates a unified csv ("boris_summary.csv") that includes numbers of occurrences of each event in each trial, trial duration, event rates, 'rebreathing period' (interval between first and last rebreathing event in trials with >1 rebreathing event) and intervals between all rebreathing events (also in trials with >1 rebreathing event). Intervals were not a part of the paper, and are saved to a separate file ("intervals.csv").

In some videos, important regions of a lizard's body were obscured (by an experimentor's hand/lizard movements/etc). Durations of these periods are calculated by the script and are removed from consideration for the relevant behaviours.

We noted some rebreathing air movements that did not extend beyond the nares. These were not thought to be a equivalent to other performances of the rebreathing behaviour and were not counted towards rebreathing rates. These events were coded 'lbn' (low-by-narial) and scores are still tabulated per trial and per individual; however they were not incorporated into the interval scores. The comments field of BORIS events was used to denote 'lbn' reinhales as well the location of bubbles during full performances of rebreathing (top-of-head, toh; side-of-head, soh; high-by-narial, hbn).
<br/>
<br/>
### Processing of BORIS summary file, generation of individual scores
Next, we summarized data by individual. This was done in two ways (results were equivalent with both methods). 
* Best performance for each individual only
* Performance averaged across all trials

We defined 'best' as having the greatest number of rebreathing events; if two (or more) trials had the same number of events (especially relevant for lizards for which no rebreathing was observed), the trial with the longest duration was selected).

Individual summaries (and downstream analyses) were done using R (version 3.6.1). All remaining scripts are best used interactively in RStudio (but can be run from the command line, i.e,
```
source("generate_individual.R")
```
<br/>

#### R Script pathway
1. generate_individual.R
   * matches BORIS summary data with individual lizard/trial metadata (e.g., mass, SVL, sex, trial ID, etc.) using the "behaviour_metadata.csv" file and a custom R function (matchstick)
   * generates matched .csv file including all individuals
   * then filters out juveniles, trials below 15 seconds in duration
   * calculates and ouputs trial-averaged individual scores and best performance individual scores ("mean_individual.csv", "bp_individual.csv" respectively)
     * a description of the variables of interest/other information found in this dataframe can be viewed in "variable_description.csv"

2. At this stage, you can now analyze trial-averaged or best performance results. 
   * Step 3 gives details on how to run the "analyze_behaviour_final.R" script, which will do the majority of analyses included in the paper and print Figures 2 and 3. The majority of model probabilties and statistic values cited under the "Underwater Rebreathing Behavior in _Anolis_" heading in results were generated by this script.
   * Step 4 discusses the "analyze_behaviour_mean_final.R" script that was used to generate supplemental results and figures figure S2a-f/associated caption results.
   * Both scripts make use of 3 custom functions: rb_assess generates species averages from an individual dataset and runs a PGLS model to test significance between aquatic and non-aquatic groups; var_endings adds endings to an input vector; and matchstick matches up metadata with a focal dataframe using IDs in a shared column

3. Best performance analysis
   * Use "analyze_behaviour_final.R"
   * Script takes "bp_individual.csv", filters out individuals from species with n \< 3, then matches up each individual with mass, snout-vent length (SVL), habitat, and other individual data
   * Script then aggregates data by species, then runs phylogenetic generalized least squares analyses using the _caper_ R package.
     * The phylogeny used for this was taken from [Poe 2017](https://doi.org/10.1093/sysbio/syx029); the tree was trimmed to just our focal taxa.
   * Script will output figures 2 and 3 to file
   * Script is extensively commented

4. Trial-averaged analysis
   * Use "analyze_behaviour_mean_final.R"
   * Conducts majority of the analyses conducted in the 'best performance' version, with the exeception of bubble and max reinhales (since both of these are based on 'best performance' only for each individual).
   * Results from these analyses are essentially fully concordant with our best performance results
   * Outputs figure S2a-f
   * Outputs results quoted in caption for figure S2

5. Bootstrap subsample analysis
   * Due to our uneven sampling, we decided to do a subsampling analysis to ensure that over-sampling of semi-aquatic species was not responsible for our results
   * We randomly subsampled all species to the lowest n we allowed into our principle analyses (n=3)
   * This was accomplished using the "subsample_bootstrap.R" script
   * Best performance and mean data sets were subsampled 10,000 times each
   * Script is split into best performance and mean results
     * Can use aq_subsample function defined in script to subsample either dataset
   * Generates figure S1, results in caption
