# Code repository for "Repeated evolution of underwater “rebreathing” in diving Anolis lizards"

This document covers how to generate the analyses and plots used in the paper
<br/>
<br/>
### There are two primary data sources:
  * Ethogram behavioural data coded in BORIS (http://www.boris.unito.it/pages/download.html)

  * Oxygen partial pressure and temperature data collected using a PyroScience OXB50 or OXR50 probe and FireSting GO2 meter (submergible O<sub>2</sub> microsensor setup)
<br/>

All intermediary files are provided in this repository. This means you skip straight to the analysis scripts if you do not intend to change/re-filter any of the raw data

### To run rebreathing behaviour analyses, open the 'behaviour' folder and follow the readme markdown file
  * This pipeline will:
   * Process raw BORIS .csv files
   * Generate trial and individual scores for all behaviours
   * Assess behaviour occurrence statistically
   * Plot behavioural scores (and generate figures 2 and 3 from the paper)

### To run oxygen trace analyses, open the 'oxygen' folder and follow its readme markdown file
  * This pipeline will:
    * Analyze and plot all rebreathing traces
    * Filter out poor quality trials (short duration, low oxygen readings)
    * Assess model fits, potential covariates
    * Plot best traces, oxygen consumption (figure 4)
    * Plot all traces with linear and exponential decline models (figure S8)
    * Plot all 'sham' trials (bubble inhalations/exhalations mimicked using a submerged syringe)