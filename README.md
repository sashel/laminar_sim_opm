Laminar inferences using OPM-MEG
=======================

Accompanying code for my simulation paper on laminar inferences with OPMs. 

> S Helbling<br>
> **Inferring laminar origins of MEG signals with optically pumped magnetometers (OPMs): a simulation study**<br>
> doi: https://doi.org/10.1101/2023.08.20.554011

Based on code by JJ Bonaiuto (https://github.com/jbonaiuto/laminar_sim) for laminar inference using high-precision forward models
and by Tim Tierney (https://github.com/tierneytim/OPM/) for OPM-MEG array simulations.

> JJ Bonaiuto, HE Rossiter, SS Meyer, N Adams, S Little, MF Callaghan, F Dick, S Bestmann, GR Barnes<br>
> **Non-invasive laminar inference with MEG: Comparison of methods and source inversion algorithms**<br>
> NeuroImage (2017), http://www.sciencedirect.com/science/article/pii/S1053811917310145<br>
> (bioRxiv link: https://www.biorxiv.org/content/early/2017/11/30/147215)

> TM Tierney, S Mellor, GC O’Neill, N Holmes, E Boto, G Roberts, RM Hill, J Leggett, R Bowtell, MJ Brookes & GR Barnes<br> 
> **Pragmatic spatial sampling for wearable MEG arrays**<br> 
> Scientific Reports (2020), https://www.nature.com/articles/s41598-020-77589-8

## Requirements

* OPM tools: https://github.com/tierneytim/OPM/
* McNemarextest: http://www.mathworks.com/matlabcentral/fileexchange/6297
* myBinomTest: https://uk.mathworks.com/matlabcentral/fileexchange/24813/
* MEGSurfer: https://github.com/jbonaiuto/MEGsurfer
* BrewerMap: https://github.com/DrosteEffect/BrewerMap

## Usage

### Running simulations

% Create OPM-MEG arrays

% Run simulations across sensor densities whole-brain and ROI-based analysis 

Uses: 

    simlayer_free_energy(subjects(1), 1, [10 30], -20);
and
    simlayer_roi(subjects(1), 1, [10 30], -20);


% Run free energy - patch size simulations

    simlayer_free_energy_patch_size(subjects(1), 1, [10 30], -20);


% Run ROI - patch size simulations
  
    simlayer_roi_patch_size(subjects(1), 1, [10 30], -20);


### Analyzing results

% Plot whole brain and ROI stats for each simulation


## Support
Email saskia.helbling@gmail.com with any questions.
