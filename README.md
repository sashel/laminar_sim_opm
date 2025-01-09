Laminar inferences using OPM-MEG
=======================

Accompanying code for my simulation paper on laminar inferences with OPMs. 

> S Helbling<br>
> **Inferring laminar origins of MEG signals with optically pumped magnetometers (OPMs): a simulation study**<br>
> Now out in Imaging Neuroscience: https://doi.org/10.1162/imag_a_00410
> bioRxiv doi: https://doi.org/10.1101/2023.08.20.554011

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

* OPM tools: https://github.com/sashel/OPM/ (forked from https://github.com/tierneytim/OPM/, contains the option to simulate multiple trials)
* McNemarextest: http://www.mathworks.com/matlabcentral/fileexchange/6297
* myBinomTest: https://uk.mathworks.com/matlabcentral/fileexchange/24813/
* MEGSurfer: https://github.com/jbonaiuto/MEGsurfer
* BrewerMap: https://github.com/DrosteEffect/BrewerMap

## Usage

### Construct the OPM-MEG arrays

Create OPM-MEG arrays across inter-sensor distances, number of measurement axes for a default scalp-sensor offset of 6.5 mm. Simulated data consists of 200 trials of 1 second each, sampled at a rate of 200/sec.

    prepare_opm(<PROJECT_DIR>/DATA_DIR>)

Create OPM-MEG arrays across varying larger scalp-sensor offsets of 20, 30 and 40 mm at an inter-sensor distance of 35 mm and a single measurement axis.  

    prepare_opm_offsets(<PROJECT_DIR>/DATA_DIR>)

### Running simulations

Run simulations across SNRs for an OPM-MEG sensor array with an inter-sensor distance of 35 mm and one measurement axis using the whole-brain and ROI-based analyses 

    run_opm_sim(<PROJECT_FOLDER>,35,1,[-5,-10,-20,-30,-40,-50)

Run simulations across inter-sensor distances ([25 35 45 55] mm) for a sensor array with radial measurement axes

    run_opm_sim(<PROJECT_FOLDER>,[25 35 45 55],1,[-5,-10,-20,-30,-40])

Run simulations across number of measurement axes ([1 2 3]) at an inter-sensor distance of 55 mm 

    run_opm_sim(<PROJECT_FOLDER>,55,[1 2 3],[-5,-10,-20,-30,-40])

Run simulations across varying scalp-sensor offsets ([20 30 40] mm) for a radial axes OPM array at an inter-sensor distance of 35 mm 

    run_opm_sim_offset(<PROJECT_FOLDER>,35,1,[-5,-10,-20,-30,-40],[20 30 40])

Run simulations across varying co-registration errors ([1 2 3 4] mm std of fiducial errors)

    run_opm_sim_coreg(<PROJECT_FOLDER>,35,1,-10,[1 2 3 4])

Run simulations for congruent and incongruent patch sizes

    run_opm_sim_patch_sizes('/data/pt_np-helbling/layer_opm_sim/',35,1,[-5,-10,-20,-30,-40])

Run laminar inference simulations in the presence of confounding internal brain noise sources

    run_opm_sim_noise_sources('/data/pt_np-helbling/layer_opm_sim/',25,1,[-5])

Run simulations with an Alternative Generative Model (AGM)

    run_opm_sim_AGM('/data/pt_np-helbling/layer_opm_sim/',55,1,[-5,-10,-20,-30,-40,-50])


### Analyzing results
     
Plot classification accuracy and bias across SNRs (Fig.2)

    plot_snr_classification_performance('correct',35,1)
    plot_snr_classification_performance('pial',35,1)

Plot classification results across inter-sensor distances (Fig.3)

    plot_space_classification_performance('correct',1)
    plot_space_classification_performance('pial',1)

Plot classification results across number of measurement axes (Fig.4)

    plot_axes_classification_performance('correct',35)
    plot_axes_classification_performance('pial',35)

Plot classification results across scalp-sensor offsets (Fig. 5)

    plot_offset_classification_performance('correct',55,1)
    plot_offset_classification_performance('pial',55,1)

Plot classification results across co-registration errors (Fig. 6)

    plot_coreg_classification_performance('correct',35)
    plot_coreg_classification_performance('pial',35)

Plot classification results for congruent and incongruent patch sizes (Fig. 7)

    plot_free_energy_patch_size_classification_performance()

## Support
Email saskia.helbling@gmail.com with any questions.
