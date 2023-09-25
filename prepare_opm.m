function [] = prepare_opm(data_dir)
% use as prepare_opm('/data/pt_np-helbling/layer_opm_sim/opm_sim_data/')
% data_dir = '/data/pt_np-helbling/layer_opm_sim/opm_sim_data/';
% Add the right OPM toolbox! Careful, the developer's SPM version already
% contains an OPM toolbox, but without the option for multiple trials
addpath('/data/pt_np-helbling/layer_opm_sim/OPM')

if ~(exist(data_dir,'dir'))
    mkdir(data_dir)
end
cd(data_dir)

space = 25:10:55;
axis = [1 2 3];
for i = 1:length(space)
    for j = 1:length(axis)
        S = [];
        S.space = space(i);
        S.lead = 0;
        S.wholehead = 0;
        S.offset = 6.5; 
        S.axis = axis(j);
        S.nTrials = 200;
        S.fs = 200;
        S.nSamples = 200;
        S.cortex = fullfile('<FS_DIR>','/sub-01/','surf','ds_mid.surf.gii'); % down-sampled mid cortical surface
        S.sMRI = fullfile('<T1_DIR>','Sim_MEG_hcT1.nii');
        S.fname = sprintf('sim_opm_custom_space_%d_axis_%d',space(i),axis(j));
        spm_opm_sim(S);

        load(sprintf('sim_opm_custom_space_%d_axis_%d',space(i),axis(j)),'D')
        
        D.type = 'single';
        D.timeOnset = -0.5;
        D.condlist = 'stim';
        save(sprintf('sim_opm_custom_space_%d_axis_%d',space(i),axis(j)),'D')
    end
end

