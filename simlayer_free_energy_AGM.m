function simlayer_free_energy_AGM(SNR, varargin)
% SIMLAYER_FREE_ENERGY_AGM  Run Alternative Generative Model (AGM) simulations 
%   with whole-brain free energy analysis
%
% Use as
%   simlayer_free_energy(-20), where the argument is the SNR (db).
%
%   simlayer_free_energy(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * surf_dir - directory containing subject surfaces
%    * mri_dir - directory containing subject MRIs
%    * out_file - output file name (automatically generated if not
%    specified)
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * sim_patch_size - 5 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 5 (default) or interger - reconstruction patch size
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * invfoi - [10 30] (default) - frequency range for source inversion

% parse inputs
defaults = struct('surf_dir', '/data/pt_user-helbling_ticket017439/helbling/NormativeMEG/Data/Freesurfer6.0.0_Recons/',...
    'mri_dir', '/data/pt_user-helbling_ticket017439/helbling/NormativeMEG/Data/hcT1s/AB300686/',...
    'out_file', '', 'dipole_moment', 10, 'sim_patch_size', 5,...
    'reconstruct_patch_size', 5, 'nsims', 60, 'invfoi',[10 30]); % define default values

fid = [-12.2234 124.8297 2.6319; -82.8018 19.6687 -47.1408; 72.7654 28.4151 -36.0353]; % fiducial locations

params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

if isempty(params.rawfile)
    params.rawfile = '/data/pt_np-helbling/layer_opm_sim/opm_sim_data/sim_opm_custom_space_35_axis_1.mat';
end
% output directory
if isempty(params.out_path)
    params.out_path = '/data/pt_np-helbling/layer_opm_sim/results_opm_sim_space_35_axis_1';
end

if exist(params.out_path,'dir')~=7
    mkdir(params.out_path);
end

% file name for new file to work with
newfile = fullfile(params.out_path, 'opm_sim.mat');

% construct name for out file
if isempty(params.out_file)
    params.out_file=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
        invfoi(1),invfoi(2),SNR,params.dipole_moment);
end

spm('defaults', 'EEG');
spm_jobman('initcfg');

% make new file to work with
clear jobs
matlabbatch = [];
matlabbatch{1}.spm.meeg.other.copy.D = {params.rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = newfile;
spm_jobman('run', matlabbatch);

% load meshes for the AGM source space
white_AGM_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds11_white.surf.gii');
pial_AGM_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds11_pial.surf.gii');

white_AGM_mesh_gifti = gifti(white_AGM_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds11_white.surf_gifti.gii'), white_AGM_mesh_gifti.vertices, white_AGM_mesh_gifti.faces);
white_AGM_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds11_white.surf_gifti.gii');

pial_AGM_mesh_gifti = gifti(pial_AGM_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds11_pial.surf_gifti.gii'), pial_AGM_mesh_gifti.vertices, pial_AGM_mesh_gifti.faces);
pial_AGM_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds11_pial.surf_gifti.gii');

% load meshes for the GM source space
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');

white_mesh_gifti = gifti(white_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf_gifti.gii'), white_mesh_gifti.vertices, white_mesh_gifti.faces);
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf_gifti.gii');

pial_mesh_gifti = gifti(pial_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf_gifti.gii'), pial_mesh_gifti.vertices, pial_mesh_gifti.faces);
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf_gifti.gii');

Nmesh = 2;
simmeshes = {white_AGM_mesh, pial_AGM_mesh};
allmeshes = {white_mesh,pial_mesh};

% create smoothed meshes
for meshind = 1:length(allmeshes)
    spm_eeg_smoothmesh_mm(allmeshes{meshind},params.sim_patch_size);
end
% setup simulation - number of sources, list of vertices to simulate on
mesh_one = gifti(allmeshes{1});
nverts = size(mesh_one.vertices,1);
rng(0,'twister')
simvertind = randperm(nverts); %% random list of vertex indices to simulate sources on
Nsim = params.nsims; %% number of simulated sources

% for MSP
patchfilename = fullfile(params.out_path, 'temppatch.mat');
% number of patches as priors
Npatch = round(Nsim*1.5);
% so use all vertices that will be simulated on (plus a few more) as MSP priors
Ip = simvertind(1:Npatch);

% save priors
save(patchfilename,'Ip','simvertind');
% inversion method to use
methodnames = {'EBB','MSP'};

Nmeth = length(methodnames);

% inversion parameters
invwoi = [100 500];
% number of cross validation folds - just 1 here as we do not perform any
% cross validation here
Nfolds = 1;
% percentage of test channels in cross validation
ideal_pctest = 0;
% use all available spatial modes
ideal_Nmodes = [];


% all F values
allcrossF = zeros(Nmesh,Nsim,Nmesh,Nmeth);
allcrossVE = zeros(Nmesh,Nsim,Nmesh,Nmeth);

regfiles = {};
spatialmodesnames = {};

for meshind = 1:Nmesh
    regfile = fullfile(params.out_path, sprintf('opm_sim_%dcoreg.mat',meshind)); % this is the regfile name
    regfiles{meshind} = regfile;
    clear jobs
    matlabbatch = [];
    matlabbatch{1}.spm.meeg.other.copy.D = {params.rawfile};
    matlabbatch{1}.spm.meeg.other.copy.outfile = regfile;
    spm_jobman('run', matlabbatch);

    % add sensor displacements
    matlabbatch = [];
    matlabbatch{1}.spm.meeg.source.sensorshift.D = {regfile};
    matlabbatch{1}.spm.meeg.source.sensorshift.val = 1;
    matlabbatch{1}.spm.meeg.source.sensorshift.movewhat = 'independent';
    matlabbatch{1}.spm.meeg.source.sensorshift.meanshift = [0 0 0];
    matlabbatch{1}.spm.meeg.source.sensorshift.sdshift = [1 1 1];
    matlabbatch{1}.spm.meeg.source.sensorshift.meanangshift = [0 0 0];
    matlabbatch{1}.spm.meeg.source.sensorshift.sdangshift = [0 0 0];
    matlabbatch{1}.spm.meeg.source.sensorshift.outprefix = '';
    spm_jobman('run',matlabbatch);

    % coregister simulated dataset to reconstruction mesh
    matlabbatch = [];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {regfile};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,'Sim_MEG_hcT1.nii')};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {allmeshes{meshind}};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fid(1,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fid(2,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fid(3,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run', matlabbatch);

    % setup spatial modes for cross validation
    spatialmodesname = fullfile(params.out_path, sprintf('%d_testmodes.mat',meshind));
    [spatialmodesname,Nmodes,pctest] = spm_eeg_inv_prep_modes_xval(regfile, ideal_Nmodes, spatialmodesname, Nfolds, ideal_pctest);
    spatialmodesnames{meshind} = spatialmodesname;
end

% simulate sources on each mesh
for simmeshind = 1:Nmesh % choose mesh to simulate on
    simmesh = simmeshes{simmeshind};

    % coregister to correct mesh
    matlabbatch = [];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {newfile};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,'Sim_MEG_hcT1.nii')};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {simmesh};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fid(1,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fid(2,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fid(3,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run', matlabbatch);

    Dmesh = spm_eeg_load(filename);

    % now simulate sources on this mesh
    for s = 1:Nsim
        % get location to simulate dipole on this mesh
        simpos = Dmesh.inv{1}.mesh.tess_mni.vert(simvertind(s),:);
        prefix = sprintf('sim_mesh%d_source%d',simmeshind,s);

        % simulate source
        matlabbatch = [];
        matlabbatch{1}.spm.meeg.source.simulate.D = {filename};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = prefix;
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = invwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment params.sim_patch_size];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end

        spm_jobman('run', matlabbatch); 

        % now reconstruct onto all the meshes and look at cross val and F vals
        for meshind = 1:Nmesh

            % copy forward model from pial or white coregistered file
            simfilename = fullfile(params.out_path,sprintf('%sopm_sim.mat',prefix));
            sim = load(simfilename);
            reconcoreg = load(regfiles{meshind});
            sim.D.other = reconcoreg.D.other;
            sim.D.other.inv{1,1}.forward.loc = 0;
            D = sim.D;
            copyfile(fullfile(params.out_path, sprintf('SPMgainmatrix_opm_sim_%dcoreg_1.mat', meshind)), fullfile(params.out_path, sprintf('SPMgainmatrix_%sopm_sim_1.mat', prefix)));
            D.other.inv{1}.gainmat = sprintf('SPMgainmatrix_%sopm_sim_1.mat', prefix);
            save(simfilename,'D');

            % reconstruct using each method
            for methind = 1:Nmeth

                % do inversion of simulated data with this surface
                matlabbatch = [];
                matlabbatch{1}.spm.meeg.source.invertiter.D = {simfilename};
                matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = methodnames{methind};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = invwoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = invfoi;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = [1 1];
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = -params.reconstruct_patch_size;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesnames{simmeshind}};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 4;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
                matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
                matlabbatch{1}.spm.meeg.source.invertiter.crossval = [pctest Nfolds];
                spm_jobman('run', matlabbatch);

                % load inversion - get cross validation error end F
                Drecon = spm_eeg_load(simfilename);
                allcrossF(simmeshind,s,meshind,methind) = Drecon.inv{1}.inverse.crossF;
                allcrossVE(simmeshind,s,meshind,methind) = Drecon.inv{1}.inverse.VE;
            end
        end
        close all;
        delete(fullfile(params.out_path,'sim_mesh*'))
    end
end
if exist(fullfile(params.out_path,params.out_file),'file')
    save(fullfile(params.out_path,params.out_file),'allcrossF','allcrossVE','-append');
else
    save(fullfile(params.out_path,params.out_file),'allcrossF','allcrossVE');
end

% delete(fullfile(params.out_path,'SPMg*'))
% rmdir(fullfile(params.out_path, 'simp*'),'s')

