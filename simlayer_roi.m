function simlayer_roi(SNR, varargin)
% SIMLAYER_ROI  Run simulations with ROI analysis
%
% Use as
%   simlayer_roi(-20), where the argument is the SNR (db).
%
%   simlayer_roi(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * surf_dir - directory containing subject surfaces
%    * mri_dir - directory containing subject MRIs
%    * out_path - output file path (automatically generated if not
%    specified)
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * sim_patch_size - 5 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 5 (default) or interger - reconstruction patch size
%    * invfoi - [10 30] (default) - frequency range for source inversion

% parse inputs
defaults = struct('surf_dir', '<FS_DIR>'...
    'mri_dir', '<T1_DIR>'...
    'rawfile', '', 'out_path', '', 'dipole_moment', 10, 'sim_patch_size', 5,...
    'out_path', '/data/pt_np-helbling/layer_opm_sim/results_opm_sim_space_55_axis_1',...
    'reconstruct_patch_size', 5, 'nsims', 60, 'npatch_factor',1.25, 'invfoi',[10 30]);  % define default values

params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

fid = [-12.2234 124.8297 2.6319; -82.8018 19.6687 -47.1408; 72.7654 28.4151 -36.0353]; % fiducial locations

if isempty(params.rawfile)
    params.rawfile = '/data/pt_np-helbling/layer_opm_sim/opm_sim_data/sim_opm_custom_space_55_axis_1.mat';
end

% output directory
params.out_path = fullfile(params.out_path,...
sprintf('f%d_%d_SNR%d_dipolemoment%d',params.invfoi(1),params.invfoi(2),SNR,params.dipole_moment));

if exist(params.out_path,'dir')~= 7
    mkdir(params.out_path);
end

% file names for the new files to work with
newfile = fullfile(params.out_path, 'opm_sim.mat');
greyregfile = fullfile(params.out_path, 'opm_sim_greycoreg.mat');

spm('defaults', 'EEG');
spm_jobman('initcfg');

% make new files to work with
clear jobs
matlabbatch = [];
matlabbatch{1}.spm.meeg.other.copy.D = {params.rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = newfile;
spm_jobman('run', matlabbatch);

% make new files to work with
clear jobs
matlabbatch = [];
matlabbatch{1}.spm.meeg.other.copy.D = {params.rawfile};
matlabbatch{1}.spm.meeg.other.copy.outfile = greyregfile;
spm_jobman('run', matlabbatch);

% load meshes
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf.gii');
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf.gii');
combine_surfaces(white_mesh, pial_mesh,fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2-ds_pial.hc_PDw2.surf.gii'))
pialwhite_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2-ds_pial.hc_PDw2.surf.gii');

white_mesh_gifti = gifti(white_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf_gifti.gii'), white_mesh_gifti.vertices, white_mesh_gifti.faces);
white_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_white.hc_PDw2.surf_gifti.gii');

pial_mesh_gifti = gifti(pial_mesh);
write_surf_gifti(fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf_gifti.gii'), pial_mesh_gifti.vertices, pial_mesh_gifti.faces);
pial_mesh = fullfile(params.surf_dir,'/sub-01/','surf','ds_pial.hc_PDw2.surf_gifti.gii');

simmeshes = {white_mesh,pial_mesh};
allmeshes = {white_mesh,pial_mesh,pialwhite_mesh};

% create smoothed meshes
for meshind = 1:length(allmeshes)
    spm_eeg_smoothmesh_mm(allmeshes{meshind},params.sim_patch_size);
end

% setup simulation - number of sources, list of vertices to simulate on
white = gifti(white_mesh); % to get number of vertices
nverts = size(white.vertices,1);
rng(1,'twister')
simvertind = randperm(nverts); % random list of vertex indices to simulate sources on
Nsim = 60; % number of simulated sources on each surface

% for MSP
Npatch = round(Nsim*params.npatch_factor);
Ip = [simvertind(1:Npatch) nverts+simvertind(1:Npatch)];
% save priors
patchfilename = fullfile(params.out_path, 'temppatch.mat');
save(patchfilename,'Ip');
methodnames = {'EBB','MSP'};
Nmeth = length(methodnames);

% inversion parameters
invwoi = [-500 500];
simwoi = [100 500];
baselinewoi = [-500 -100];
% number of cross validation folds
Nfolds = 1;
% percentage of test channels in cross validation
ideal_pctest = 0;
% use all available spatial modes
ideal_Nmodes = [];

% coregister simulated dataset to combined pial/white mesh
matlabbatch = [];
matlabbatch{1}.spm.meeg.source.headmodel.D = {greyregfile};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,'Sim_MEG_hcT1.nii')};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {pialwhite_mesh};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fid(1,:);
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fid(2,:);
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fid(3,:);
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
spm_jobman('run',matlabbatch);

% setup spatial modes for cross validation
spatialmodesname = fullfile(params.out_path, 'testmodes.mat');
[spatialmodesname,Nmodes,pctest] = spm_eeg_inv_prep_modes_xval(greyregfile, ideal_Nmodes, spatialmodesname, Nfolds, ideal_pctest);

greycoreg = load(greyregfile);

for simmeshind = 1:length(simmeshes)
    simmesh = simmeshes{simmeshind};

    % coregister to correct mesh
    spm_jobman('initcfg');
    matlabbatch = [];
    matlabbatch{1}.spm.meeg.source.headmodel.D = {newfile};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(params.mri_dir,'Sim_MEG_hcT1.nii')};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {simmesh};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fid(1,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fid(2,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fid(3,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    spm_jobman('run',matlabbatch);

    Dmesh = spm_eeg_load(newfile);

    % now simulate sources on this mesh
    for s = 1:Nsim
        % get location to simulate dipole on this mesh
        simpos = Dmesh.inv{1}.mesh.tess_mni.vert(simvertind(s),:);
        prefix = sprintf('sim_mesh%d_source%d',simmeshind, s);

        % simulate source
        matlabbatch = [];
        matlabbatch{1}.spm.meeg.source.simulate.D = {newfile};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = prefix;
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = simwoi;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [params.dipole_moment, params.sim_patch_size];
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
        if abs(params.dipole_moment)>0
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;
        else
            matlabbatch{1}.spm.meeg.source.simulate.isSNR.whitenoise = 100;
        end
        spm_jobman('run', matlabbatch);


        % copy forward model from pial/white coregistered file
        simfilename = fullfile(params.out_path,sprintf('%sopm_sim.mat',prefix));
        sim = load(simfilename);
        sim.D.other = greycoreg.D.other;
        sim.D.other.inv{1,1}.forward.loc = 0;

        D = sim.D;
        copyfile(fullfile(params.out_path, 'SPMgainmatrix_opm_sim_greycoreg_1.mat'), fullfile(params.out_path, sprintf('SPMgainmatrix_%sopm_sim_1.mat', prefix)));
        D.other.inv{1}.gainmat = sprintf('SPMgainmatrix_%sopm_sim_1.mat', prefix);
        save(simfilename,'D');

        % reconstruct using each method
        for methind = 1:Nmeth
            method = methodnames{methind};

            % do inversion of simulated data with this surface
            matlabbatch = [];
            matlabbatch{1}.spm.meeg.source.invertiter.D = {simfilename};
            matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = method;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = invwoi;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = invfoi;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {patchfilename};
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = [1 1]; 
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm  = -params.reconstruct_patch_size; 
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 4;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
            matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
            matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
            matlabbatch{1}.spm.meeg.source.invertiter.crossval = [pctest Nfolds];
            spm_jobman('run',matlabbatch);

            D = spm_eeg_load(simfilename);
            goodchans = D.indchantype('MEG','good');
            M = D.inv{1}.inverse.M;
            U = D.inv{1}.inverse.U{1};
            T = D.inv{1}.inverse.T;
            It   = D.inv{1}.inverse.It;
            times  = D.inv{1}.inverse.pst;
            Dgood = squeeze(D(goodchans,It,:));
            ntrials = size(Dgood,3);
            m = export(gifti(D.inv{1}.mesh.tess_ctx),'patch');
            GL      = spm_mesh_smooth(m);

            wois = [baselinewoi; simwoi];
            for w = 1:size(wois,1)
                woi = wois(w,:);

                fwhm = max(diff(woi),8);
                t    = exp(-4*log(2)*(times(:) - mean(woi)).^2/(fwhm^2));
                t    = t/sum(t);

                % get frequency space and put PST subspace into contrast (W -> T*T'*W)
                %--------------------------------------------------------------------------
                wt = 2*pi*times(:)/1000;
                W  = [];
                for f = invfoi(1):invfoi(end)
                    W = [W sin(f*wt) cos(f*wt)];
                end
                W  = diag(t)*W;
                W  = spm_svd(W,1);
                TW     = T'*W;
                TTW    = T*TW;

                woi_vals = zeros(nverts*2,ntrials);
                MU = M*U;
                for i = 1:ntrials
                    MUd1 = MU*squeeze(Dgood(:,:,i));
                    Y     = sum((MUd1*TTW).^2,2);
                    woi_vals(:,i) = spm_mesh_smooth(GL,Y,8);
                end

                woi_dir = fullfile(params.out_path, ['t' num2str(woi(1)) '_' num2str(woi(2))]);
                if exist(woi_dir,'dir')~= 7
                    mkdir(woi_dir);
                end
                delete(fullfile(woi_dir,'*'));
                out_filename = fullfile(woi_dir, sprintf('%sopm_sim_1_t%d_%d_f%d_%d', prefix, woi(1), woi(2), invfoi(1), invfoi(2)));
                write_metric_gifti(out_filename,woi_vals);

                % split pial and grey sources
                split_inversion_results(woi_dir);
            end

            baseline_dir = fullfile(params.out_path, ['t' num2str(baselinewoi(1)) '_' num2str(baselinewoi(2))]);
            woi_dir = fullfile(params.out_path, ['t' num2str(simwoi(1)) '_' num2str(simwoi(2))]);

            % load all pial data from wois
            pial_woi_trials = gifti(fullfile(woi_dir,sprintf('pial_%sopm_sim_1_t%d_%d_f%d_%d.gii', prefix, simwoi(1), simwoi(2), invfoi(1), invfoi(2))));
            pial_baseline_trials = gifti(fullfile(baseline_dir,sprintf('pial_%sopm_sim_1_t%d_%d_f%d_%d.gii', prefix, baselinewoi(1), baselinewoi(2), invfoi(1), invfoi(2))));

            % load all white matter data from wois
            white_woi_trials = gifti(fullfile(woi_dir,sprintf('white_%sopm_sim_1_t%d_%d_f%d_%d.gii', prefix, simwoi(1), simwoi(2), invfoi(1), invfoi(2))));
            white_baseline_trials = gifti(fullfile(baseline_dir,sprintf('white_%sopm_sim_1_t%d_%d_f%d_%d.gii', prefix,baselinewoi(1), baselinewoi(2), invfoi(1), invfoi(2))));

            % save pial diff
            pial_diff = pial_woi_trials.cdata(:,:)-pial_baseline_trials.cdata(:,:);
            white_diff = white_woi_trials.cdata(:,:)-white_baseline_trials.cdata(:,:);

            write_metric_gifti(fullfile(params.out_path, sprintf('pial.%s.%s.gii',method,prefix)), pial_diff);
            write_metric_gifti(fullfile(params.out_path, sprintf('white.%s.%s.gii',method,prefix)), white_diff);

        end
        close all
        delete(fullfile(params.out_path, sprintf('sim_mesh%d_source%d*.*',simmeshind, s)));
    end
end
