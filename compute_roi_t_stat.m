function tstat=compute_roi_t_stat(file_prefix, pial_meshname, wm_meshname,...
    varargin)
% Compute the t-statistic for an ROI

% Parse inputs
defaults = struct('mapType', 'link', 'recompute', false, 'delete', false, 'origPial', '',...
    'origWhite', '', 'nsims', 60);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

[filepath, filename, ~]=fileparts(file_prefix);
% Files containing t-statistics and pial-wm diff
pial_t_filename=fullfile(filepath, sprintf('pial.%s.t.gii', filename));
wm_t_filename=fullfile(filepath, sprintf('white.%s.t.gii', filename));
pial_wm_diff_filename=fullfile(filepath, sprintf('pial-white.%s.gii', filename));
pial_white_map=map_pial_to_white(wm_meshname, pial_meshname, ...
        'mapType', params.mapType, 'origPial', params.origPial, ...
        'origWhite', params.origWhite);
    
% If tstat files or pial_wm file do not exist or recomputing
if exist(pial_t_filename,'file')~=2 || exist(wm_t_filename,'file')~=2 || exist(pial_wm_diff_filename, 'file')~=2 || params.recompute
    fprintf('Recomputing t-stat for %s\n', file_prefix);
    
    % Load pial and white matter surface (woi-baseline)
    pial_diff=gifti(fullfile(filepath,sprintf('pial.%s.gii',filename)));
    wm_diff=gifti(fullfile(filepath,sprintf('white.%s.gii',filename)));
            
    % Run pial surface t-test
    varpop=nanvar([pial_diff.cdata(:,:) wm_diff.cdata(:,:)],[],2);
    [tstat,~]=ttest_corrected(pial_diff.cdata(:,:)','correction',.01*max(varpop));
    pial_tvals=tstat';
    write_metric_gifti(pial_t_filename, pial_tvals);
    
    % Run wm surface t-test
    [tstat,~]=ttest_corrected(wm_diff.cdata(:,:)','correction',.01*max(varpop));
    wm_tvals=tstat';
    write_metric_gifti(wm_t_filename, wm_tvals);
            
    % Compute pial-white difference
    pial_wm_diff=abs(pial_diff.cdata(:,:))-abs(wm_diff.cdata(pial_white_map,:));
    write_metric_gifti(pial_wm_diff_filename, pial_wm_diff);
else % Otherwise load data from files
    x=gifti(pial_t_filename);
    pial_tvals=x.cdata(:);
    x=gifti(wm_t_filename);
    wm_tvals=x.cdata(:);
    x=gifti(pial_wm_diff_filename);
    pial_wm_diff=x.cdata(:,:);
end

pial_threshold=prctile(pial_tvals(~isinf(pial_tvals)),75);
% Create pial and white masks and mapped white mask
pial_mask=find(pial_tvals>pial_threshold & ~isinf(pial_tvals));

wm_threshold=prctile(wm_tvals(~isinf(wm_tvals)),75);
mapped_wm_tvals=wm_tvals(pial_white_map);
mapped_wm_mask=find(mapped_wm_tvals>wm_threshold & ~isinf(mapped_wm_tvals));
% lower threshold to 50th percentile if needed
if isempty(mapped_wm_mask)
    pial_threshold=prctile(pial_tvals(~isinf(pial_tvals)),50);
    % Create pial and white masks and mapped white mask
    pial_mask=find(pial_tvals>pial_threshold & ~isinf(pial_tvals));
    wm_threshold=prctile(wm_tvals(~isinf(wm_tvals)),50);
    mapped_wm_tvals=wm_tvals(pial_white_map);
    mapped_wm_mask=find(mapped_wm_tvals>wm_threshold & ~isinf(mapped_wm_tvals));
    sprintf('Lower threshold in file %s', filename)
end
mask=union(pial_mask, mapped_wm_mask);
        
% Get mean pial-wm in ROI
pial_wm_roi_diff=mean(pial_wm_diff(mask,:));
% Perform ROI t-stat
[tstat,~]=ttest_corrected(pial_wm_roi_diff','correction',1000*var(pial_wm_roi_diff));
fprintf('ROI size=%d, tstat=%.2f\n',length(mask),tstat);

if params.delete
    if ~isempty(dir(fullfile(filepath, 'simprior*')))
        rmdir(fullfile(filepath, 'simprior*'),'s')
    end
    delete(fullfile(filepath, sprintf('pial.%s.dat', filename)),fullfile(filepath, sprintf('white.%s.dat', filename)), fullfile(filepath, sprintf('pial-white.%s.dat', filename)))
    delete(fullfile(filepath, sprintf('pial.%s.t.dat', filename)),fullfile(filepath, sprintf('white.%s.t.dat', filename)))
    delete(fullfile(filepath, 'SPMgainmatrix*mesh*'))
end

