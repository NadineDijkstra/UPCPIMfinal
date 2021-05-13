function GroupMeanNifti(cfg)
% function GroupMeanNifti(cfg)
%
% cfg.subjects = cell array of subject names
% cfg.root     = root directory
% cfg.contrast = name of nifti contrast to average
% cfg.dir      = name of dir in root/subjectName where contrasts are
%
% Results will be saved in root/GroupResults/cfg.dir

outputDir = fullfile(cfg.root,'GroupResults',cfg.dir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

if ~isfield(cfg,'ttest'); cfg.ttest = false; end

%% Get map per subject 
nsubjects = length(cfg.subjects);

niftis = [];

for sub = 1:nsubjects
    
    [hdr,tmp] = read_nii(fullfile(cfg.root,cfg.subjects{sub},cfg.dir,cfg.contrast));
    
    niftis = cat(5,niftis,tmp); clear tmp   
    
end

% average over subjects
avgMap = squeeze(mean(niftis,5));
write_nii(hdr,avgMap,fullfile(outputDir,['GA_' cfg.contrast '.nii.gz']));

% make a mask
if cfg.mask
    mask  = sum(niftis~=0,5);
    write_nii(hdr,mask,fullfile(outputDir,['mask_' cfg.contrast '.nii.gz']));
end

% do a t-test against chance
if cfg.ttest   
    [~,p,~,stats] = ttest(niftis,0.5,'tail','right','dim',5);
    write_nii(hdr,1-p,fullfile(outputDir,['pval_' cfg.contrast '.nii.gz']));
    write_nii(hdr,stats.tstat,fullfile(outputDir,['tval_' cfg.contrast '.nii.gz']));    
end

