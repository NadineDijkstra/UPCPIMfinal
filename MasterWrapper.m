% add the relevant paths
restoredefaultpath;
addpath('/spm12')
addpath('Analyses')
addpath('Analyses/Utilities')
addpath('Analyses/Subjects')
addpath('/fieldtrip');
ft_defaults;

% root
root = '/UPCPIM';
if ~exist(root,'dir'); mkdir(root); end
cd(root)

dataDir = fullfile(root,'raw'); % raw data


subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10',...
    'S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21',...
    'S22','S23','S24','S25','S26','S27','S28','S29','S30','S31','S32','S35','S36','S37'};
nsubjects = length(subjects);


%% 0. Behavioural

for sub = 1:nsubjects
    subjectID = subjects{sub};
    
    cfg = [];
    cfg.subjectID = subjectID;
    cfg.dataDir   = dataDir;
    cfg.root      = root;
    cfg.plot      = 0;
    cfg.outputDir = 'Behaviour';
    
    Behaviour_analysis(cfg);
    
end

%% 1. Preprocessing
        
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    fprintf('Running subject %s \n',subjectID);
    
    eval(sprintf('Subject%s',subjectID))
    runs = [subjectdata.IMruns subjectdata.UPCPruns subjectdata.struct];
    
    % 1.1 performs dicom import, realignemnt, co-registration, segmentation and
    % renaming of key files
    cfg  = [];
    cfg.dicom_dir = fullfile(dataDir,subjectID,'ses-mri01');
    cfg.run_nr    = runs;
    cfg.nifti_dir = fullfile(root,subjectID,'Niftis');
    
    
    PreProcessDICOM(cfg);    
    
end

% 1.2 check movement
subjectID = 'S04';
cfg   = [];
cfg.nifti_dir = fullfile(root,subjectID,'Niftis');

CheckMovement(cfg)

%% 2. Create DARTEL template for normalisation

% make template
cfg = [];
cfg.root       = root;
cfg.nifti_dir  = 'Niftis';
cfg.dartel_dir = fullfile(root,'GroupResults','DARTEL_all');
cfg.subjects   = subjects;

DartelCreateTemplate(cfg)

%% 3. Normalise functionals using DARTEL

cfg = [];
cfg.root       = root;
cfg.smooth     = 6;
cfg.nifti_dir  = 'Niftis';
cfg.dartel_dir = fullfile(root,'GroupResults','DARTEL_all');
cfg.subjects   = subjects(15:20);

DartelNormalisation(cfg)

% reslice wm masks into normalized space
wm_mask = 'path to your csf mask';

resliceParameters = struct(...
    'prefix', 'sw',...
    'mask', 1,...
    'interp', 4, ...
    'wrap', [0 0 0], ...
    'which', [2 1]);

resFlags = struct(...
    'interp', resliceParameters.interp,... % interpolation type
    'wrap', resliceParameters.wrap,...     % wrapping info (ignore...)
    'mask', resliceParameters.mask,...     % masking (see spm_reslice)
    'which', 1,...                  % what images to reslice
    'mean', 0);                     % write mean image


mean_file = str2fullfile(fullfile(root,'S01','Niftis'),'swrf*.nii');
mean_file = mean_file{100};

spm_reslice({mean_file,wm_mask}, resFlags);


%% 4. First level

% 4.1 per trial for decoding
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    
    if ~exist(fullfile(root,subjectID,'FirstLevel_Trials',...
            'beta_0001.nii'),'file')
        
        eval(sprintf('Subject%s',subjectID))
        
        % all trials in one GLM
        cfg = [];
        cfg.TR           = 1;
        cfg.prefix       = 'swrf';
        cfg.run_nr{1}    = subjectdata.IMruns;
        cfg.run_nr{2}    = subjectdata.UPCPruns;
        cfg.logfiles{1}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*IM_' subjectID '*']);
        cfg.logfiles{2}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*UPCP_' subjectID '*']);
        cfg.identifier   = ['sub0' subjectID(2:end)];%'945416';
        cfg.nifti_dir = fullfile(root,subjectID,'Niftis');
        
        cfg.outputDir    = fullfile(root,subjectID,'FirstLevel_Trials');
        
        FirstLevelTrials(cfg)
    end
end

% 4.2 per condition for PPI
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    
    if ~exist(fullfile(root,subjectID,'GLM_Conditions',...
            'beta_0001.nii'),'file')
        
        eval(sprintf('Subject%s',subjectID))
        
        % conditions
        cfg = [];
        cfg.prefix       = 'swrf';
        cfg.TR           = 1;
        cfg.run_nr{1}    = subjectdata.IMruns; % which runs are IM
        cfg.run_nr{2}    = subjectdata.UPCPruns; % which runs are UPCP
        cfg.identifier   = ['sub0' subjectID(2:end)]; %'945416';
        cfg.nifti_dir    = fullfile(root,subjectID,'Niftis');
        cfg.logfiles{1}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*IM_' subjectID '*']);
        cfg.logfiles{2}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*UPCP_' subjectID '*']);
        cfg.outputDir    = fullfile(root,subjectID,'GLM_Conditions');
        
        GLMConditions(cfg)
    end
end

% 4.3 make contrasts
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    cfg = [];
    cfg.root      = root;
    cfg.subjectID = subjectID;
    cfg.model_dir = 'GLM_Conditions';
    cfg.tcontrasts = {[-1 1 0],... % cp vs up
        [0 1 -1],... % cp vs im
        [-1 0 1],... % im vs up
        [1 1 1]}; % average over all
    cfg.tnames     = {'CP vs UP','CP vs IM','IM vs UP','all'};
    cfg.fcontrasts = {[1 0 0; 0 1 0; 0 0 1]};
    cfg.fnames     = {'effects of interest'};
    
    MakeContrasts(cfg)
end


%% 5. Searchlight decoding analyses within conditions
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    
    if ~exist(fullfile(root,subjectID,'Within','consciousDecoding.nii.gz'),'file')
        % 5.1 within CP
        cfg = [];
        cfg.gamma = 0.01;
        cfg.radius = 4;
        cfg.betas = fullfile(root,subjectID,'FirstLevel_Trials');
        cfg.nifti_dir = fullfile(root,subjectID,'Niftis');
        cfg.condition = 'conscious';
        cfg.outputDir = 'Within';
        
        DecodingSearchlight(cfg)
    end    
    
    if ~exist(fullfile(root,subjectID,'Within','unconsciousDecoding.nii.gz'),'file')
        % 5.2 within UP
        cfg = [];
        cfg.gamma = 0.01;
        cfg.radius = 4;
        cfg.betas = fullfile(root,subjectID,'FirstLevel_Trials');
        cfg.nifti_dir = fullfile(root,subjectID,'Niftis');
        cfg.condition = 'unconscious';
        cfg.outputDir = 'Within';
        
        DecodingSearchlight(cfg)
    end
    
    
    if ~exist(fullfile(root,subjectID,'Within','imageryDecoding.nii.gz'),'file')
        % 5.3 within IM
        cfg = [];
        cfg.gamma = 0.01;
        cfg.radius = 4;
        cfg.betas = fullfile(root,subjectID,'FirstLevel_Trials');
        cfg.nifti_dir = fullfile(root,subjectID,'Niftis');
        cfg.condition = 'imagery';
        cfg.outputDir = 'Within';
        
        DecodingSearchlight(cfg)
    end
end

%% 6. Searchlight cross-decoding analyses
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    
    % 6.1a conscious to imagery
    cfg = [];
    cfg.gamma = 0.01;
    cfg.radius = 4;
    cfg.dataDir      = 'Within'; % data
    cfg.outputDir    = 'Cross';
    cfg.condition{1} = 'conscious';
    cfg.condition{2} = 'imagery';
    cfg.betas        = fullfile(root,subjectID,'FirstLevel_Trials');
    cfg.nifti_dir    = fullfile(root,subjectID,'Niftis');
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'conscious_to_imagery.nii.gz'),'file')
        CrossDecodingSearchlight(cfg)
    end
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'imagery_to_conscious.nii.gz'),'file')
        cfg.condition{1} = 'imagery';
        cfg.condition{2} = 'conscious'; % 6.1b and vice versa
        
        CrossDecodingSearchlight(cfg)
    end
    
    % 6.2a unconscious to imagery
    cfg = [];
    cfg.gamma = 0.01;
    cfg.radius = 4;
    cfg.dataDir      = 'Within'; % data
    cfg.outputDir    = 'Cross';
    cfg.condition{1} = 'unconscious';
    cfg.condition{2} = 'imagery';
    cfg.betas        = fullfile(root,subjectID,'FirstLevel_TrialsFinDARTEL');
    cfg.nifti_dir    = fullfile(root,subjectID,'Niftis');
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'unconscious_to_imagery.nii.gz'),'file')
        CrossDecodingSearchlight(cfg)
    end
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'imagery_to_unconscious.nii.gz'),'file')
        cfg.condition{1} = 'imagery';
        cfg.condition{2} = 'unconscious'; % 6.2b and vice versa
        
        CrossDecodingSearchlight(cfg)
    end
    
    % 6.3a conscious to unconscious
    cfg = [];
    cfg.gamma = 0.01;
    cfg.radius = 4;
    cfg.dataDir      = 'Within'; % data
    cfg.outputDir    = 'Cross';
    cfg.condition{1} = 'conscious';
    cfg.condition{2} = 'unconscious';
    cfg.betas        = fullfile(root,subjectID,'FirstLevel_TrialsFinDARTEL');
    cfg.nifti_dir    = fullfile(root,subjectID,'Niftis');
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'conscious_to_unconscious.nii.gz'),'file')
        CrossDecodingSearchlight(cfg)
    end
    
    if ~exist(fullfile(root,subjectID,cfg.outputDir,'unconscious_to_conscious.nii.gz'),'file')
        cfg.condition{1} = 'unconscious';
        cfg.condition{2} = 'conscious'; % 6.3b and vice versa
        
        CrossDecodingSearchlight(cfg)
    end
end

%% 7. Psycho-physiological interaction

% 7.1 get the VOI
VOIcentre  = [-56 -61 -6]; % in MNI coordinates point of overlap within-decoding all three conditions
fCon       = 1; % effects of interest F-contrast
tCon       = 4; % T-contrast to find peak activity (avg over conds)
VOIname    = 'LOC';
pthresh    = 0.05;

for sub = 1:nsubjects
    spm_file = fullfile(root,subjects{sub},'GLM_Conditions','SPM.mat');
    getVOI(spm_file,fCon,tCon,pthresh,VOIname,VOIcentre) ;
end

% 7.2. create VOI heat-map
VOI = [];
for sub = 1:nsubjects
    [hdr,voi] = read_nii(fullfile(root,subjects{sub},'GLM_Conditions','VOI_LOC_mask.nii'));
    VOI = cat(4,VOI,voi);
    clear voi
end
VOI = sum(VOI,4);
write_nii(hdr,VOI,fullfile(root,'GroupResults','GLM_Conditions',['VOI_' VOIname '_heat.nii']))


% 7.3 create PPI variable
VOIname  = 'LOC';
contrast = [1 1 -2; 2 1 1; 3 1 1];%[1 1 1; 2 1 1; 3 1 -2];
conName  = 'conscious';%'input';
PPIname = [VOIname '_' conName];

for sub = 1:nsubjects
    
    spm_file = fullfile(root,subjects{sub},'GLM_Conditions','SPM.mat');
    voi      = fullfile(root,subjects{sub},'GLM_Conditions',['VOI_' VOIname '_1.mat']);
    
    if ~exist(fullfile(root,subjects{sub},'GLM_Conditions',['PPI_' PPIname '.mat']),'file')
        getPPI{1}.spm.stats.ppi.spmmat = {spm_file};
        getPPI{1}.spm.stats.ppi.type.ppi.voi = {voi};
        getPPI{1}.spm.stats.ppi.type.ppi.u = contrast;
        getPPI{1}.spm.stats.ppi.name = PPIname;
        getPPI{1}.spm.stats.ppi.disp = 1;
        
        spm_jobman('run',getPPI)
    end
end

% 7.4 create PPI GLM folder and copy PPI variable there
folderName = ['GLM_Conditions_PPI_' conName];
PPIname    = ['PPI_' VOIname '_' conName '.mat'];
for sub = 1:nsubjects
    
    PPIdir = fullfile(root,subjects{sub},folderName);
    if ~exist(PPIdir,'dir')
        mkdir(PPIdir);
    end
    
    PPIfile = fullfile(root,subjects{sub},'GLM_Conditions',PPIname);
    
    copyfile(PPIfile,fullfile(PPIdir,PPIname));
end

% 7.5 define and estimate PPI GLM
cfg = [];
cfg.PPIname    = PPIname;
cfg.prefix     = 'swrf';
cfg.TR         = 1;
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    eval(sprintf('Subject%s',subjectID))
    
    cfg.run_nr{1}    = subjectdata.IMruns; % which runs are IM
    cfg.run_nr{2}    = subjectdata.UPCPruns; % which runs are UPCP
    cfg.identifier   = ['sub0' subjectID(2:end)]; %'945416';
    cfg.nifti_dir    = fullfile(root,subjectID,'Niftis');
    cfg.logfiles{1}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*IM_' subjectID '*']);
    cfg.logfiles{2}  = str2fullfile(fullfile(dataDir,subjectID,'Behaviour'),['*UPCP_' subjectID '*']);
    cfg.outputDir    = fullfile(root,subjectID,folderName);
    
    if ~exist(fullfile(cfg.outputDir,'beta_0001.nii'),'file')
        GLMConditionsPPI(cfg);
    end
end

% 7.6 make contrasts
conName  = 'UP';
folderName = ['GLM_Conditions_PPI_' conName];
for sub = 1:nsubjects
    
    subjectID = subjects{sub};
    cfg = [];
    cfg.root      = root;
    cfg.subjectID = subjectID;
    cfg.model_dir = folderName;
    cfg.tcontrasts = {[0 0 0 0 0 0 0 0 1]}; 
    cfg.tnames     = {'PPI interaction'};
    
    MakeContrasts(cfg)
end


%% %%%%%%%%%% GROUP ANALYSES %%%%%%%%%% %%

%% 8. Behavioural results
cfg = [];
cfg.root     = root;
cfg.subjects = subjects;
cfg.dir      = 'Behaviour';
cfg.dataName = 'behaviourData';

GroupBehaviour(cfg);


% 8.1 show accuracy and visibility per stim
V = []; A = [];
for sub = 1:nsubjects
    load(fullfile(root,subjects{sub},'Behaviour','behaviourData'),'V_stim','A_stim');
    
    V = cat(4,V,V_stim);
    A = cat(3,A,A_stim); 
    
    clear V_stim A_stim;   
end

figure;
subplot(2,1,1);
bar(squeeze(mean(A,3)))

%% 9. Average maps
%subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09',...
%    'S11','S12','S13','S14','S15','S16','S17','S18'};

cfg = [];
cfg.root     = root;
cfg.subjects = subjects;
cfg.mask     = true;
cfg.contrast = 'unconsciousDecoding.nii.gz';
cfg.dir      = 'Within';

GroupMeanNifti(cfg)


%% 10. Statistics decoding
% add decoding path

% 10.1 create WITHIN permutation maps for each subject
for sub = 1:nsubjects
    cfg = [];
    cfg.subjectID = subjects{sub};
    cfg.nPerm = 25;
    cfg.root = root;
    cfg.data  = fullfile('Within','unconsciousDecoding');
    if ~exist(fullfile(root,subjects{sub},[cfg.data '_perm.mat']),'file')
        DecodingSearchlightPermutation(cfg)
    end
end

% 10.2 create CROSS permutation maps for each subject
for sub = 1:nsubjects
    cfg = [];
    cfg.subjectID = subjects{sub};
    cfg.dataDir      = 'Within'; % data
    cfg.condition{1} = 'unconscious';
    cfg.condition{2} = 'imagery';
    cfg.nPerm = 25;
    cfg.root = root;
    cfg.data  = fullfile('Cross',[cfg.condition{1} '_to_' cfg.condition{2}]);
    
    if ~exist(fullfile(root,subjects{sub},[cfg.data '_perm.mat']),'file')
        CrossDecodingSearchlightPermutation(cfg)
    end
end

% 10.3 bootstrap permutations to create null-distributions
cfg.root       = root;
cfg.subjects   = subjects;
cfg.contrast   = 'unconsciousDecoding';
cfg.dir        = 'Within';
cfg.nBootstrap = 10000;
cfg.pair       = 5;

DecodingSearchlightBootstrap(cfg)


% 10.4 FDR correction for multiple comparisons
cfg = [];
pair = 7;
cfg.contrast          = 'conscious_to_unconscious';
cfg.dir               = 'CrossRun';
cfg.inputfile         = fullfile(root,'GroupResults',cfg.dir,['fpVal_' int2str(pair) '_' cfg.contrast '_btstrp.nii']);
cfg.qvalue            = 0.05; % FDR threshold (default = 0.05)
cfg.mask              = fullfile(root,'GroupResults','fmask.nii');

[~,map] = read_nii(cfg.inputfile); [~,mask] = read_nii(cfg.mask);
[~,fdr_threshold] = fdr_bh(1-map(mask==1),cfg.qvalue);

% 10.5 seperate clusters
cfg.inputfile = cfg.inputfile; 
cfg.threshold = 1-fdr_threshold; % corrected p value
cfg.numVox    = 50;
bb_separate_clustersND(cfg)

% 10.6 get significant accuracy clusters
cfg = [];
contrast = 'unconsciousDecoding_animacy';
cfg.root = root;
cfg.dir  = 'Within';
cfg.map  = sprintf('GA_%s.nii.gz',contrast);
cfg.pair = 7;
cfg.clusters = 'rpval_unconsciousDecoding_btstrp_clusters.nii.gz';
GetAccuracyClusters(cfg)

%% 11. Psychophysiological Interaction

subjects(15) = []; VVIQ(15)=[];
nsubjects = length(subjects);

% 11.1 correlations per condition
seed = 'LOC';
ROIs = {'V1','ldlPFC','rdlPFC'};
conds = {'UP','CP','IM'};

rAll = nan(nsubjects,length(conds),length(ROIs));
for sub = 1:nsubjects
    
    seed_PPI = cell(length(conds),1);
    ROI_PPIs = cell(length(conds),length(ROIs));
    
    for c = 1:length(conds)
        % PPI seed
        seed_file = str2fullfile(fullfile(root,subjects{sub},'GLM_Conditions'),['PPI_*' seed '_' conds{c} '.mat']);
        seed_PPI{c} = load(seed_file);
        
        % get other ROIs
        for r = 1:length(ROIs)
            ROI_file = str2fullfile(fullfile(root,subjects{sub},'GLM_Conditions'),['PPI_*' ROIs{r} '_' conds{c} '.mat']);
            ROI_PPIs{c,r} = load(ROI_file);
        end
    end
    
    % plot
    plotting = true;
    if plotting
        for c = 1:length(conds)
            plot(seed_PPI{c}.PPI.ppi(:),ROI_PPIs{1,c}.PPI.ppi(:),'.'); hold on
            
        end
    end
    
    % calculate correlations
    for c = 1:length(conds)
        for r = 1:length(ROIs)
            rAll(sub,c,r) = corr(seed_PPI{c}.PPI.ppi(:),ROI_PPIs{c,r}.PPI.ppi(:));
        end
    end
    
end

% 11.2 plot the contrasts
contrasts = {'CP min UP','CP min IM','UP min IM'};
contrastVs = [-1 1 0; 0 1 -1; 1 0 -1];
colors = ['b','r','g']; pVals = zeros(length(ROIs),length(conds));
figure; count = 1;
for r = 1:length(ROIs)
    
    for c = 1:length(conds)
        
        subplot(length(ROIs),length(conds),count);
        count = count+1;
        
        % calculate contrast
        con = squeeze(rAll(:,:,r))*contrastVs(c,:)';
        [~,pVals(r,c)] = ttest(con,0);
        
        % plot boxplot
        boxplot(con,'colors',colors(c),'labels',sprintf('LOC to %s',ROIs{r})); hold on;
        
        % add random jitter to dots
        tmp = ones(length(squeeze(rAll(:,c,r))),1); 
        tmp = tmp+(rand(size(tmp))-0.5)*0.1;          
        
        % plot the dots
        scatter(tmp,con,40,'MarkerEdgeColor',colors(c),'MarkerFaceColor',colors(c),'MarkerFaceAlpha',0.5','MarkerEdgeAlpha',0.5);
        hold on; ylim([-0.55 0.55]);
        ylabel(contrasts(c))
        
        title(sprintf('p-value %.3f',pVals(r,c)));
        
    end
end
