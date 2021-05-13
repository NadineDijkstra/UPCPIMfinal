function DecodingSearchlight(cfg)
% function DecodingSearchlight(cfg)

% add decoding path
addpath('/vol/ccnlab1/naddij/Analyses/Decoding')

% set random generator for repeatability
rng(1,'twister')

outputDir = fullfile(fileparts(cfg.nifti_dir),cfg.outputDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Get the searchlight indices
slradius = cfg.radius;

% load grey matter mask and functional
[~,GM]     = read_nii('/vol/ccnlab1/naddij/Templates/rmni_icbm152_gm_tal_nlin_sym_09a.nii');
[~,func]   = read_nii(fullfile(cfg.betas,'beta_0001.nii'));
mask       = ~isnan(func) & GM > 0.1;

% infer the searchlight indices
[vind,~] = searchlightIndices(mask,slradius);

nSearchlights = length(vind);

%% Get the betas
% identify beta numbers
load(fullfile(cfg.betas,'SPM'));
betaIDs = find(strncmp(SPM.xX.name,['Sn(1) ' cfg.condition],12));

% load the betas
nTrials = length(betaIDs);
hdr     = read_nii(fullfile(cfg.betas,sprintf('beta_%04d.nii',betaIDs(1))));
betas   = zeros([nTrials,hdr.dim]);

for t = 1:nTrials
    if mod(t,5) == 0
        fprintf('Reading in betas %d out of %d \n',t,nTrials)
    end
    [~,betas(t,:,:,:)] = read_nii(fullfile(cfg.betas,sprintf('beta_%04d.nii',betaIDs(t))));
end

%% Balance the trials - for the few wrongly balances subs

% get the labels
load(fullfile(cfg.betas,'labels'))
labels = eval(cfg.condition);

% extra balancing in the case of some imagery
[~,subj] = fileparts(fileparts(cfg.betas));
if strcmp(cfg.condition,'imagery') && ...
        any(strcmp(subj,{'S01','S02','S03','S04','S05','S06','S07','S08'}))
    
    % get the trial matrix
    load(fullfile('/vol/ccnlab1/naddij/UPCPIM',subj,'Behaviour','trialMatrixIM.mat'))
    nRunTrls = length(trialMatrix)/4;
    trialMatrix(:,end+1) = [ones(nRunTrls,1); ones(nRunTrls,1)*2; ones(nRunTrls,1)*3; ones(nRunTrls,1)*4];
    
    
    % check out balancing
    idxTrials = cell(4,4,2);
    numTrials = nan(4,4,2);
    for i = 1:4 % stim 1
        for j = 1:4 % stim 2
            for c = 1:2 % cue
                if i ~= j                    
                    idxTrials{i,j,c} = find(trialMatrix(:,1)==i & ....
                        trialMatrix(:,2) == j & trialMatrix(:,3) == c);
                    
                    numTrials(i,j,c) = length(idxTrials{i,j,c,1});
                end
            end
        end
    end   
    
    % actually balance them
    nTrials = min(min(min(min(numTrials))));
    selectedTrials = [];
    
    for i = 1:numel(idxTrials)
        if ~isempty(idxTrials{i})
            
            ind = idxTrials{i};
            while numel(ind)>nTrials
                ind(randi(numel(ind))) = []; % randomly delete one trial
            end
            selectedTrials = cat(1,selectedTrials,ind);
            clear ind
        end
    end
    
    % select only those trials
    betas = betas(selectedTrials,:,:,:);
    labels = labels(selectedTrials);    
end

%% Balance the trials - for subs with missing trials
% get the trial matrix
if strcmp(cfg.condition,'imagery')
    load(fullfile('/vol/ccnlab1/naddij/UPCPIM',subj,'Behaviour','trialMatrixIM.mat'))
    for t = 1:length(trialMatrix)
        tmp(t) = trialMatrix(t,trialMatrix(t,3));
    end
    trialMatrix = tmp'; clear tmp
else
    load(fullfile('/vol/ccnlab1/naddij/UPCPIM',subj,'Behaviour','trialMatrixUPCP.mat'))
    if strcmp(cfg.condition,'conscious')
        trialMatrix = trialMatrix(trialMatrix(:,2)==2,1);
    else
        trialMatrix = trialMatrix(trialMatrix(:,2)==1,1);
    end
end

% assign RUN ID's
nRunTrls = length(trialMatrix)/4;
trialMatrix(:,2) = [ones(nRunTrls,1); ones(nRunTrls,1)*2; ones(nRunTrls,1)*3; ones(nRunTrls,1)*4];

% for the extra balancing subs - select those trials
if exist('selectedTrials','var')
    trialMatrix = trialMatrix(selectedTrials,:); 
end

% get rid of missing trials from trial matrix
while ~isempty(find(labels(:,1) ~= trialMatrix(1:length(labels),1), 1))
    idx = find(labels(:,1) ~= trialMatrix(1:length(labels),1));
    trialMatrix(idx(1),:) = [];
end
labels = trialMatrix;

%% Balance the trials - for everybody, order stim class per run

% balance per run
nRuns = 4;
idx = cell(nRuns,1); betas2 = [];
Y = []; run_idx = [];
for r = 1:nRuns
    ind = find(labels(:,2) == r); 
    idx{r} = balance_trials(labels(ind,1),'downsample');
    
    Y = [Y; labels(ind(cell2mat(idx{r})),1)];
    run_idx = [run_idx; labels(ind(cell2mat(idx{r})),2)];
    betas2 = [betas2; betas(ind(cell2mat(idx{r})),:,:,:)];
end

betas = betas2; clear betas2

% save the betas and the labels
save(fullfile(outputDir,[cfg.condition 'Decoding']),'Y','betas','run_idx','mask','vind','-v7.3')

%% Do decoding per searchlight
% decoding settings
cfgD.gamma = cfg.gamma;
ind        = find(mask);

% pairwise decoding
pairs      = [1,2;1,3;1,4;2,3;2,4;3,4];
accuracy   = cell(size(pairs,1),1);
for p = 1:size(pairs,1)
    
    accuracy{p} = zeros(hdr.dim);
    fprintf('Decoding %d versus %d \n',pairs(p,1),pairs(p,2));
    
    % select the two classes
    idx        = ismember(Y,pairs(p,:));
    r_idx      = run_idx(idx);
    y          = Y(idx); y = y == pairs(p,1);
    X          = betas(idx,:,:,:);
    nTrials    = length(y);
    
    % create folds - leave one run out
    folds = cell(nRuns,1); 
    for r = 1:nRuns
        folds{r} = find(r_idx==r);
    end
    nFolds = length(folds);
    
    % run over searchlights
    for s = 1:nSearchlights
        
        if s >= (nSearchlights/10) && mod(s,round((nSearchlights/10))) == 0
            fprintf('Progress: %d percent of searchlights \n',round((s/nSearchlights)*100))
        end
        
        % mask the betas
        x = X(:,vind{s});
        
        % decoding
        Xhat = zeros(nTrials,1);
        for f = 1:nFolds
            testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
            labels = y(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
            
            % train
            decoder = train_LDA(cfgD,labels,trainX');
            %decoder = fitcsvm(trainX,labels);
            
            % decode
            Xhat(testidx) = decode_LDA(cfgD,decoder,testX');
            %Xhat(testidx) = predict(decoder,testX);
        end
        
        % determine accuracy
        accuracy{p}(ind(s)) = mean(Xhat > 0' == y);
        clear x
        
    end
end

% recode into matrix
acc = zeros([hdr.dim size(pairs,1)]);
for p = 1:size(pairs,1)
    acc(:,:,:,p) = accuracy{p};
end
clear accuracy
acc(:,:,:,p+1) = squeeze(mean(acc,4)); % add average

% write results
hdr.dim = size(acc);
write_nii(hdr, acc, fullfile(outputDir,[cfg.condition 'Decoding.nii.gz']))

clear acc
