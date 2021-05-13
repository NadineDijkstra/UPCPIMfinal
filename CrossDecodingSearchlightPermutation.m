function CrossDecodingSearchlightPermutation(cfg)
% function CrossDecodingSearchlight(cfg)

% add decoding path
addpath('/vol/ccnlab1/naddij/Analyses/Decoding')

% set random generator for repeatability
rng(1,'twister')

%% Get the data
Betas = cell(2,1); Y_labels = cell(2,1);
Run_idx = cell(2,1);
for d = 1:2
    % load
    load(fullfile(cfg.root,cfg.subjectID,cfg.dataDir,[cfg.condition{d} 'Decoding.mat']),...
        'betas','Y','run_idx','mask','vind');
    
    % put in cell
    Betas{d} = betas;
    Y_labels{d} = Y;
    Run_idx{d} = run_idx;
end

nRuns = length(unique(Run_idx{2}));

%% Do decoding per searchlight
% decoding settings
hdr = read_nii(fullfile(cfg.root,cfg.subjectID,cfg.dataDir,[cfg.condition{1} 'Decoding.nii.gz']));
cfgD.gamma = 0.01;
ind        = find(mask);

% pairwise decoding
pairs      = [1,2;1,3;1,4;2,3;2,4;3,4];
accuracy   = cell(size(pairs,1),cfg.nPerm);
nSearchlights = length(vind);

for p = 1:size(pairs,1)
    
    
    fprintf('Decoding %d versus %d \n',pairs(p,1),pairs(p,2));
    
    % select the two classes
    y = cell(2,1); X = cell(2,1); r_idx = cell(2,1);
    for d = 1:2
        idx        = ismember(Y_labels{d},pairs(p,:));
        y{d}       = Y_labels{d}(idx); y{d} = y{d} == pairs(p,1);
        X{d}       = Betas{d}(idx,:,:,:);
        r_idx{d}   = Run_idx{d}(idx);
        clear idx
    end
    
    for per = 1:cfg.nPerm
        accuracy{p,per} = zeros(hdr.dim(1:3));
        fprintf('\t Permutation %d out of %d \n',per,cfg.nPerm)
        
        % permuting test labels per run
        yPtest = zeros(size(y{2}));
        for r = 1:nRuns
            idx = find(r_idx{d}==r);
            yPtest(idx) = y{2}(idx(randperm(length(idx))));
        end
        
        % run over searchlights
        for s = 1:nSearchlights
            
            if s >= (nSearchlights/10) && mod(s,(nSearchlights/10)) == 0
                fprintf('Progress: %d percent of searchlights \n',round((s/nSearchlights)*100))
            end
            
            % mask the betas
            x = cell(2,1);
            for d = 1:2
                x{d} = X{d}(:,vind{s});
            end            
            
            % decoding
            if any(strcmp(cfg.condition,'imagery'))
                % train decoder
                decoder = train_LDA(cfgD,y{1},x{1}');
                
                % decode
                Xhat    = decode_LDA(cfgD,decoder,x{2}');
                
            else  % leave one run out decoding
                Xhat = zeros(1,length(y{2}));
                for r = 1:nRuns
                    testidx = r_idx{2} == r; % this run
                    trainidx = r_idx{1} ~= r; % other runs
                    
                    ytest = y{2}(testidx); xtest = x{2}(testidx,:);
                    ytrain = y{1}(trainidx); xtrain = x{1}(trainidx,:);
                    
                    % train decoder
                    decoder = train_LDA(cfgD,ytrain,xtrain');
                    
                    % decode
                    Xhat(testidx)    = decode_LDA(cfgD,decoder,xtest');
                end
            end
            
            % determine accuracy
            accuracy{p,per}(ind(s)) = mean(Xhat > 0 == yPtest');
            clear x
        end
        
    end
end

% save
outputDir = fullfile(cfg.root,cfg.subjectID,fileparts(cfg.data));
save(fullfile(outputDir,[cfg.condition{1} '_to_' cfg.condition{2} '_perm']),'accuracy');

rmpath('/vol/ccnlab1/naddij/Analyses/Decoding')
