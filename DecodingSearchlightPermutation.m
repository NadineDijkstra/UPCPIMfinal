function DecodingSearchlightPermutation(cfg)
% function DecodingSearchlightPermutation(cfg)
addpath('/vol/ccnlab1/naddij/Analyses/Decoding')

% load the data
load(fullfile(cfg.root,cfg.subjectID,cfg.data))

%% Do decoding per searchlight
% decoding settings
nRuns      = 4;
cfgD.gamma = 0.01;
ind        = find(mask);

% pairwise decoding
pairs      = [1,2;1,3;1,4;2,3;2,4;3,4];
accuracy   = cell(size(pairs,1),cfg.nPerm);
nSearchlights = length(vind);

for p = 1:size(pairs,1)    
    
    fprintf('Decoding %d versus %d \n',pairs(p,1),pairs(p,2));
    
    % select the two classes
    idx        = ismember(Y,pairs(p,:));
    r_idx      = run_idx(idx);
    L          = Y(idx); L = L == pairs(p,1);
    X          = betas(idx,:,:,:);
    nTrials    = length(L);
    
    % create folds - leave one run out
    folds = cell(nRuns,1); 
    for r = 1:nRuns
        folds{r} = find(r_idx==r);
    end
    nFolds = length(folds);    
    
    
    for per = 1:cfg.nPerm
        fprintf('\t Permutation %d out of %d \n',per,cfg.nPerm)
        
        % permute labels per run
        y = zeros(size(L));
        for r = 1:nRuns
            idx = find(r_idx==r);
            y(idx) = L(idx(randperm(length(idx))));            
        end        
               
        % run over searchlights
        accuracy{p,per} = zeros(size(mask));   
        for s = 1:nSearchlights
            
            if s >= (nSearchlights/10) && mod(s,(nSearchlights/10)) == 0
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
            accuracy{p,per}(ind(s)) = mean(Xhat > 0' == y);
        end
        clear x
    end
end

% save
save(fullfile(cfg.root,cfg.subjectID,[cfg.data '_perm']),'accuracy');

rmpath('/vol/ccnlab1/naddij/Analyses/Decoding')

