function DecodingSearchlightBootstrap(cfg)
% function DecodingSearchlightBootstrap(cfg)
% 
% cfg.root       = root;
% cfg.subjects   = subjects;
% cfg.contrast   = 'consciousDecoding';
% cfg.dir        = 'WithinFinDARTEL';
% cfg.nBootstrap = 10000;
% cfg.pair = 7;
%

% get some settings
nsubjects     = length(cfg.subjects);
load(fullfile(cfg.root,cfg.subjects{1},cfg.dir,[cfg.contrast '_perm.mat']),'accuracy')
nPermutations = size(accuracy,2); 
dim           = size(accuracy{1});
nPairs        = size(accuracy,1); clear accuracy

% load all permutations of right pair in memory
permutations  = zeros([nsubjects,nPermutations,dim]);
for sub = 1:nsubjects
    
    fprintf('Loading permutations for subject %d out of %d \n',sub,nsubjects)
    load(fullfile(cfg.root,cfg.subjects{sub},cfg.dir,'unconsciousDecoding_perm_animacy.mat'),'accuracy');%[cfg.contrast '_perm.mat']),'accuracy')
    
    for p = 1:nPermutations
%         if ismember(cfg.pair,1:nPairs) % select that pair
            permutations(sub,p,:,:,:) = accuracy{p};%accuracy{cfg.pair,p};
%         elseif cfg.pair == nPairs+1
%             tmp = zeros([nPairs,dim]);
%             for pa = 1:nPairs
%                 tmp(pa,:,:,:) = accuracy{pa,p};
%             end
%             permutations(sub,p,:,:,:) = mean(tmp,1); % average over pairs
%         end
    end
    
    clear accuracy
    
end

% bootstrap
bAccuracy     = zeros([cfg.nBootstrap,dim]);
for b = 1:cfg.nBootstrap
    
    fprintf('Running bootstrap %d of %d \n',b,cfg.nBootstrap)
    
    pAccuracy = zeros([nsubjects,dim]);
    for sub = 1:nsubjects
        
        per = randi(nPermutations); % random permutation
        pAccuracy(sub,:,:,:) = squeeze(permutations(sub,per,:,:,:));
        
    end
    bAccuracy(b,:,:,:) = mean(pAccuracy,1); % average over subs
    clear pAccuracy
end
clear permutations
%save(fullfile(cfg.root,'GroupResults',cfg.dir,sprintf('%s_bootstrap_pair%d.mat',cfg.contrast,cfg.pair)),'bAccuracy','-v7.3');

% compare with empirical distribution
%[V,Y] = read_nii(fullfile(cfg.root,'GroupResults',cfg.dir,['GA_' cfg.contrast '.nii.gz']));
[V,Y] = read_nii(fullfile(cfg.root,'GroupResults',cfg.dir,'GA_unconsciousDecoding_animacy.nii.gz'));
[~,mask] = read_nii(fullfile(cfg.root,'GroupResults',cfg.dir,['mask_unconsciousDecoding_animacy.nii.gz']));
%mask = squeeze(mask(:,:,:,cfg.pair)); 
mAcc = Y;%squeeze(Y(:,:,:,cfg.pair)); clear Y

% plot the distributions
subplot(3,2,1:2)
histogram(bAccuracy(:,mask==nsubjects)); xlim([0.42 0.58])
title('Permuted distribution')
subplot(3,2,3:4)
histogram(mAcc(mask==nsubjects)); xlim([0.42 0.58])
title('Empirical distribution')

% calculate p-values
pVal = nan(dim); tVal = nan(dim);
for x = 1:dim(1)
    for y = 1:dim(2)
        for z = 1:dim(3)
            if mask(x,y,z) == nsubjects % only look at voxels with all subs
                pVal(x,y,z) = sum(squeeze(bAccuracy(:,x,y,z)) > mAcc(x,y,z))/cfg.nBootstrap;
                tVal(x,y,z) = tinv(1-pVal(x,y,z),nsubjects-1);
            else  
                pVal(x,y,z) = NaN;
                tVal(x,y,z) = NaN;
            end
        end
    end
end

subplot(3,2,5:6); histogram(pVal);
title('p-values (# perm > emp)')

% GM mask
[~,gm] = read_nii('/vol/ccnlab1/naddij/Templates/rmni_icbm152_gm_tal_nlin_sym_09a.nii');
rpVal = 1-pVal;
rpVal(gm < 0.3) = NaN;
rpVal(mask ~= nsubjects) = NaN;

pVal(gm < 0.3) = NaN;
pVal(mask ~= nsubjects) = NaN;

tVal(gm < 0.3) = NaN;
tVal(mask ~= nsubjects) = NaN;

% write results
V.dim = dim;
cfg.contrast = 'unconsciousDecoding_animacy';
write_nii(V,pVal,fullfile(cfg.root,'GroupResults',cfg.dir,'pval_unconsciousDecoding_animacy_btstrp.nii.gz')); %['pVal_' int2str(cfg.pair) '_' cfg.contrast '_btstrp.nii.gz']))
write_nii(V,tVal,fullfile(cfg.root,'GroupResults',cfg.dir,'tval_unconsciousDecoding_animacy_btstrp.nii.gz')); %['tVal_' int2str(cfg.pair) '_'  cfg.contrast '_btstrp.nii.gz']))
write_nii(V,rpVal,fullfile(cfg.root,'GroupResults',cfg.dir,'rpval_unconsciousDecoding_animacy_btstrp.nii.gz')); %['rpVal_' int2str(cfg.pair) '_'  cfg.contrast '_btstrp.nii.gz']))



