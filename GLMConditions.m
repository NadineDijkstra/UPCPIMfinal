function GLMConditions(cfg)
% function FirstLevelConditions(cfg)

% output dir
if ~exist(cfg.outputDir,'dir'); mkdir(cfg.outputDir); end

% get the niftis 
[allRuns,idx] = sort([cfg.run_nr{1} cfg.run_nr{2}]);
conRuns = ismember(allRuns,cfg.run_nr{2})+1; % which con is this run
niftis  = []; nr_scans = zeros(length(allRuns),1);
for r = 1:length(allRuns)    
    tmp = str2fullfile(cfg.nifti_dir,sprintf([cfg.prefix '*%s-%04d*'],cfg.identifier,allRuns(r)));
    niftis = [niftis, tmp];
    nr_scans(r) = length(tmp); clear tmp
end

% sort the logfiles
tmp = [cfg.logfiles{:}];
logfiles = tmp(idx)';

%% Get the regressors
IM = []; 
CP = []; 
UP = []; 

animacyTextTime    = []; animacyResponseTime    = []; animacyTextDur = [];
visibilityTextTime = []; visibilityResponseTime = []; visibilityTextDur = [];
stim1Time = []; stim2Time = []; 

for l = 1:length(logfiles)
   
    load(logfiles{l},'A','V','P','T')
    V(isnan(V(:,2)),2) = 3; % replace nans
    A(isnan(A(:,2)),2) = 3;
    
    totalRunTime = nr_scans(l)*cfg.TR;
    if l == 1
        prevTime = 0;
    else
        prevTime = sum(nr_scans(1:l-1))*cfg.TR;
    end
    
    if conRuns(l) == 1 % IM     
        
        % fill in the timings per stimulus category
        IM = [IM; T.realtimings(:,7)+prevTime];
        
        % fill in the other timings
        animacyTextTime     = [animacyTextTime; T.realtimings(:,8)+prevTime];
        animacyResponseTime = [animacyResponseTime; T.realtimings(:,8)+A(:,2)+prevTime];
        visibilityTextTime  = [visibilityTextTime; T.realtimings(:,10)+prevTime];
        visibilityResponseTime = [visibilityResponseTime;T.realtimings(:,10)+V(:,2)+prevTime]; 

        stim1Time           = [stim1Time; T.realtimings(:,2)+prevTime];
        stim2Time           = [stim2Time; T.realtimings(:,4)+prevTime];
        
        aIdx = ones(length(T.realtimings(:,8)),1);
        vIdx = ones(length(T.realtimings(:,8)),1);
        
    elseif conRuns(l) == 2 % UPCP
        
        ind                 = P.trialMatrix(:,2) == 2; % split conscious-unconscious
        
        % fill in the timings per stimulus category
        cIdx                = ind & T.realtimings(:,2) < totalRunTime;        
        CP                  = [CP; T.realtimings(cIdx,2)+prevTime];
        
        uIdx                = ~ind & T.realtimings(:,2) < totalRunTime;        
        UP                  = [UP; T.realtimings(uIdx,2)+prevTime];
            
        % fill in the other timings
        
        aIdx                = T.realtimings(:,5)+A(:,2) < totalRunTime;
        
        animacyTextTime     = [animacyTextTime; T.realtimings(aIdx,5)+prevTime];
        animacyResponseTime = [animacyResponseTime; T.realtimings(aIdx,5)+A(aIdx,2)+prevTime];
        
        vIdx                = T.realtimings(:,7)+V(:,2) < totalRunTime;
        
        visibilityTextTime  = [visibilityTextTime; T.realtimings(vIdx,7)+prevTime];
        visibilityResponseTime = [visibilityResponseTime;T.realtimings(vIdx,7)+V(vIdx,2)+prevTime]; 
        
        % how many deleted because number of scans was less than time      
        fprintf('Deleted %d conscious trials and %d unconscious trials from block %d \n'...
            ,length(cIdx)/2-sum(cIdx), length(uIdx)/2-sum(uIdx),l)

    end     
    
    animacyTextDur      = [animacyTextDur; A(aIdx,2)];    
    visibilityTextDur   = [visibilityTextDur; V(vIdx,2)];    
    
    clear A V P T ind
end

%% Model specification
if ~exist(fullfile(cfg.outputDir,'SPM.mat'),'file')
    
    specification{1}.spm.stats.fmri_spec.dir = {cfg.outputDir};
    specification{1}.spm.stats.fmri_spec.timing.units = 'secs';
    specification{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    specification{1}.spm.stats.fmri_spec.sess.scans = niftis';
    
    
    % unconscious regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(1).name = 'UP';
    specification{1}.spm.stats.fmri_spec.sess.cond(1).onset = UP;
    specification{1}.spm.stats.fmri_spec.sess.cond(1).duration = ones(length(UP),1)*0.5;
        
    % conscious regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(2).name = 'CP';
    specification{1}.spm.stats.fmri_spec.sess.cond(2).onset = CP;
    specification{1}.spm.stats.fmri_spec.sess.cond(2).duration = ones(length(CP),1)*0.5;
      
    % imagery regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(3).name = 'IMstim1';
    specification{1}.spm.stats.fmri_spec.sess.cond(3).onset = IM;
    specification{1}.spm.stats.fmri_spec.sess.cond(3).duration = ones(length(IM),1)*4;
    counter = 4;
    
    % animacy text regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'animacy text'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = animacyTextTime;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = animacyTextDur;  
    counter = counter + 1;
    
    % animacy response regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'animacy response'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = animacyResponseTime;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = zeros(length(animacyResponseTime),1); % model as spike  
    counter = counter + 1;
    
    % visibility text regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'visibility text'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = visibilityTextTime;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = visibilityTextDur;  
    counter = counter + 1;
    
    % visibility response regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'visibility response'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = visibilityResponseTime;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = zeros(length(visibilityResponseTime),1); % model as spike  
    counter = counter + 1;
    
    % stim1 regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'stim1'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = stim1Time;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = zeros(length(stim1Time),1)+0.5;  
    counter = counter + 1;
    
    % stim2 regressor
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).name = 'stim2'; 
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).onset = stim2Time;  
    specification{1}.spm.stats.fmri_spec.sess.cond(counter).duration = zeros(length(stim2Time),1)+0.5;  
    
    % add WM and CSF regressors   
    [~,wm_mask] = read_nii(str2fullfile(cfg.nifti_dir,'rc2*.nii'));
    [~,csf_mask] = read_nii(str2fullfile(cfg.nifti_dir,'rc3*.nii'));
    wm = zeros(length(niftis),1); csf = zeros(length(niftis),1);
    for n = 1:length(niftis)        
        if mod(n,10) == 0
            fprintf('Calculating wm and csf for scan %d out of %d \n',n,length(niftis))
        end
       [~,scan] = read_nii(niftis{n});
       wm(n)    = mean(scan(wm_mask>0.9));  
       csf(n)   = mean(scan(csf_mask>0.9));
       clear scan       
    end
    
    specification{1}.spm.stats.fmri_spec.sess.regress(1).name = 'wm';
    specification{1}.spm.stats.fmri_spec.sess.regress(1).val = wm;
    specification{1}.spm.stats.fmri_spec.sess.regress(2).name = 'csf';
    specification{1}.spm.stats.fmri_spec.sess.regress(2).val = csf;
    
    % movement regressor
    specification{1}.spm.stats.fmri_spec.sess.multi = {''}; 
    specification{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    specification{1}.spm.stats.fmri_spec.sess.multi_reg = {str2fullfile(cfg.nifti_dir,'/rp*.txt')}; 
    specification{1}.spm.stats.fmri_spec.sess.hpf = 128;
    
    specification{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    specification{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    specification{1}.spm.stats.fmri_spec.volt = 1;
    specification{1}.spm.stats.fmri_spec.global = 'None';
    specification{1}.spm.stats.fmri_spec.mask = {''};
    specification{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % run the model
    spm_jobman('run',specification) 
    clear niftis   
  
end

if ~exist(fullfile(cfg.outputDir,'beta_0001.nii'),'file')
    
    % concatenate the different sessions
    if ~exist(fullfile(cfg.outputDir,'SPM_backup.mat'),'file')
        fprintf('Concatenating scan sessions \n');
        spm_fmri_concatenate(fullfile(cfg.outputDir,'SPM.mat'),nr_scans');
    end
    
    % estimate the model
    spm_file = str2fullfile(cfg.outputDir,'SPM.mat');
    estimation{1}.spm.stats.fmri_est.spmmat = {spm_file};
    estimation{1}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('run',estimation)
    
    %fprintf('Running model estimation for subject %d',sub)
    %qsubfeval('runBatch',estimation,'timreq',60*60*4,'memreq',1024^3*8,'options','-V','batchid',sprintf('Estimation %s',subjectname))
    
end