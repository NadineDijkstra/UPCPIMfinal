function Behaviour_analysis(cfg)
% function Behaviour_analysis(cfg)
%
% cfg.subjectID = subject name
% cfg.root      = root directory
% cfg.plot      = plot results or not (true/false)
% cfg.outputDir = name of output folder

% data paths
dataPath  = fullfile(cfg.dataDir,cfg.subjectID,'Behaviour');
outputDir = fullfile(cfg.root,cfg.subjectID,cfg.outputDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Get the data for UPCP
data   = str2fullfile(dataPath,'*_UPCP_*');
if ~iscell(data); tmp{1} = data; data = tmp; end
animacy = []; visibility = []; trialMatrix = [];
% load the data
for d = 1:length(data)
    load(data{d},'A','V','P');
    
    % if wrong buttons were used
    tmp = unique(A(:,1)); tmp(isnan(tmp)) = [];
    if tmp == 2 % shifted buttons to the right  
       fprintf('Wrong buttons were used, adjust them now. \n')
       A(A(:,1)==2,1) = 1; % shift them back
       A(isnan(A(:,1)),1) = 2;        
    end
    
    animacy = [animacy; A];
    visibility = [visibility; V];
    trialMatrix = [trialMatrix; P.trialMatrix];
    
    clear A B T P
end

% Reverse scales
nTrials = length(animacy); 
for t = 1:nTrials
    
    if animacy(t,3)
        animacy(t,1) = 3-animacy(t,1);
    end
   
    if visibility(t,3)
        visibility(t,1) = 5-visibility(t,1);
    end   
    
end

% Calculate animacy scores per SOA
SOA = trialMatrix(:,2);
nSOAs = length(unique(SOA));

% recode animacy 
At = trialMatrix(:,1) < 3;
Ab = animacy(:,1) == 1;

% get accuracies
Ar = zeros(nSOAs,1);
for s = 1:nSOAs    
    Ar(s) = sum(Ab(SOA==s)==At(SOA==s))/(sum(SOA==s));      
end

% calculate d prime
H = zeros(nSOAs,1); Fa = zeros(nSOAs,1); 
D = zeros(nSOAs+1,1); 
for s = 1:nSOAs    
    
    H(s) = sum(Ab(At == 1 & SOA == s))/length(Ab(At == 1 & SOA == s));
    if H(s) == 0; H(s) = 0.00001; elseif H(s) == 1; H(s) = 1-0.00001; end % can't handle 0 or 1
    
    Fa(s) = sum(Ab(At == 0 & SOA == s))/length(Ab(At == 0 & SOA == s));    
    if Fa(s) == 0; Fa(s) = 0.00001; elseif Fa(s) == 1; Fa(s) = 1-0.00001; end
    
    D(s) = dprime(H(s),Fa(s),sum(At==1 & SOA == s),sum(At==0 & SOA == s));    
end

% Calculate visbility scores per SOA
nVis = 4;

% get scores
Vr = zeros(nSOAs,nVis);
mVis = zeros(nSOAs,1);
for s = 1:nSOAs
    mVis(s) = nanmean(visibility(SOA==s,1),1);
    for v = 1:nVis
        Vr(s,v) = sum(visibility(SOA==s)==v)/sum(SOA==s);
    end
end

% calculate accuracy and visibility per stim
nStim  = length(unique(trialMatrix(:,1)));
A_stim = zeros(nStim,nSOAs+1);
V_stim = zeros(nStim,nSOAs+1,4);
for st = 1:nStim    
    for s = 1:nSOAs
        idx = (trialMatrix(:,1)==st) & (trialMatrix(:,2) == s);
        
        if st < 3 % for animate
        A_stim(st,s) = sum(animacy(idx,1)==1)/sum(idx);
        elseif st > 2
        A_stim(st,s) = sum(animacy(idx,1)==2)/sum(idx);
        end  
        
        for v = 1:4
            V_stim(st,s,v) = sum(visibility(idx,1)==v)/sum(idx);
        end
    end    
end

%% Get the data for IM
data   = str2fullfile(dataPath,'*_IM_*');
if ~iscell(data); tmp{1} = data; data = tmp; end
animacyIM = []; visibilityIM = []; trialMatrixIM = [];

% load the data
for d = 1:length(data)
    load(data{d},'A','V','P');
    
    % if wrong buttons were used
    tmp = unique(A(:,1)); tmp(isnan(tmp)) = [];
    if tmp == 2 % shifted buttons to the right
       fprintf('Wrong buttons were used, adjust them now. \n')
       A(A(:,1)==2,1) = 1; % shift them back
       A(isnan(A(:,1)),1) = 2;        
    end
    
    animacyIM = [animacyIM; A];
    visibilityIM = [visibilityIM; V];
    trialMatrixIM = [trialMatrixIM; P.trialMatrix];
    
    clear A B T P
end

% recode respones
nTrials = length(animacyIM);
imagined = zeros(nTrials,1);
for t = 1:nTrials    
    imagined(t) = trialMatrixIM(t,trialMatrixIM(t,3));
    
    if animacyIM(t,3)
        animacyIM(t,1) = 3-animacyIM(t,1);
    end
   
    if visibilityIM(t,3)
        visibilityIM(t,1) = 5-visibilityIM(t,1);
    end       
end

% recode animacy 
At = imagined < 3;
Ab = animacyIM(:,1) == 1;

% get accuracy
Ar(nSOAs+1) =  sum(Ab==At)/length(At);

% get D prime
H = sum(Ab(At == 1))/length(Ab(At == 1));
if H == 0; H = 0.00001; elseif H == 1; H = 1-0.00001; end % can't handle 0 or 1

Fa = sum(Ab(At == 0))/length(Ab(At == 0));
if Fa == 0; Fa = 0.00001; elseif Fa == 1; Fa = 1-0.00001; end

D(nSOAs+1) = dprime(H,Fa,sum(At==1),sum(At==0));


% get score visibility
for v = 1:nVis
    Vr(nSOAs+1,v) = sum(visibilityIM(:,1)==v)/nTrials;
end
mVis(nSOAs+1) = nanmean(visibilityIM(:,1),1); 

% calculate accuracy and visibility per stim
for st = 1:nStim
   
    idx = imagined==st;
    if st < 3
        A_stim(st,nSOAs+1) = sum(animacyIM(idx,1)==1)/sum(idx);
    elseif st > 2
        A_stim(st,nSOAs+1) = sum(animacyIM(idx,1)==2)/sum(idx);
    end
    
    for v = 1:4
        V_stim(st,nSOAs+1,v) = sum(visibilityIM(idx,1)==v)/sum(idx);
    end
    
end

%% Save everything
save(fullfile(outputDir,'behaviourData'),'V_stim','A_stim','Vr','D','Ar','mVis')

%% Plot the things
if cfg.plot
subplot(4,3,1:3); 
bar(D)
set(gca,'XTicklabel',{'0 ms','67 ms','Imagined'})
xlim([0 4]); 
xlabel('Stimulus Onset Asynchrony')
ylabel('D prime')

subplot(4,3,4:6);
bar(Vr); 
ylabel('Counts'); legend('Visbility 1','Visbility 2','Visibility 3','Visibility 4')
set(gca,'XTicklabel',{'0 ms','67 ms','Imagined'})

subplot(4,3,7:9);
bar(A_stim'); legend('Stim 1','Stim 2','Stim 3','Stim 4');
set(gca,'XTicklabel',{'0 ms','67 ms','Imagined'})
ylabel('Accuracy'); ylim([0 1])

name = {'0 ms','67 ms','Imagined'};
for s = 1:3
   subplot(4,3,s+9)
   bar(squeeze(V_stim(:,s,:))); 
   xlabel('Stimulus'); ylabel('Counts')
   ylim([0 1]); 
   title(name{s}); legend('Visbility 1','Visbility 2','Visibility 3','Visibility 4')
end
end


