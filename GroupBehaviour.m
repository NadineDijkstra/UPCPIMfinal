function GroupBehaviour(cfg)

outputDir = fullfile(cfg.root,'GroupResults',cfg.dir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end
addpath(genpath('/vol/ccnlab1/naddij/Analyses/RainCloudPlots-master'))

%% get the data for all subjects
nsubjects = length(cfg.subjects);
nCond     = 3;
acc       = zeros(nsubjects,nCond);
Dp        = zeros(nsubjects,nCond);
vis       = zeros(nsubjects,nCond,4);
MVis      = zeros(nsubjects,nCond);

for sub = 1:nsubjects
    
    load(fullfile(cfg.root,cfg.subjects{sub},cfg.dir,[cfg.dataName '.mat']),...
        'Vr','Ar','D','mVis')
    
    acc(sub,:) = Ar; vis(sub,:,:) = Vr; Dp(sub,:) = D;
    MVis(sub,:) = mVis;
    clear D Vr Ar mVis
    
end

%% plot the results
conditions = {'Unconscious','Conscious','Imagery'};

subplot(3,3,1:3); % d prime
colors = ['b','r','g'];
boxplot(Dp,'positions',[1 2 3],'labels',conditions,'Widths',0.2,'colors',colors);
hold on
for c = 1:nCond

    % add random jitter to dots
    tmp = repmat(c,1,length(Dp(:,c)))';
    tmp = tmp+(rand(size(tmp))-0.5)*0.05;
    
    % plot the dots
    scatter(tmp,Dp(:,c),40,'MarkerEdgeColor',colors(c),'MarkerFaceColor',colors(c),'MarkerFaceAlpha',0.5','MarkerEdgeAlpha',0.5);
    hold on  
end
hold on; plot(xlim,[0 0],'k--'); ylabel('D prime');

subplot(3,3,4:6); % accuracy
colors = ['b','r','g'];
boxplot(acc,'positions',[1 2 3],'labels',conditions,'Widths',0.2,'colors',colors);
hold on
for c = 1:nCond

    % add random jitter to dots
    tmp = repmat(c,1,length(acc(:,c)))';
    tmp = tmp+(rand(size(tmp))-0.5)*0.05;
    
    % plot the dots
    scatter(tmp,acc(:,c),40,'MarkerEdgeColor',colors(c),'MarkerFaceColor',colors(c),'MarkerFaceAlpha',0.5','MarkerEdgeAlpha',0.5);
    hold on  
end
hold on; plot(xlim,[0.5 0.5],'k--'); ylabel('Accuracy')

% visibility ratings
for c = 1:nCond
    subplot(3,3,6+c)
    colors = ['b','r','g'];
    boxplot(squeeze(vis(:,c,:)),'positions',[1 2 3 4],'labels',{'1','2','3','4'},'Widths',0.2,'colors',colors(c));
    hold on
    for r = 1:4
        
        % add random jitter to dots
        tmp = repmat(r,1,length(squeeze(vis(:,c,r))))';
        tmp = tmp+(rand(size(tmp))-0.5)*0.05;
        
        % plot the dots
        scatter(tmp,squeeze(vis(:,c,r)),40,'MarkerEdgeColor',colors(c),'MarkerFaceColor',colors(c),'MarkerFaceAlpha',0.5','MarkerEdgeAlpha',0.5);
        hold on
    end
    hold on; ylabel('Count'); xlabel('Visibility')
end

%% Raincloud plots
X{1} = Dp(:,1); X{2} = Dp(:,2); X{3} = Dp(:,3);
figure;
colours = [1 0 0; 0 1 0; 0 0 1];
rm_raincloud(X,colours,1)

% d prime
close; figure; lims = [-1 7];
for con = 1:3
    c = zeros(1,3);
    c(con) = 0.5;
    
    subplot(3,1,con);  
    h = raincloud_plot(Dp(:,con),'color',c,'box_on',1);
    h{1,2}.SizeData = 50;
    h{1,2}.YData = h{1,2}.YData -1;
    %h{1,3}.Position([1,2]) = h{1,3}.Position([1,2]) + 0.1
    
    xlim(lims); ylim([-2.5 2.5])
    title(conditions{con})
    
end


% visibility
figure; lims = [-0.1 1.1];
counter = 1;
for con = 1:3
    
    visibility = squeeze(vis(:,con,:))*[1:4]';
%     V = [];
%     if con < 3
%         tmp = round(squeeze(vis(:,con,:))*184);
%     else
%         tmp = round(squeeze(vis(:,con,:))*144);
%     end
% 
%     for sub = 1:nsubjects
%         for v = 1:4
%             V = [V; ones(tmp(sub,v),1)*v];
%         end
%     end
%     
    
    c = zeros(1,3);
    c(con) = 0.5;
    for v = 1:4
        
        subplot(3,4,counter);
        h = raincloud_plot(squeeze(vis(:,con,v)),'color',c,'box_on',1);
        h{1,2}.SizeData = 50;
        %h{1,2}.YData = h{1,2}.YData -1;
        %h{1,3}.Position([1,2]) = h{1,3}.Position([1,2]) + 0.1
        
        xlim(lims); %ylim([-2.5 2.5])
        %title(conditions{con})
        counter = counter + 1;
    end
    
end

