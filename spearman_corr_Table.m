%close all;
clear vars;
clc;
sca;
global sub;
global diff;
%% Load Data

data = LoadData(); % load data
data.cue = categorical(data.cue, [-1 0 1], {'L', 'N', 'R'}); % relabel cues

% determine conditionsx, 
contLo = data.contrast==min(data.contrast);
contHi = data.contrast==max(data.contrast);
varLo = data.variance==min(data.variance);
varHi = data.variance==max(data.variance);

% 1:baseline 2:low-c 3:hi-v
condition = categorical((contHi & varLo) + 2*(contLo & varLo) + 3*(contHi & varHi), ...
    [0 1 2 3], {'other', 'baseline', 'low-c', 'hi-v'});

data.condition = condition; % store conditions

conditions = categorical((contHi & varLo) + 2*(contLo & varLo) + 3*(contHi & varHi), ...
    [1 2 3], {'baseline', 'low-c', 'hi-v'});



%% Accuracy

% Within Subject
[accGroup, subjID, condID] = findgroups(data.subject, data.condition);
[accBySubjCond, accSC_CI] = splitapply(@MeanCI, data.accuracy, accGroup);

accTable = table(subjID, condID, accBySubjCond, accSC_CI);

numSubj = size(unique(subjID), 1);

% Plot Within Subject
figure(1);
sgtitle("Accuracy per Condition by Subject")

nRows = 3;

for i=1:numSubj
    subplot(ceil(numSubj/nRows), nRows, i)
    
    x = accTable{(4*i-3):(4*i), 'condID'};
    y = accTable{(4*i-3):(4*i), 'accBySubjCond'};
    ci = accTable{(4*i-3):(4*i), 'accSC_CI'};
    subj = accTable.subjID(4*i);
    
    hold on
    bar(x,y);
    errorbar(x, y, ci(:,1), ci(:,2), 'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');
    
    title(sprintf("Subject %d", subj));
    ylim([0.4 1]);
    xlim({'baseline', 'hi-v'});
    ylabel("proportion correct");
    xlabel("condition");
end

% Across Subjects
[accGroupAll, condID] = findgroups(accTable.condID);
[accAllSubj, accAllSubj_CI] = splitapply(@MeanCI, accTable.accBySubjCond, accGroupAll);

% Plot Across Subject
figure(2);
hold on

title("Accuracy per Condition");

bar(x, accAllSubj);
errorbar(x, accAllSubj, accAllSubj_CI(:,1), accAllSubj_CI(:,2), ...
    'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');
scatter(accTable.condID, accTable.accBySubjCond, 50, 'red', 'filled', 'jitter','on', 'jitterAmount',0.25);

ylim([0.5 1]);
xlim({'baseline', 'hi-v'});
ylabel("proportion correct");
xlabel("condition");

%% conf V acc(by subj, by cond)
%this chunk does a few preliminary stuff, then creates a relevant table,
%and sends each subject's data @subject_overconf
clc;

numSubj = size(unique(subjID), 1);

% divide data by subject
[subject_groups] = findgroups(data.subject);

%next, since condition is stored as a string, we convert it into double and store it
% 2 = baseline; 3 = low-c; 4= high-v
conds = double(condition(:,1)); 

% matrix with all attributes
combine = [conds, data.accuracy, data.confidence, data.orientMean, subject_groups];
%sort by signal (mean orientation) to facilitate binning
combine = sortrows(combine,4);
%retaining only baseline, lowc & highv conditions
combine(combine(:,1) < 2, :) = [];
subject_groups = combine(:,5);
conds = combine(:,1);

% this section is for the ksdensity plot
setGlobalx(0)
setGlobaly([])
% stores the Probability distribution for all subjects' behavior
% figure('Name','Confidence-Accuracy alignment')
% 'p_val for \Delta spearman corr under the null hypothesis'
% 'Subjects'
spearmans = splitapply(@(x1){subject_overconf(x1)}, combine, subject_groups);

for i = 1:10
    sp_sub = spearmans{i};
    %get the median and CIs
    for j = 1:3
        sp_sub_cond = sp_sub{j};
        med = nanmedian(sp_sub_cond);
        sorted_sp = sort(sp_sub_cond);
        lb_CI = sorted_sp(25);
        ub_CI = sorted_sp(975);
        if (i == 1 && j == 1) % first table entry
            sp_tbl = table(med,lb_CI,ub_CI);
        else
            sp_tbl = [sp_tbl; table(med,lb_CI,ub_CI)];
        end
    end
end

save('sp_tbl.mat', 'sp_tbl')

% combined_differences = cell2mat(combined_differences);

% plot PD
% [fi,xi] = ksdensity(combined_differences);
% 'Combined'
% PD = [fi;xi];
% indices = find(PD(2,:) > 0);
% p_val = sum(PD(1,indices))/sum(PD(1,:))
% plot(xi,fi,'LineWidth',2.5,'Color','black');
% ax = gca;
% ax.FontSize = 20;
% xlim([-.5 1]);
% ylim([0 10]);
% xlabel('Probability distribution (bootstrap:1000)')
% ylabel('density')
% title({'Difference in Conf-Acc correlation', 'between','low contrast & high variability'})
% get X at which density is max
[maxYVal, indexAtMaxY] = max(fi);
xValueAtMaxYVal = xi(indexAtMaxY(1));
%% Data Analysis Functions

function [all_btstrps] = subject_overconf(combine)
% this function essentially takes each suject's data, and splits it by condition
    cond = combine(:,1);
    acc = combine(:,2);
    conf = combine(:,3);
    orientMean = combine(:,4);

%split data by condition
    [condition_groups] = findgroups(cond); 
%form a table with attributes
    combine = [acc, conf, orientMean, condition_groups];
    
    %calling bootstrapping function
	all_btstrps = splitapply(@(x1){bootstrp(1000,@condition_overconf,x1)}, combine, condition_groups);
	
    % once we receive by-condition data, we compute PD differences between lowc and hi-v conditions
%     lowc = all_btstrps{2};
%     highv = all_btstrps{3};
% 	differences = lowc - highv;% - all_btstrps{3};
  
% %this section is for the ksdensity plots
%     [fi,xi] = ksdensity(differences);
%     PD = [fi;xi];
%     indices = find(PD(2,:) < 0);
%     p_val = sum(PD(1,indices))/sum(PD(1,:));
%     [maxYVal, indexAtMaxY] = max(fi);
%     xValueAtMaxYVal = xi(indexAtMaxY(1));
%     % adding peak of PD for each subject
%     setGlobaly([getGlobaly xValueAtMaxYVal]);
%     hold on
%     % globalx gets subject number so we can use appropriate color
%     setGlobalx(getGlobalx+1);
%     plot(xi,fi,':','LineWidth',1,'Color',[0+(getGlobalx)*.07  .8-(getGlobalx)*.06  1+(getGlobalx)*-.08]);   
end

function [spearman] = condition_overconf(combine)
% receives data by condition, by subject
% extract data from table 'combine'
    acc = combine(:,1);
    conf = combine(:,2);
    [spearman,p_val] = corr(acc,conf,'Type','Spearman');
end

% get mean and 95% confidence intervals
function [xMean, xCI] = MeanCI(x)
    N = size(x,1); % determine size of data
    
    xMean = mean(x); % calculate mean
    xSEM = std(x)/sqrt(N); % calculate standard error
    
    CI = tinv([0.025 0.975], N-1); % calculate confidence probability intervals
    xCI = bsxfun(@times, xSEM, CI(:))'; % calculate confidence intervals 
end

%% Aux Functions

% for stripping data to completed trials
function [n] = LastTrial(data)
    n=1;

    while n<=size(data.response.correct,2)
        if isnan(data.response.correct(n))
            n = n - 1;
            break
        end
        n = n + 1;
    end

    n = n-1;
end

% load data as a table
function [tbl] = LoadData()
    files = dir(fullfile(pwd, 'data')); % get directory name
    tbl = nan; % init table as nan
    
    % loop over files
    for i = 1:length(files)
        
        % regex to get filename
        match = cell2mat(regexp(files(i).name,'subject\d{2}-\d{1,2}', 'match'));
        
        % if data file
        if ~isempty(match)
            subjNum = str2num(match(8:9)); % get subject number
            seshNum = str2num(match(11:end)); % get session number
            
            % which subjects to exclude
            % subject01 for testing--not real data
            % subj 2-4: informed, less data
            if any(subjNum==[1]) 
                continue
            end
            
            % load file
            var = load(fullfile(pwd, 'data', files(i).name));
            data = var.data.response;
            
            n = LastTrial(var.data); % index for completed trials
            
            %if subjNum==2
            %    data.responseRight = abs(1-data.responseRight);
            %    data.accuracy = abs(1-data.accuracy);
            %end
            
            % for table input
            subject = repmat(subjNum,n,1);
            session = repmat(seshNum,n,1);
            
            if ~istable(tbl) % first table entry
                tbl = table(...
                    data.correct(1:n)', ...
                    data.responseRight(1:n)', ...
                    data.accuracy(1:n)',...
                    data.reactionTime(1:n)',...
                    data.confidence(1:n)',...
                    data.isCuedBlock(1:n)',...
                    data.cue(1:n)',...
                    data.orientationMean(1:n)',...
                    data.trueOrientaions(1:n)',...
                    data.contrast(1:n)',...
                    data.variance(1:n)',...
                    subject,...
                    session,...
                    'VariableNames',...
                    {'correct','responseR','accuracy',...
                    'reactTime','confidence','isCued',...
                    'cue','orientMean','allOrientations','contrast',...
                    'variance','subject','session'});
            else
                tbl_ = table(... % append to table
                    data.correct(1:n)', ...
                    data.responseRight(1:n)', ...
                    data.accuracy(1:n)',...
                    data.reactionTime(1:n)',...
                    data.confidence(1:n)',...
                    data.isCuedBlock(1:n)',...
                    data.cue(1:n)',...
                    data.orientationMean(1:n)',...
                    data.trueOrientaions(1:n)',...
                    data.contrast(1:n)',...
                    data.variance(1:n)',...
                    subject,...
                    session,...
                    'VariableNames',...
                    {'correct','responseR','accuracy',...
                    'reactTime','confidence','isCued',...
                    'cue','orientMean','allOrientations','contrast',...
                    'variance','subject','session'});
                
                tbl = [tbl;tbl_]; % append to big table
            end          
        end % end if match
    end % end loop over files
end

function setGlobalx(val)
    global sub
    sub = val;
end
function sub = getGlobalx
    global sub;
    sub = sub;
end
function setGlobaly(val)
    global diff
    diff = val;
end
function diff = getGlobaly
    global diff;
    diff = diff;
end