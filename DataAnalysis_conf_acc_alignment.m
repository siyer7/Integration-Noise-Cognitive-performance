%close all;
clear vars;
clc;
sca;
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
%retaining only lowc & highv conditions
combine(combine(:,1) < 3, :) = [];
subject_groups = combine(:,5);
conds = combine(:,1);

% this section is for the identity line plot

points = (splitapply(@(x1){subject_overconf(x1)}, combine, subject_groups));
points = cell2mat(points);

% colors = [[0.1 .1 .1]; [0.1 .1 .1]; [0.3 .3 .3]; [0.3 .3 .3];...
%           [0.5 .5 .5]; [0.5 .5 .5]; [0.7 .7 .7]; [0.7 .7 .7];...
%           [0.9 .9 .9]; [0.9 .9 .9]; [0.1 .3 .5]; [0.1 .3 .5];...
%           [0.6 .4 .2]; [0.6 .4 .2]; [0.3 .5 .7]; [0.3 .5 .7];...
%           [0.8 .6 .4]; [0.8 .6 .4]; [0.5 .7 .9]; [0.5 .7 .9]];
% points(:,(3:5)) = colors;
% x = points(:,1);
% y = points(:,2);

for i = 1:2:20
    plot(points(i:i+1,1),points(i:i+1,2),':x','LineWidth',1,'MarkerIndices',(2));
%     errorbar(points(i:i+1,1), points(i:i+1,2), points(i:i+1,5), points(i:i+1,6), points(i:i+1,3), points(i:i+1,4));
    text(.35, 0.25, '. = low contrast')
	text(.35, 0.23, 'x = high variability')
    hold on
end
xlabel('\Delta accuracy between high & medium confidence')
ylabel('\Delta accuracy between medium & low confidence')
xlim([-.05 .5])
ylim([-.05 .5])
title('')
refline(1,0)

% figure('Name', 'Confidence-Accuracy alignment')
% scatter(x, y);
% hold on
% errorbar(x, y, se_x, 'horizontal','LineStyle', 'none');
% errorbar(x, y, se_y, 'vertical', 'LineStyle', 'none');
% xlabel('Confidence-Accuracy correlation in low-c')
% ylabel('Confidence-Accuracy correlation in high-v')
% title('')
% refline(1,0)

%% Data Analysis Functions

function [combo] = subject_overconf(combine)
% this function essentially takes each suject's data, and splits it by condition
% extract data from the table 'combine'    
    cond = combine(:,1);
    acc = combine(:,2);
    conf = combine(:,3);
    orientMean = combine(:,4);
%split data by condition
    [condition_groups] = findgroups(cond); 
%form a table with attributes
    combine = [acc, conf, orientMean, condition_groups];
% 	combo = splitapply(@condition_overconf, combine, condition_groups);
	combo = splitapply(@(x1){bootstrp(100,@condition_overconf,x1)}, combine, condition_groups);
    lowc = combo{1};
    highv = combo{2};
    [m_lowc,CI_lowc] = MeanCI(lowc);
    [m_highv,CI_highv] = MeanCI(highv);
%     combo = [m_lowc;m_highv];
    lowc_x_intervals = CI_lowc(:,1);
    lowc_y_intervals = CI_lowc(:,2);
    new_lowc_x_intervals = [lowc_x_intervals(1,:), lowc_x_intervals(2,:)];
    new_lowc_y_intervals = [lowc_y_intervals(1,:), lowc_y_intervals(2,:)];
    lowc_intervals = [new_lowc_x_intervals, new_lowc_y_intervals];
    highv_x_intervals = CI_highv(:,1);
    highv_y_intervals = CI_highv(:,2);
    new_highv_x_intervals = [highv_x_intervals(1,:), highv_x_intervals(2,:)];
    new_highv_y_intervals = [highv_y_intervals(1,:), highv_y_intervals(2,:)];
    highv_intervals = [new_highv_x_intervals, new_highv_y_intervals];
    combo = [m_lowc , lowc_intervals ; m_highv, highv_intervals]
end

function [combo] = condition_overconf(combine)
% receives data by condition, by subject
% extract data from table 'combine'
    acc = combine(:,1);
    conf = combine(:,2);
    combo = [acc,conf];
    acc_lowConf = mean(combo(combo(:,2)==-1,1));
    acc_mediumConf = mean(combo(combo(:,2)==0,1));
    acc_highConf = mean(combo(combo(:,2)==1,1));
    x = acc_highConf - acc_mediumConf;
    y = acc_mediumConf - acc_lowConf;
    combo = [x,y];
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