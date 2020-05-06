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

%% Bias Index
%this chunk does a few preliminary stuff, then creates a relevant table,
%and sends each subject's data @subject_overconf
clc;

numSubj = size(unique(subjID), 1);
[subject_groups] = findgroups(data.subject);
% 2 = baseline; 3 = low-c; 4= high-v
conds = double(condition(:,1)); 

% matrix with all attributes
data.Cue = double(data.cue);
data.Cue(data.Cue==1) = -1;
data.Cue(data.Cue==2) = 0;
data.Cue(data.Cue==3) = 1; % L = -1; N = 0; R = 1
data.orientMean(data.orientMean > 0) = 1; % correct answer = CW
data.orientMean(data.orientMean < 0) = -1; % correct answer = CCW

combine = [conds, data.Cue, data.orientMean, data.responseR, subject_groups];

%retaining only lowc & highv conditions
combine(combine(:,1) < 2, :) = [];
%retaining only informative cued trials
combine(combine(:,2) == 0, :) = [];
subject_groups = combine(:,5);

% combine(:,1) = conds = 2 (bsline), 3(lowC), 4(highV)
% combine(:,1) = cue = -1 (L), 1(R)
% combine(:,1) = orientMean = -1 (CCW), 1(CW)
% combine(:,1) = responseR = 0 (no, i.e., CCW), 1 (yes, i.e., CW)
% combine(:,1) = subject_groups = 1:10

[bySubjs] = splitapply(@(x1){subject_overconf(x1)}, combine, subject_groups);

bInd_bsline = zeros(10,3);
bInd_lowC = zeros(10,3);
bInd_highV = zeros(10,3);
for i = 1:10
    byConds = cell2mat(bySubjs{i});
    [bInd_bsline(i,:)] = byConds(1,:);
    [bInd_lowC(i,:)] = byConds(2,:);
    [bInd_highV(i,:)] = byConds(3,:);
end

x = [1 2 3];
% subject color and jitter for plotting
subjCol = [abs(linspace(0,.7, numSubj))', abs(linspace(-.8,.2, numSubj))',...
    abs(linspace(1, .2, numSubj))'];  % plotting colors for each subject
subjJit = linspace(-0.15,0.15,numSubj); % plotting jitter for each subject

for i = 1:10
    y = [bInd_bsline(i,1), bInd_lowC(i,1), bInd_highV(i,1)];
    y_lb = [bInd_bsline(i,2), bInd_lowC(i,2), bInd_highV(i,2)];
    y_ub = [bInd_bsline(i,3), bInd_lowC(i,3), bInd_highV(i,3)];
    
    if any(i==1 | i==2 | i==1) % informed subjects
        plot(x+subjJit(i), y, 's:', 'Color', subjCol(i,:),...
        'MarkerSize', 4, 'MarkerFaceColor', subjCol(i,:),'LineWidth', 1);
    else % naive subjects
        plot(x+subjJit(i), y, '.:', 'Color', subjCol(i,:),'MarkerSize', 20, 'LineWidth', 1);
    end
    hold on
    errorbar(x+subjJit(i), y, y_lb, y_ub, '.',...
        'MarkerSize', 0.1, 'Color', subjCol(i,:), 'CapSize', 0);
    ylim([-2 5]);
    xlim([0,4]);
end

hold on
bsline = mean(bInd_bsline)
lowC = mean(bInd_lowC)
highV = mean(bInd_highV)
y = [bsline(1), lowC(1), highV(1)];
y_lb = [bsline(2), lowC(2), highV(2)];
y_ub = [bsline(3), lowC(3), highV(3)];

plot(x, y, '.-', 'Color', 'black','MarkerSize', 15);
hold on
errorbar(x, y, y_lb, y_ub,'.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'black','CapSize', 0);

ylim([-5 5]);
xlim([0 4]);
xticks(x)
xticklabels({'baseline','low contrast','high variability'})
xlabel('conditions');
ylabel('bias index');
title('Cue Reliance');
%% Data Analysis Functions

function [c_bySubj] = subject_overconf(combine)
% this function essentially takes each suject's data, and splits it by condition
    % combine = [conds, cue, meanOrient, response, subj]
    % split data by condition
    cond = combine(:,1);
    [condition_groups] = findgroups(cond); 
	[c_bySubj] = splitapply(@(x1){byCond(x1)}, combine, condition_groups);
 end

function [c_test] = byCond(combine)
% receives data by condition, by subject
    % combine = [conds, cue, meanOrient, response, subj]
    %split data by cue
    cue = combine(:,2);
    [cue_groups] = findgroups(cue);
    
    boot = 100;
    c_byCond = splitapply(@(x1){bootstrp(boot,@byCue,x1)}, combine, cue_groups);
    cue_CW = c_byCond{2};
	cue_CCW = c_byCond{1};
    
    c_byCond = cue_CCW; 
    c_sorted = sortrows(c_byCond);
    c_med = nanmedian(c_sorted);
    c_lb = c_sorted(round(.025*boot));
    c_ub = c_sorted(round(.975*boot));
    c_test = [c_med c_lb c_ub];
end

function [c] = byCue(combine)
% receives data by cue, by condition, by subject
    % combine = [conds, cue, meanOrient, response, subj]
    meanOrient = combine(:,3);
    responseR = combine(:,4);
    HR = sum(meanOrient == 1 & responseR == 1)/sum(meanOrient == 1);
    FAR = sum(meanOrient == -1 & responseR == 1)/sum(meanOrient == -1);
    if(HR==1)
        HR = .99;
    elseif(HR==0)
        HR = .01;
    end
    if(FAR==1)
        FAR = .99;
    elseif(FAR==0)
        FAR = .01;
    end
    c = .5 * (norminv(HR) + norminv(FAR));
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