%close all;
clear vars;
clc;
sca;
global count;
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
numSubj = size(unique(subjID), 1);
[subject_groups] = findgroups(data.subject);
conds = double(condition(:,1));
% 1 = other; 2 = baseline; 3 = low-c; 4= high-v
combine = [conds, data.accuracy, data.confidence, data.orientMean, subject_groups];

figure (5)
setGlobalx(2); % set position of plot of 1st subject
splitapply(@subject_overconf, combine, subject_groups);
legend('','red = baseline', '', 'blue = low-c', '', 'green = high-v', 'Location', 'northwest','Position',[0.2 .8 0.1 0.1])
%% Psychometric Curves
% per subject
% generate psychometric points
[psyGroup, subjID, condID, cueID] = findgroups(data.subject, data.condition, data.cue);
[binMeans, psychResp, psychErr] = splitapply(@psychometric, data.orientMean, data.responseR, psyGroup);

% get confidence intervals
errMin = -psychErr;
errMax = psychErr;

% fit psychometric curve
[alpha, beta, sz] = splitapply(@psychometricFit, data.orientMean, data.responseR, psyGroup);

% create tabe
pct = table(subjID, cueID, condID, binMeans, psychResp, errMin, errMax, alpha, beta); % psychometric curve table

keys = {'L', 'N', 'R'};
vals = {'red', 'black', 'blue'};
clrs = containers.Map(keys, vals);

% plot
figure(3)
sgtitle("Psychometric Curves");

for s=1:numSubj
    for i=1:4 % conditions
        
        if i==1
            continue
        end
        
        subplot(3, 1, i-1);
        hold on

        p = [];

        cond = string(pct{12*s-12+3*i, 'condID'});

        for c=0:2 % cues
            cue = string(pct{12*s-12+3*i+c-2, 'cueID'});

            bin = pct{12*s-12+3*i+c-2, 'binMeans'};
            resp = pct{12*s-12+3*i+c-2, 'psychResp'};

            a = pct{12*s-12+3*i+c-2, 'alpha'};
            b = pct{12*s-12+3*i+c-2, 'beta'};

            x = linspace(-20, 20, 100);
            y = glmval([a;b], x, 'logit');

            p = [p plot(x,y, ':', 'Color', clrs(cue), 'LineWidth', 2)];
            %plot(bin, resp, 'o', 'MarkerSize', 5, 'Color', clrs(cue));
        end

        legend(p, {'L', 'N', 'R'}, 'Location', 'best');
        title(cond);
    end
end

% example subject
figure(4)
hold on

s = 5; % subject number
subj = pct{12*s-1, 'subjID'};
sgtitle(strcat("Psychometric Curves Example Subject - ", string(subj)));

for i=1:4 % cond
    if i==1
        continue
    end
        
    subplot(3,1,i-1);
    hold on
    
    cond = string(pct{12*s-12+3*i, 'condID'});
    
    p = [];
    
    for c=0:2 % cue
        cue = string(pct{12*s-12+3*i+c-2, 'cueID'});
        
        a = pct{12*s-12+3*i+c-2, 'alpha'};
        b = pct{12*s-12+3*i+c-2, 'beta'};
        
        bin = pct{12*s-12+3*i+c-2, 'binMeans'};
        resp = pct{12*s-12+3*i+c-2, 'psychResp'};

        x = linspace(-40, 40, 100);
        y = glmval([a;b], x, 'logit');
        
        p = [p plot(x, y, 'Color', clrs(cue), 'LineWidth', 2)];
        scatter(bin, resp, clrs(cue), 'filled');
    end
    
    legend(p, {'L', 'N', 'R'}, 'Location', 'best');
    title(cond);
end

%% Data Analysis Functions

function [] = subject_overconf(combine) % split data by subject; called for each subject
    
    cond = combine(:,1);
    acc = combine(:,2);
    conf = combine(:,3);
    orientMean = combine(:,4);
    
    [condition_groups] = findgroups(cond); % split data by condition
    combine = [acc, conf, orientMean, condition_groups];
    count = getGlobalx; % gets location of plot for each subject
    subplot(ceil(4), 3, count)
    splitapply(@condition_overconf, combine, condition_groups);
	setGlobalx(count + 1); % shifts plot location by 1 for next subject
end

function [] = condition_overconf(combine)% called for each condition (for each subject)
    
    acc = combine(:,1);
    conf = combine(:,2);
    orientMean = combine(:,3);
    cond = combine(1,4);
    
    switch cond
        case 2
            color = 'red';
        case 3
            color = 'blue';
        case 4
            color = 'green';
    end
    if cond > 1
        
        %orientMean =sort(orientMean);
        y = quantile(orientMean,6);
        [N,~,bins] = histcounts(orientMean,6);
        % N = how many orientMeans in each bin; edges = the cutoffs for each bin; bins = which bin each orientMean is located in, therefore its size is the size of the actual data
        %bins = unique(bins)
        
        meanAcc = accumarray(bins(:),acc,[],@mean);
        meanConf = accumarray(bins(:),conf,[],@mean);
        stddev = accumarray(bins(:),conf,[],@std);
        sqrtTrials = sqrt(accumarray(bins(:),orientMean,[],@length));
        
        err = stddev./sqrtTrials;
        err = transpose(err);
        
%         plot(meanAcc,orientMean);
%         hold on
%         plot(meanConf,orientMean);
        
        scatter(meanAcc,meanConf,color,'filled')
        errorbar(meanAcc, meanConf, err, color,'LineStyle','none');
        hold on
        
        xlim([0.5 1]);
        Fit = polyfit(meanAcc,meanConf,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
        plot(meanAcc,polyval(Fit,meanAcc),color,'LineWidth',1);
        xlim([.5 1])        
        ylim([-1 1])        
        xlabel('Mean Accuracy')
        ylabel('Mean Confidence')
        title('conf vs acc')
    end
end

% get mean and 95% confidence intervals
function [xMean, xCI] = MeanCI(x)
    N = size(x,1); % determine size of data
    
    xMean = mean(x); % calculate mean
    xSEM = std(x)/sqrt(N); % calculate standard error
    
    CI = tinv([0.025 0.975], N-1); % calculate confidence probability intervals
    xCI = bsxfun(@times, xSEM, CI(:))'; % calculate confidence intervals 
end

% fit psychometric curve
function [t1, t2, sz] = psychometricFit(x, y)
    theta = glmfit(x, y, 'binomial', 'link', 'logit');
    t1 = theta(1);
    t2 = theta(2);
    
    sz = size(x, 1);
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

% calculate psychometric curve
function [binMeans, psychResp, psychErr] = psychometric(orientMean, responses)  
    quants = quantile(orientMean, 5, 1); % determine quantiles for binning
    quants = [-100; quants; 100]; % captural entire range
    
    bins = discretize(orientMean, quants); % index by bin
    
    % init return variables
    psychResp = zeros(1,6);
    psychErr = zeros(1,6);
    binMeans = zeros(1,6);
    
        % loop over bin index
    for i=1:6
        tempResp = []; % temporary varariables
        tempBins = []; % for data collection
        
        for j=1:length(bins) % loop over bin indeces
            if bins(j)==i
                tempResp = [tempResp responses(j)];
                tempBins = [tempBins orientMean(j)];
            end
        end
        
        % store data
        psychResp(i) = mean(tempResp); % p(response=1)
        psychErr(i) = sqrt(mean(tempResp)*(1-mean(tempResp))/length(tempResp)); % std dev
        binMeans(i) = mean(tempBins); % bin mean
    end 
end
function setGlobalx(val)
    global count
    count = val;
end
function count = getGlobalx
    global count;
    count = count;
end