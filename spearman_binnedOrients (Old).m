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
%this chunk does a few preliminary stuff, then creates a relevant table,
%and sends each subject's data @subject_overconf
clc;
setGlobalx(0);

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
%keeping only those trial data where the condition is not 'other'
combine(combine(:,1) == 1, :) = [];
subject_groups = combine(:,5);
conds = combine(:,1);

%figure ('Name', 'Confidence vs Accuracy')

% to plot by subject [REMOVE 'test' AS ARGUMENT FROM subject_overconf()]
setGlobalx(2); % set position of plot for 1st subject

combined_differences = splitapply(@(x1){subject_overconf(x1)}, combine, subject_groups);
combined_differences = cell2mat(combined_differences);
%combined_differences
[fi,xi] = ksdensity(combined_differences);
figure()
plot(xi,fi);
str = sprintf('PD for differences in spearman correlation coefficient of Conf v Acc\n between the low-c & high-v conditions');
title(str);
xlim([-2 2]);
% to plot across subjects [ADD 'test' AS ARGUMENT TO subject_overconf()]
% setGlobalx(100);
% subject_overconf(combine);

%legend('red = baseline', 'blue = low-c', 'green = high-v', 'Location', 'northwest','Position',[0.2 .8 0.1 0.1])
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
figure(5)
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

function [differences] = subject_overconf(combine)
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
    count = getGlobalx; % gets location of plot for each subject
%     if(getGlobalx ~= 100)
%         subplot(ceil(4), 3, count)
%     end
%for each condition in each subject, we call @condition_overconf
    all_btstrps = splitapply(@(x1){boot(x1)}, combine, condition_groups);
    differences = all_btstrps{2} - all_btstrps{3};
%     [fi,xi] = ksdensity(differences);
%     figure()
%     plot(xi,fi);
%     xlim([-1.5 1.5]);

% shifts plot location by 1 for next subject
    setGlobalx(count + 1); 
end

function [all_btstrps] = boot(combine)
    %test = condition_overconf(combine);
    all_btstrps = [];
	all_btstrps = [all_btstrps; bootstrp(100,@condition_overconf,combine)];
%     um = bootstrp(100,@condition_overconf,combine);
%     all_btstrps;
%     [fi,xi] = ksdensity(um);
%     figure()
%     plot(xi,fi);
%     xlim([-1.5 1.5]);
end

function [spearman] = condition_overconf(combine)
% receives data by condition, by subject
% extract data from table 'combine'
    acc = combine(:,1);
    conf = combine(:,2);
    orientMean = combine(:,3);
    cond = combine(1,4);
% assigns colors to the different conditions
    color = 'black';
    switch cond
        case 1 % baseline
            color = 'red';
        case 2 % low-c
            color = 'blue';
        case 3 % high-v
            color = 'green';
    end
% dividing the signal (mean orientation) into 6* equi-depth bins
    sz = size(orientMean,1);
    buckets = 6;
    cutoffs = sz/buckets;
    iteration = 0;
    edges = [];
    for i = 1:length(orientMean(:,1))
    	iteration = iteration + 1;
        remainder = floor(mod(iteration,cutoffs));
        if(remainder == 0)
            edges = [edges round(orientMean(i,1))];
        end
    end
    edges = unique(edges); % important
    [N,~,bins] = histcounts(orientMean,edges);
    bins = bins + 1; % important for accumarray to work
        
% N = how many orientMeans in each bin;
% bins = which bin each orientMean is located in, therefore its size is the size of the actual data
                
    % mean accuracy, mean confidence for each bin
    meanAcc = accumarray(bins(:),acc,[],@mean);
    meanConf = accumarray(bins(:),conf,[],@mean);
% calculations for Confidence Intervals
    sqrtTrials = sqrt(accumarray(bins(:),orientMean,[],@length));
% confidence
    stddev_C = accumarray(bins(:),conf,[],@std);
% std error of confidence
    err_C = stddev_C./sqrtTrials;
% accuracy
    stddev_A = accumarray(bins(:),acc,[],@std);
% std error of accuracy
    err_A = stddev_A./sqrtTrials;

%     scatter(meanAcc,meanConf,color,'filled')
%     hold on
%     errorbar(meanAcc, meanConf, err_C, color,'LineStyle','none');
%     errorbar(meanAcc, meanConf, err_A, color,'horizontal', 'LineStyle','none');
%     hold on
%         
%     Fit = polyfit(meanAcc,meanConf,1);
    [spearman,PVAL] = corr(meanAcc,meanConf,'Type','Spearman');
%     plot(meanAcc,polyval(Fit,meanAcc),color,'LineWidth',1);
%     xlim([.5 1])        
%     ylim([-1 1])        
%     xlabel('Mean Accuracy')
%     ylabel('Mean Confidence')
%     if(getGlobalx < 5) % subjects were the authors
%     	title('Author ')% + (getGlobalx-1))
% 	elseif(getGlobalx ~= 100)
%         title('Subject ')% + str(getGlobalx-1))
%     else
%     	title('Confidence vs Accuracy')
%     end
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