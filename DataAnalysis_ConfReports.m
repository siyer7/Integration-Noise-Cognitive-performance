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

% useful vars
numSubj = size(unique(subjID), 1);
subjCol = [ones(numSubj,1)-0.8, linspace(.1,1, numSubj)',...
    linspace(.7, .9, numSubj)'];
infsub = [2 3 4]; % informed subjects


%% Confidence (explicit)
clc
% Within Subject
[confGroup, subjID, condID] = findgroups(data.subject, data.condition);
[confBySubjCond, confSC_CI] = splitapply(@MeanCI, data.confidence, confGroup);

confTable = table(subjID, condID, confBySubjCond, confSC_CI);

% Plot Within Subject
figure(1);
hold on;

lgnd = []; % legend labels
plts = []; % corresponding plots (avoid doubles due to errorbars)

for i=1:numSubj
    x = confTable{(4*i-2):(4*i), 'condID'};
    y = confTable{(4*i-2):(4*i), 'confBySubjCond'};
    ci = confTable{(4*i-2):(4*i), 'confSC_CI'};
    subj = confTable.subjID(4*i);
    
    if any(infsub==subj)
        p = plot(x, y, 's-', 'Color', subjCol(i,:),...
        'MarkerSize', 10, 'MarkerFaceColor', subjCol(i,:));
    else
        p = plot(x, y, '.-', 'Color', subjCol(i,:),...
        'MarkerSize', 30);
    end
    
    errorbar(x, y, ci(:,1), ci(:,2), '.-',...
        'MarkerSize', 0.1, 'LineWidth', 1, 'Color', subjCol(i,:));
    
    plts = [plts p];
    lgnd = [lgnd sprintf("Subject %d", subj)];
end

% Across Subjects
[confGroupAll, confID] = findgroups(confTable.condID);
[confAllSubj, confAllSubj_CI] = splitapply(@MeanCI, confTable.confBySubjCond, confGroupAll);

p = plot(x, confAllSubj(2:4), '.-', 'Color', 'black',...
        'MarkerSize', 40);
errorbar(x, confAllSubj(2:4), confAllSubj_CI(2:4,1), confAllSubj_CI(2:4,2), ...
     '.-', 'MarkerSize', 40, 'LineWidth', 2, 'Color', 'black');
 
plts = [plts p];
lgnd = [lgnd "Subject Mean"];

title("Confidence per Condition by Subject");
ylim([-1 1]);
xlim({'baseline', 'hi-v'});
ylabel("confidence reported");
xlabel("condition");
legend(plts, lgnd, 'Location', 'best', 'NumColumns', 2);
legend('boxoff');

% % Plot Across Subject
% figure(2);
% hold on
% 
% title("Accuracy per Condition");
% 
% bar(x, accAllSubj);
% errorbar(x, accAllSubj, accAllSubj_CI(:,1), accAllSubj_CI(:,2), ...
%     'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');
% scatter(accTable.condID, accTable.accBySubjCond, 50, [.3 .3 .3], 'filled', 'jitter','on', 'jitterAmount',0.25);
% 
% ylim([0.5 1]);
% xlim({'baseline', 'hi-v'});
% ylabel("proportion correct");
% xlabel("condition");

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

%% example subject
figure(4)
hold on

s = 4; % subject number
subj = pct{12*s-1, 'subjID'};
sgtitle(strcat("Psychometric Curves Example Subject - ", string(subj)));

for i=1:4 % cond
    if i==1
        continue
    end
        
    subplot(1,3,i-1);
    hold on
    
    cond = string(pct{12*s-12+3*i, 'condID'});
    
    p = [];
    
    for c=0:2 % cue
        cue = string(pct{12*s-12+3*i+c-2, 'cueID'});
        
        a = pct{12*s-12+3*i+c-2, 'alpha'};
        b = pct{12*s-12+3*i+c-2, 'beta'};
        
        bin = pct{12*s-12+3*i+c-2, 'binMeans'};
        resp = pct{12*s-12+3*i+c-2, 'psychResp'};

        x = linspace(-20, 20, 41);
        y = glmval([a;b], x, 'logit');

        %plot([round(a), round(a)],[0, y(x==round(a))], 'Color', clrs(cue), 'LineWidth', 2);
        p = [p plot(x, y, 'Color', clrs(cue), 'LineWidth', 2)];
        scatter(bin, resp, clrs(cue), 'filled');
        xlabel("orientation mean");
        ylabel("probability response CW");
    end
    
    legend(p, {'L', 'N', 'R'}, 'Location', 'best');
    title(cond);
end

%% Bias Index

% per subject
[biasGroup, subjID, condID] = findgroups(pct.subjID, pct.condID);
biasIndx = splitapply(@biasIndex, pct.alpha, pct.cueID, biasGroup);

biasTable = table(subjID, condID, biasIndx);

figure(5)
sgtitle("Bias Index per Subject");
nRows = 3;

for s=1:numSubj
    subplot(ceil(numSubj/nRows), nRows, s);
    
    subj = biasTable{s*4, 'subjID'};
    
    x = biasTable{s*4-2:s*4, 'condID'};
    y = biasTable{s*4-2:s*4, 'biasIndx'};
    
    bar(x,y)
    
    grid on
    
    xlim({'baseline', 'hi-v'});
    ylim([min(min(y)-0.5,0), max(y)+0.5]);
    
    ylabel("bias index");
    xlabel("condition");
    
    title(strcat("Subject ", string(subj)));
end

% across subject
[biasAllGroup, condID] = findgroups(biasTable.condID);
[biasIndxAll, errBiasAll] = splitapply(@MeanCI, biasTable.biasIndx, biasAllGroup);

biasAllTable = table(condID, biasIndxAll, errBiasAll);

figure(6)
hold on
title("Bias Index by Condition");

bar(biasAllTable.condID, biasAllTable.biasIndxAll);
errorbar(biasAllTable.condID, biasAllTable.biasIndxAll, errBiasAll(:,1), errBiasAll(:,2),...
    'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');
scatter(biasTable.condID, biasTable.biasIndx, 50, [.3 .3 .3], 'filled', 'jitter','on', 'jitterAmount',0.25);

xlim({'baseline', 'hi-v'});

xlabel("condition");
ylabel("bias index");

%% Overconfidence
[confGroup, subjID, condID, confID] = findgroups(data.subject, data.condition, data.confidence);
accByConf = splitapply(@mean, data.accuracy, confGroup);

confTable = table(subjID, condID, confID, accByConf)

[confAll, condID, confID] = findgroups(confTable.condID, confTable.confID);
[accByConfAll, accConfErr] = splitapply(@MeanCI, confTable.accByConf, confAll);

confAllTable = table(condID, confID, accByConfAll, accConfErr);

% plot
hold on

for i=1:4
    if i==1
        continue
    end
    
    y = confAllTable{3*i-2:3*i, 'accByConfAll'};
    x = confAllTable{3*i-2:3*i, 'confID'};
    err = confAllTable{3*i-2:3*i, 'accConfErr'};
    
    p = [p plot(x, y, 'o--', 'LineWidth', 4)];
    %errorbar(x, y, err(1), err(2), )
end

title("Mean Accuracy by Confidence")
xline(-1);xline(0);xline(1);
ylim([0.5, 1]);
xlim([-1.1, 1.1]);
xlabel("confidence rating");
ylabel("mean accuracy");
legend({'BL', 'lo-c', 'hi-v'}, 'Location', 'best');    
    
%% Data Analysis Functions

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

% calculate bias index
function [bi] = biasIndex(a, cue)
    bi = a(cue==categorical({'R'})) - a(cue==categorical({'L'}));
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
                    data.contrast(1:n)',...
                    data.variance(1:n)',...
                    subject,...
                    session,...
                    'VariableNames',...
                    {'correct','responseR','accuracy',...
                    'reactTime','confidence','isCued',...
                    'cue','orientMean','contrast',...
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
                    data.contrast(1:n)',...
                    data.variance(1:n)',...
                    subject,...
                    session,...
                    'VariableNames',...
                    {'correct','responseR','accuracy',...
                    'reactTime','confidence','isCued',...
                    'cue','orientMean','contrast',...
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