%close all;
clear vars;
clc;
sca;

%% Load Data

nsamp = 1000; % bootstrap samples

data = LoadData(); % load data

% determine conditions 
contLo = data.contrast==min(data.contrast);
contHi = data.contrast==max(data.contrast);
varLo = data.variance==min(data.variance);
varHi = data.variance==max(data.variance);

% 1:baseline 2:low-c 3:hi-v
condition = (contHi & varLo) + 2*(contLo & varLo) + 3*(contHi & varHi);

data.condition = condition; % store conditions

% condition labels for plotting
num2cond = containers.Map([0 1 2 3],...
    {'other', 'baseline', 'low-c', 'hi-v'});
conds = {'baseline', 'low contrast', 'high variance'};

% cue labels for plotting
num2cue = containers.Map([-1 0 1], {'L', 'N', 'R'});
cues = {'L', 'N', 'R'};
    
% useful subj vars
numSubj = size(unique(data.subject), 1); % number of subjects
infsub = [2 3 4]; % informed subject nums

% subject color and jitter for plotting
subjCol = [abs(linspace(0,.7, numSubj))', abs(linspace(-.8,.2, numSubj))',...
    abs(linspace(1, .2, numSubj))'];  % plotting colors for each subject
subjJit = linspace(-0.15,0.15,numSubj); % plotting jitter for each subject


%% Accuracy

[accGroup, subjID, condID] = findgroups(data.subject, data.condition);
[accBySubjCond] = splitapply(@(x){bootstrp(nsamp,@mean,x)},...
    data.accuracy, accGroup);
[accCIBySubjCond] = splitapply(@(x){bootci(nsamp,@mean,x)},...
    data.accuracy, accGroup);

acc = arrayfun(@(x)mean(cell2mat(x)), accBySubjCond);
accCI = reshape(cell2mat(accCIBySubjCond), 2, [])';
accTable = table(subjID, condID, acc, accCI);

% Plot Within Subject
figure;
hold on;

lgnd = []; % legend labels
plts = []; % corresponding plots (avoid doubles due to errorbars)

diff = [];

for i=1:numSubj
    x = accTable{(4*i-2):(4*i), 'condID'};
    y = accTable{(4*i-2):(4*i), 'acc'};
    ci = accTable{(4*i-2):(4*i), 'accCI'} - y;
    subj = accTable.subjID(4*i);
    
    diff = [diff y(2)-y(3)];
    
    if any(infsub==subj)
        p = plot(x+subjJit(i), y, 's:', 'Color', subjCol(i,:),...
        'MarkerSize', 10, 'MarkerFaceColor', subjCol(i,:),...
        'LineWidth', 1.5);
    else
        p = plot(x+subjJit(i), y, '.:', 'Color', subjCol(i,:),...
        'MarkerSize', 30, 'LineWidth', 1.5);
    end
    
    errorbar(x+subjJit(i), y, ci(:,1), ci(:,2), '.',...
        'MarkerSize', 0.1, 'Color', subjCol(i,:), 'CapSize', 0);
    
    plts = [plts p];
    lgnd = [lgnd sprintf("Subject %d", subj)];
end

yline(0.5, '--', 'Color', 'black');

% Across Subjects
[accGroupAll, condID] = findgroups(accTable.condID);
[accAllSubj] = splitapply(@(x){bootstrp(nsamp,@mean,x)},...
    accTable.acc, accGroupAll);
[accCIByAllSubj] = splitapply(@(x){bootci(nsamp,@mean,x)},...
    accTable.acc, accGroupAll);

accAll = arrayfun(@(x)mean(cell2mat(x)), accAllSubj);
accAllCI = reshape(cell2mat(accCIByAllSubj), 2, [])'- accAll;

p = plot(x, accAll(2:4), '.-', 'Color', 'black',...
        'MarkerSize', 40);
errorbar(x, accAll(2:4), accAllCI(2:4,1), accAllCI(2:4,2), ...
     '.-', 'MarkerSize', 40, 'LineWidth', 2, 'Color', 'black',...
        'CapSize', 0);
 
plts = [plts p];
lgnd = [lgnd "Subject Mean"];

ylim([0.45 1]);
xlim([0.5 3.5]);

ax = gca;
ax.TickDir = 'out';
ax.FontSize = 16;
ax.YTick = [0.5 0.75 1];
ax.XTick = [1 2 3];
ax.XTickLabel = conds;

%title("Accuracy");
ylabel("proportion correct");
xlabel("condition");
%legend(plts, lgnd, 'Location', 'best', 'NumColumns', 2);
%legend('boxoff');

% plot differences inline
axes('Position',[.7 .8 .25 .15])
hold on

for s=1:numSubj    
    bar(s, diff(s), 'FaceColor', subjCol(s,:));
end

ax = gca;
ax.XTick = [1:10];
ax.XTickLabels = [2,3,4,5,8,9,10,11,12,13];
ylabel("\Delta accuracy");xlabel("subject");

yline(0);
title("Low Contrast - High Variability");

hold off

%% Accuracy Significance Test
accSamp = arrayfun(@(x){cell2mat(x)}, accBySubjCond);
accTable({"accSamp"}) = accSamp;

BL_LC = [];
BL_HV = [];
HV_LC = [];

for s=1:numSubj
    
end

%% Confidence

[confGroup, subjID, condID] = findgroups(data.subject, data.condition);
[confBySubjCond] = splitapply(@(x){bootstrp(nsamp,@mean,x)},...
    data.confidence, confGroup);
[confCIBySubjCond] = splitapply(@(x){bootci(nsamp,@mean,x)},...
    data.confidence, confGroup);

conf = arrayfun(@(x)mean(cell2mat(x)), confBySubjCond);
confCI = reshape(cell2mat(confBySubjCond), 2, [])';
confTable = table(subjID, condID, conf, confCI);

% Plot Within Subject
figure;
hold on;

lgnd = []; % legend labels
plts = []; % corresponding plots (avoid doubles due to errorbars)

diff = [];

for i=1:numSubj
    x = confTable{(4*i-2):(4*i), 'condID'};
    y = confTable{(4*i-2):(4*i), 'conf'};
    ci = confTable{(4*i-2):(4*i), 'confCI'} - y;
    subj = confTable.subjID(4*i);
    
    diff = [diff y(2)-y(3)];
    
    if any(infsub==subj)
        p = plot(x+subjJit(i), y, 's:', 'Color', subjCol(i,:),...
        'MarkerSize', 10, 'MarkerFaceColor', subjCol(i,:),...
        'LineWidth', 1.5);
    else
        p = plot(x+subjJit(i), y, '.:', 'Color', subjCol(i,:),...
        'MarkerSize', 30, 'LineWidth', 1.5);
    end
    
    errorbar(x+subjJit(i), y, ci(:,1), ci(:,2), '.',...
        'MarkerSize', 0.1, 'Color', subjCol(i,:), 'CapSize', 0);
    
    plts = [plts p];
    lgnd = [lgnd sprintf("Subject %d", subj)];
end

yline(0.5, '--', 'Color', 'black');

% Across Subjects
[confGroupAll, condID] = findgroups(confTable.condID);
[confAllSubj] = splitapply(@(x){bootstrp(nsamp,@mean,x)},...
    accTable.acc, accGroupAll);
[confCIByAllSubj] = splitapply(@(x){bootci(nsamp,@mean,x)},...
    confTable.conf, confGroupAll);

confAll = arrayfun(@(x)mean(cell2mat(x)), confAllSubj);
confAllCI = reshape(cell2mat(confCIByAllSubj), 2, [])'- confAll;

p = plot(x, confAll(2:4), '.-', 'Color', 'black',...
        'MarkerSize', 40);
errorbar(x, confAll(2:4), confAllCI(2:4,1), confAllCI(2:4,2), ...
     '.-', 'MarkerSize', 40, 'LineWidth', 2, 'Color', 'black',...
        'CapSize', 0);
 
plts = [plts p];
lgnd = [lgnd "Subject Mean"];

ylim([0.45 1]);
xlim([0.5 3.5]);

ax = gca;
ax.TickDir = 'out';
ax.FontSize = 16;
ax.YTick = [0.5 0.75 1];
ax.XTick = [1 2 3];
ax.XTickLabel = conds;

%title("Accuracy");
ylabel("mean confidence");
xlabel("condition");

%% Psychometric Curves and Bias Index
[psyGroup, subjID, condID, cueID] = findgroups(...
    data.subject, data.condition, data.cue);
[psyBySubjCond] = splitapply(@(x,y){bootstrp(nsamp,@psychofit,x,y)},...
    data.orientMean, data.responseR, psyGroup);
[psyCIBySubjCond] = splitapply(@(x,y){bootci(nsamp,@psychofit,x,y)},...
    data.orientMean, data.responseR, psyGroup);

psyTable = table(subjID, condID, cueID, psyBySubjCond);
save('psy_tbl','psyTable');

%% Bias Index
load('psy_tbl.mat')
[biGroup, subjID, condID] = findgroups(psyTable.subjID, psyTable.condID);
[biBySubjCond] = splitapply(@(x,y){biasIndex(x,y)},...
    psyTable.psyBySubjCond, psyTable.cueID, biGroup);

bi = arrayfun(@(x)median(cell2mat(x)), biBySubjCond);
biTable = table(subjID, condID, bi);

% exclude the first 10 ('other' condition) rows
% BI_tbl = table('Size',30);
for i = 1:40
    if(mod(i,4) == 1)
        continue
    end
    %get the median and CIs
    BI_sub_cond = biBySubjCond{i};
    med = median(BI_sub_cond);
    sortedBI = sort(BI_sub_cond);
    lb_CI = sortedBI(25);
    ub_CI = sortedBI(975);
    if (i == 2) % first table entry
        BI_tbl = table(med,lb_CI,ub_CI);
    else
        BI_tbl = [BI_tbl; table(med,lb_CI,ub_CI)];
    end
end

save('BI_tbl.mat', 'BI_tbl')

% plot
figure
hold on

lgnd = [];
plts = [];

for i=1:numSubj
    x = biTable{(4*i-2):(4*i), 'condID'};
    y = biTable{(4*i-2):(4*i), 'bi'};
    %ci = biTable{(4*i-2):(4*i), 'accCI'} - y;
    subj = biTable.subjID(4*i);
    
    if any(infsub==subj)
        p = plot(x+subjJit(i), y, 's:', 'Color', subjCol(i,:),...
        'MarkerSize', 10, 'MarkerFaceColor', subjCol(i,:),...
        'LineWidth', 1.5);
    else
        p = plot(x+subjJit(i), y, '.:', 'Color', subjCol(i,:),...
        'MarkerSize', 30, 'LineWidth', 1.5);
    end
    
%     errorbar(x+subjJit(i), y, ci(:,1), ci(:,2), '.',...
%         'MarkerSize', 0.1, 'Color', subjCol(i,:), 'CapSize', 0);
    
    plts = [plts p];
    lgnd = [lgnd sprintf("Subject %d", subj)];
end

ax = gca;
%ax.YScale = 'log';
ax.TickDir = 'out';
ax.FontSize = 14;
%ax.YTick = [0.5 0.75 1];
ax.XTick = [1 2 3];
ax.XTickLabel = conds;

%ylim([-1, 4]);

title("Bias Index per Condition");
ylabel("bias index");
xlabel("condition");
legend(plts, lgnd, 'Location', 'best', 'NumColumns', 2);
legend('boxoff');
hold off

%% Confidence/Accuracy Alignment
[caGroup, subjID] = findgroups(data.subject);
[caBySubjCond] = splitapply(@(x,y,z){confAccAlign(x,y,z)},...
    data.condition, data.confidence, data.accuracy, caGroup);
[caCIBySubjCond] = splitapply(@(x){bootci(nsamp,@mean,x)},...
    data.accuracy, accGroup);

acc = arrayfun(@(x)mean(cell2mat(x)), accBySubjCond);
accCI = reshape(cell2mat(accCIBySubjCond), 2, [])';
accTable = table(subjID, condID, acc, accCI);


%% Analysis Functions

function [] = alignment(conf, acc)
    corr(a',b','Type','Spearman');
end

% calculate logistic psychometric fit
function [fit] = psychofit(om, rr)
    sz = size(om, 1);
    
    theta = glmfit(om, rr, 'binomial', 'link', 'logit');
    t1 = theta(1);
    t2 = theta(2);
    
    fit = [t1, t2, sz];
end

function [bi] = biasIndex(fits, cue)
    cl = fits{1};
    ch = fits{3};
    
    bi = (ch(:,1)) - (cl(:,1));
end

% calculate psychometric curve points
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