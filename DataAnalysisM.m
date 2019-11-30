%close all;
clear vars;
clc;
sca;

%% Load Data

data = LoadData(); % load data

data.cue = categorical(data.cue, [-1 0 1], {'L', 'N', 'R'}); % relabel cues

% determine conditions
contLo = data.contrast==min(data.contrast);
contHi = data.contrast==max(data.contrast);
varLo = data.variance==min(data.variance);
varHi = data.variance==max(data.variance);

% 1:baseline 2:low-c 3:hi-v
condition = categorical((contHi & varLo) + 2*(contLo & varLo) + 3*(contHi & varHi), ...
    [0 1 2 3], {'other', 'baseline', 'low-c', 'hi-v'});

data.condition = condition; % store conditions

%% Accuracy

% Within Subject
[accGroup, subjID, condID] = findgroups(data.subject, data.condition);
[accBySubjCond, accSC_CI] = splitapply(@MeanCI, data.accuracy, accGroup);

accTable = table(subjID, condID, accBySubjCond, accSC_CI);

numSubj = size(unique(subjID), 1);

% Plot Within Subject
figure(1);
sgtitle("Accuracy per Condition by Subject")

for i=1:numSubj
    subplot(round(numSubj/2), 2, i)
    
    x = accTable{(4*i-3):(4*i), 'condID'};
    y = accTable{(4*i-3):(4*i), 'accBySubjCond'};
    ci = accTable{(4*i-3):(4*i), 'accSC_CI'};
    subj = accTable.subjID(4*i);
    
    bar(x,y);
    hold on

    errorbar(x, y, ci(:,1), ci(:,2), 'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');
    title(sprintf("Subject %d", subj));
    ylim([0 1]);
end

% Across Subjects
[accGroupAll, condID] = findgroups(accTable.condID);
[accAllSubj, accAllSubj_CI] = splitapply(@MeanCI, accTable.accBySubjCond, accGroupAll);

% Plot Across Subject
figure(2);
hold on

title("Accuracy per Condition");

bar(x, accAllSubj)
errorbar(x, accAllSubj, accAllSubj_CI(:,1), accAllSubj_CI(:,2), 'o', 'MarkerSize', 1, 'LineWidth', 2, 'Color', 'black');


%f2a = gramm('x', data.condition, 'y', data.accuracy);
%f2a.facet_wrap(data.subject);
%f2a.stat_summary('type', 'bootci', 'geom', 'bar');
%f2a.stat_summary('type', 'bootci', 'geom', 'black_errorbar');

%f2a.set_title("Accuracy by Condition");
%f2a.set_names('x', 'condition', 'y', 'proportion of trials correct', 'column', 'Subject');
%f2a.axe_property('ylim', [0 1]);

%f2a.draw();

%% Psychometric Curves

% generate psychometric curves
[psyGroup, subjID, cueID] = findgroups(data.subject, data.cue);
[binMeans, psychResp, psychErr] = splitapply(@psychometric, data.orientMean, data.responseR, psyGroup);

% get confidence intervals
errMin = psychResp - psychErr;%./2;
errMax = psychResp + psychErr;%./2;

pct = table(subjID, cueID, binMeans, psychResp, errMin, errMax); % psychometric curve table

% plot
figure

f2b = gramm('x', pct.binMeans, 'y', pct.psychResp, 'color', pct.cueID, 'ymin', errMin, 'ymax', errMax);
f2b.facet_wrap(pct.subjID, 'ncols', 2);
f2b.geom_point();

f2b.geom_line();
%f2b.stat_smooth();
f2b.geom_interval('geom', 'black_errorbar');

f2b.geom_vline('xintercept', 0, 'style', ':');
f2b.geom_hline('yintercept', 0.5, 'style', ':');

f2b.set_title("Psychometric Curves per Cue by Subject");
f2b.set_names('x', 'mean orientation', 'y', 'proportion of response CW',...
    'column', 'Subject', 'color', 'Cues');
f2b.draw();

%% Reported Confidence
figure

f3a = gramm('x', data.condition, 'y', data.confidence);
f3a.facet_wrap(data.subject);
f3a.stat_summary('type', 'bootci', 'geom', 'bar');
f3a.stat_summary('type', 'bootci', 'geom', 'black_errorbar');

f3a.set_title("Reported Confidence by Condition");
f3a.set_names('x', 'condition', 'y', 'reported confidence', 'column', 'Subject');

f3a.axe_property('ylim', [-.7 .9])

f3a.draw();

%% Overconfidence

% calculate overconfidence
oc = @(x, y) mean(x) - mean(y);

[ocGroup, subjID, condID] = findgroups(data.subject, data.condition);
overConf = splitapply(oc, data.accuracy, data.confidence, ocGroup);

oct = table(overConf, subjID, condID); % create table

% plot
figure

f3c = gramm('x', oct.condID, 'y', oct.overConf);
f3c.facet_wrap(oct.subjID);

f3c.stat_summary('type', 'bootci', 'geom', 'bar');

f3c.set_title("Over-Confidence by Condition");
f3c.set_names('x', 'condition', 'y', 'overconfidence', 'column', 'Subject');

f3c.draw();

%% Confidence and Reaction Time
figure

g1a = gramm('y', data.reactTime, 'x', data.confidence);
g1a.facet_grid(data.accuracy, data.subject);

g1a.stat_summary('type', 'bootci', 'geom', 'bar');
g1a.stat_summary('type', 'bootci', 'geom', 'black_errorbar');
g1a.stat_fit('fun', @(m,b,x)m*x+b, 'StartPoint', [0 0], 'geom', 'line');

g1a.axe_property('xlim', [-1.5 1.5], 'ylim', [0 1.5]);
g1a.set_title("Reaction Time by Confidence Level");
g1a.set_names('x', 'confidence judgement', 'y', 'reaction time', 'row', 'Accuracy', 'column', 'Subject');

g1a.draw();

%% Data Analysis Functions

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
            
            if subjNum==1 % subject01 for testing--not real data
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