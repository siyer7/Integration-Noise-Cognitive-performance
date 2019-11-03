%close all;
clear vars;
%clc;
sca;

%% Load Data

var = load('data/uncertaintyV1-subject18-1-EarlyQuit.mat');
data = var.data;

trial=1;
while trial<=size(data.response.correct,2)
    if isnan(data.response.correct(trial))
         break
    end
    trial = trial + 1;
end

correct = data.response.correct(1:trial);
responseR = data.response.responseRight(1:trial);
accuracy = data.response.accuracy(1:trial);
reactionTime = data.response.reactionTime(1:trial);
confidence = data.response.confidence(1:trial);
isCued = data.response.isCuedBlock(1:trial);
cue = data.response.cue(1:trial);
orientMean = data.response.orientationMean(1:trial);
contrast = data.response.contrast(1:trial);
variance = data.response.variance(1:trial);

correctMean = round((sign(orientMean)+1)/2);
responseR = abs(responseR - 1); % if ...Subject24-1-EarlyQuit.mat

contrastList = data.stimuli.contrastVal;
varianceList = data.stimuli.variabilityVal;

indxCueL = find(cue==-1);
indxCueR = find(cue==0);
indxCueN = find(cue==1);

indxContLo = find(contrast==contrastList(1));
indxContHi = find(contrast==contrastList(end));

indxVarLo = find(variance==varianceList(1));
indxVarHi = find(variance==varianceList(end));

indxBase = intersect(indxContLo, indxVarLo);
indxHCLV = intersect(indxContHi, indxVarLo);
indxLCHV = intersect(indxContLo, indxVarHi);

%% Accuracy
figure

mean(correct==correctMean)

accBase = mean(correct(indxBase));
accHCLV = mean(correct(indxHCLV));
accLCHV = mean(correct(indxLCHV));

accBaseM = mean(correctMean(indxBase)==responseR(indxBase));
accHCLVM = mean(correctMean(indxHCLV)==responseR(indxHCLV));
accLCHVM = mean(correctMean(indxLCHV)==responseR(indxLCHV));

accByCond = [accBase accHCLV accLCHV; accBaseM accHCLVM accLCHVM];
accByCondLabel = [reordercats(categorical({'BL', 'HCLV', 'LCHV'})); reordercats(categorical({'BL2', 'HCLV2', 'LCHV2'}))];

bar(accByCondLabel, accByCond);
title("Accuracy by conditition (1:genDist; 2:orientMean)");
ylabel("Proportion of trials correct");
xlabel("Condition");

%% Psychometric Curves
figure
hold on

orientMeanArr = [-100 -10 -7 -5 -2 0 2 5 7 10 100];
orientMeanBuckets = zeros(3,size(orientMeanArr, 2)-1);

respMeanBuckets = zeros(3,size(orientMeanArr, 2)-1);
respStdBuckets = zeros(3,size(orientMeanArr, 2)-1);

indxCue = [indxCueL indxCueN indxCueR];

for c=1:3 % loop over cues
    for i=1:size(orientMeanArr, 2)-1
        indxBucket = find(orientMean>=orientMeanArr(i) & orientMean<orientMeanArr(i+1));
        indxBucketCue = intersect(indxBucket, indxCue(c));
        
        if ~indxBucketCue
            orientMeanBuckets(c,i) = mean(bucket);
            respMeanBuckets(c,i) = 0;
        end

        bucket = orientMean(indxBucket); % orientation mean values
        resp = responseR(indxBucket); % response right values

        orientMeanBuckets(c,i) = mean(bucket);
        respMeanBuckets(c,i) = mean(resp);        
    end
    
    plot(orientMeanBuckets(c,:), respMeanBuckets(c,:), 'LineWidth', 2);
    
end

%legend('L', 'N', 'R')