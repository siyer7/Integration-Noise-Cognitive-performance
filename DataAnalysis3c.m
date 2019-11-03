% Clear Workspace
close all;
sca;
clear vars;

var= load('data/uncertaintyV1-subject18-1-EarlyQuit.mat');
data = var.data;
contrast= data.response.contrast;
variability= data.response.variance;

lowc_count=1; lowc_condition=[];
highv_count=1; highv_condition=[];
baseline_count=1; baseline_condition=[];


%segregating trials into the 3 conditions
for i = 1:863
    if contrast(i)==.15 && variability(i)==0 % low contrast
        lowc_condition(lowc_count)=i;
        lowc_count = lowc_count+1;
    elseif variability(i)==10 && contrast(i)==.6 % high variability
        highv_condition(highv_count) = i;
        highv_count = highv_count+1;
    elseif contrast(i)==.6 && variability(i)==0 %what other conditions besides these 3 do we include?
        baseline_condition(baseline_count) = i;
        baseline_count = baseline_count+1;
    end
end

sum_confidence_lowc=0;
sum_confidence_highv=0;
sum_confidence_baseline=0;
accuracy= data.response.accuracy;
sum_accuracy_lowc=0;
sum_accuracy_highv=0;
sum_accuracy_baseline=0;

% adding up confidece in each trial for each condition
% Low, Medium & High confidence reports are chosen to represent .33,.66,.99, which is completely arbitrary since the paper does not measure confidence in a similar way 
for i = lowc_condition
    switch(data.response.confidence(i))
        case -1
            sum_confidence_lowc = sum_confidence_lowc + 30;
        case 0
            sum_confidence_lowc = sum_confidence_lowc + 60;
        case 1
            sum_confidence_lowc = sum_confidence_lowc + 90;
    end
    sum_accuracy_lowc = sum_accuracy_lowc + accuracy(i);
end
for i = highv_condition
    switch(data.response.confidence(i))
        case -1
            sum_confidence_highv= sum_confidence_highv + 30;
        case 0
            sum_confidence_highv = sum_confidence_highv + 60;
        case 1
            sum_confidence_highv = sum_confidence_highv + 90;
    end
    sum_accuracy_highv = sum_accuracy_highv + accuracy(i);
end
for i = baseline_condition
    switch(data.response.confidence(i))
        case -1
            sum_confidence_baseline = sum_confidence_baseline + 30;
        case 0
            sum_confidence_baseline = sum_confidence_baseline + 60;
        case 1
            sum_confidence_baseline = sum_confidence_baseline + 90;
    end
    sum_accuracy_baseline = sum_accuracy_baseline + accuracy(i);
end

overconfidence_lowc = sum_confidence_lowc - sum_accuracy_lowc;
overconfidence_highv = sum_confidence_highv - sum_accuracy_highv;
overconfidence_baseline = sum_confidence_baseline - sum_accuracy_baseline;

% normalizing summed confidence out of 1 for each condition
ratio_lowc = overconfidence_lowc./size(lowc_condition,2);
ratio_highv = overconfidence_highv./size(highv_condition,2);
ratio_baseline = overconfidence_baseline./size(baseline_condition,2);

y = [ratio_lowc ratio_highv ratio_baseline];
x={'low contrast' 'high variability' 'baseline'};
bar(y)
ylabel('Overconfidence') 
set(gca,'xticklabel',x)
set(gca,'YLim',[0 100])
