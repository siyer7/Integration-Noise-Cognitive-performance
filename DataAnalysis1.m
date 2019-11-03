% Clear Workspace
close all;
sca;
clear vars;

var= load('data/uncertaintyV1-subject18-1-EarlyQuit.mat');
data = var.data;
contrast= data.response.contrast;
variability= data.response.variance;

lowc_count=1;
highv_count=1;
baseline_count=1;
lowc_condition=[];
highv_condition=[];
baseline_condition=[];

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

sum_lowc=0;
sum_highv=0;
sum_baseline=0;
accuracy= data.response.accuracy;

% adding up the number of trials with the correct response in each
% condition
for i = lowc_condition
    sum_lowc = sum_lowc + accuracy(i);
end
for i = highv_condition
    sum_highv= sum_highv+ accuracy(i);
end
for i = lowc_condition
    sum_baseline= sum_baseline + accuracy(i);
end

% dividing sum of correct trials by total number of trials in each
% condition
ratio_lowc = sum_lowc./size(lowc_condition,2);
ratio_highv = sum_highv./size(highv_condition,2);
ratio_baseline = sum_baseline./size(baseline_condition,2);

y = [ratio_lowc ratio_highv ratio_baseline];
x={'low contrast' 'high variability' 'baseline'};
bar(y)
ylabel('Proportion correct') 
set(gca,'xticklabel',x)
set(gca,'YLim',[0 1])
