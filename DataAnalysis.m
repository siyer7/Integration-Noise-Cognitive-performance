var= load('data/uncertaintyV1-subject18-1-EarlyQuit.mat');
contrast= var.data.stimuli.contrastVal;
variability= var.data.stimuli.variabilityVal;
data = var.data;
count1=0;
count2=0;
count3=0;
lowc_condition= NaN(1, data.exp.numTrials); %stores trials/indices of this conditions
highv_condition= NaN(1, data.exp.numTrials);
baseline_condition= NaN(1, data.exp.numTrials);

for i = 0:863
    if contrast(i)==.15 % low contrast
        lowc_condition(count1)=i;
        count1 = count1+1;
    elseif variability(i)==10 % high variability
        highv_condition(count2) = i;
        count2 = count2+1;
    else %what other conditions besides these 3 do we include?
        baseline_condition(count3) = i;
        count3 = count3+1;
    end
end

sum_lowc=0;
sum_highv=0;
sum_baseline=0;
accuracy= var.data.response.accuracy;
for i = lowc_condition
    sum_lowc = sum_lowc + accuracy(i)
end
for i = highv_condition
    sum_highv= sum_highv+ accuracy(i)
end
for i = lowc_condition
    sum_baseline= sum_baseline + accuracy(i)
end

ratio_lowc = sum_lowc/size(lowc_condition);
ratio_highv = sum_highv/size(highv_condition);
ratio_baseline = sum_baseline/size(baseline_condition);

graph = [ratio_lowc ratio_highv ratio_basline];
bar(graph)