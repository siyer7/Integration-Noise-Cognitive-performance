% Clear Workspace
close all;
sca;
clear vars;
clc;

var= load('data/uncertaintyV1-subject18-1-EarlyQuit.mat');
data = var.data;
contrast= data.response.contrast;
variability= data.response.variance;

bucket1=[]; buck1len=1;
bucket2=[]; buck2len=1;
bucket3=[]; buck3len=1;
bucket4=[]; buck4len=1;
bucket5=[]; buck5len=1;
bucket6=[]; buck6len=1;

% segregating trials by angle
for i = 1:863
    angle= data.response.orientationMean(i);
    if angle <= -14
        bucket1(buck1len)=i;
        buck1len = buck1len+1;
    elseif angle > -14 && angle <= -7
        bucket2(buck2len)=i;
        buck2len = buck2len+1;
    elseif angle > -7 && angle <= 0
        bucket3(buck3len)=i;
        buck3len = buck3len+1;
    elseif angle >= 0 && angle <= 7
        bucket4(buck4len)=i;
        buck4len = buck4len+1;
    elseif angle >= 7 && angle <= 14
        bucket5(buck5len)=i;
        buck5len = buck5len+1;
    elseif angle > 14
        bucket6(buck6len)=i;
        buck6len = buck6len+1;
    end
end

L_buck1_sum = 0; L_buck2_sum = 0; L_buck3_sum = 0; L_buck4_sum = 0; L_buck5_sum = 0; L_buck6_sum = 0;
N_buck1_sum = 0; N_buck2_sum = 0; N_buck3_sum = 0; N_buck4_sum = 0; N_buck5_sum = 0; N_buck6_sum = 0;
R_buck1_sum = 0; R_buck2_sum = 0; R_buck3_sum = 0; R_buck4_sum = 0; R_buck5_sum = 0; R_buck6_sum = 0;

% summing up number of times subject chooses CW, for each cued condition:
% L/N/R, within each bucket
for i = bucket1
    cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L          
            L_buck1_sum= L_buck1_sum + 1;
        elseif cueLet == 0 % N
            N_buck1_sum= N_buck1_sum + 1
        elseif cueLet == 1 % R
            R_buck1_sum= R_buck1_sum + 1
        end
    end
end
L_buck1_proportion = L_buck1_sum./size(bucket1,2);
N_buck1_proportion = N_buck1_sum./size(bucket1,2);
R_buck1_proportion = R_buck1_sum./size(bucket1,2);


for i = bucket2
	cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L  
            L_buck2_sum= L_buck2_sum + 1;
        elseif cueLet == 0 % N 
            N_buck2_sum= N_buck2_sum + 1
        elseif cueLet == 1 % R 
            R_buck2_sum= R_buck2_sum + 1
        end
    end
end
L_buck2_proportion = L_buck2_sum./size(bucket2,2);
N_buck2_proportion = N_buck2_sum./size(bucket2,2);
R_buck2_proportion = R_buck2_sum./size(bucket2,2);


for i = bucket3
	cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L 
            L_buck3_sum= L_buck3_sum + 1;
        elseif cueLet == 0 % N 
            N_buck3_sum= N_buck3_sum + 1
        elseif cueLet == 1 % R 
            R_buck3_sum= R_buck3_sum + 1
        end
    end
end
L_buck3_proportion = L_buck3_sum./size(bucket3,2);
N_buck3_proportion = N_buck3_sum./size(bucket3,2);
R_buck3_proportion = R_buck3_sum./size(bucket3,2);


for i = bucket4
    cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L 
            L_buck4_sum= L_buck4_sum + 1;
        elseif cueLet == 0 % N 
            N_buck4_sum= N_buck4_sum + 1
        elseif cueLet == 1 % R 
            R_buck4_sum= R_buck4_sum + 1
        end
    end
end
L_buck4_proportion = L_buck4_sum./size(bucket4,2);
N_buck4_proportion = N_buck4_sum./size(bucket4,2);
R_buck4_proportion = R_buck4_sum./size(bucket4,2);


for i = bucket5
	cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L
            L_buck5_sum= L_buck5_sum + 1;
        elseif cueLet == 0 % N 
            N_buck5_sum= N_buck5_sum + 1
        elseif cueLet == 1 % R 
            R_buck5_sum= R_buck5_sum + 1
        end
    end
end
L_buck5_proportion = L_buck5_sum./size(bucket5,2);
N_buck5_proportion = N_buck5_sum./size(bucket5,2);
R_buck5_proportion = R_buck5_sum./size(bucket5,2);


for i = bucket6
	cueLet= data.response.cue(i);
    if data.response.responseRight(i)==1
        if cueLet == -1 % L 
            L_buck6_sum= L_buck6_sum + 1;
        elseif cueLet == 0 % N 
            N_buck6_sum= N_buck6_sum + 1
        elseif cueLet == 1 % R 
            R_buck6_sum= R_buck6_sum + 1
        end
    end
end
L_buck6_proportion = L_buck6_sum./size(bucket6,2);
N_buck6_proportion = N_buck6_sum./size(bucket6,2);
R_buck6_proportion = R_buck6_sum./size(bucket6,2);

y1 = [L_buck1_proportion L_buck2_proportion L_buck3_proportion L_buck4_proportion L_buck5_proportion L_buck6_proportion]
y2 = [N_buck1_proportion N_buck2_proportion N_buck3_proportion N_buck4_proportion N_buck5_proportion N_buck6_proportion];
y3 = [R_buck1_proportion R_buck2_proportion R_buck3_proportion R_buck4_proportion R_buck5_proportion R_buck6_proportion];
x = [-14 -9 -3 3 9 14];
figure, 
err1 = .01*ones(size(y1));
err2 = .01*ones(size(y1));
err3 = .01*ones(size(y1));

%= sqrt((y*[size(bucket1,2) size(bucket2,2) size(bucket3,2) size(bucket4,2) size(bucket5,2) size(bucket6,2)]-y)/
%errorbar(x_smooth,y_smooth,err1)
scatter(x,y1)
hold on
scatter(x,y2)
hold on
scatter(x,y3)
hold on
%errorbar(x,y2,err2)
%errorbar(x,y3,err3)

x_smooth= linspace(-14,14,50)
y_smooth1 = spline(x, y1, x_smooth);
y_smooth2 = spline(x, y2, x_smooth);
y_smooth3 = spline(x, y3, x_smooth);

line(x_smooth,y_smooth1)
hold on
line(x_smooth,y_smooth2)
hold on
line(x_smooth,y_smooth3)

set(gca,'YLim',[0 1])
set(gca,'XLim',[-15 15])