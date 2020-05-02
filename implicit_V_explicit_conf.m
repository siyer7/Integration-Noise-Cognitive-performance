clc;
load('sp_tbl.mat')
load('BI_tbl.mat')
load('conf_tbl')

table2array(conf_tbl);
arr = [ table2array(conf_tbl), table2array(sp_tbl)];
% BI    CI_lb   CI_ub   spearman    CI_lb   spearman_ub

figure('Name', 'Conf-Acc corr VS Confidence')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low-c
subplot(1, 3, 1)
x = zeros([10 1]);
y = zeros([10 1]);
index = 1;
for i = 2:3:30
    x(index) = arr(i,1);
    y(index) = arr(i,4);
    index = index + 1;
    scatter(arr(i,1),arr(i,4),60,'filled');
    hold on
    errorbar(arr(i,1),arr(i,4), arr(i,5),arr(i,6), arr(i,2),arr(i,3),'LineStyle','none','CapSize', 0);
end
corr(x,y);
vals = [x,y];
Fits = bootstrp(1000,@polyfits,vals);
med = median(Fits);
med(1)
p_val = sum(Fits(:,1) > 0)/1000
plot(x,polyval(med,x),'LineWidth',2);
title('low-c')
xlim([-.8 .8])
ylim([-.1 1])
xlabel('confidence reports')
ylabel('conf-acc correlation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% high-v
subplot(1, 3, 2)
x = zeros([10 1]);
y = zeros([10 1]);
index = 1;
for i = 3:3:30
    x(index) = arr(i,1);
    y(index) = arr(i,4);
    index = index + 1;
    scatter(arr(i,1),arr(i,4),60,'filled');
    hold on
    errorbar(arr(i,1),arr(i,4), arr(i,5),arr(i,6), arr(i,2),arr(i,3), 'LineStyle','none','CapSize', 0);
end
corr(x,y);
vals = [x,y];
Fits = bootstrp(1000,@polyfits,vals);
med = median(Fits);
med(1)
p_val = sum(Fits(:,1) > 0)/1000
plot(x,polyval(med,x),'LineWidth',2);
title('high-v')
xlim([-.8 .8])
ylim([-.1 1])
xlabel('confidence reports')
ylabel('conf-acc correlation')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline
subplot(1, 3, 3)
x = zeros([10 1]);
y = zeros([10 1]);
index = 1;
for i = 1:3:30
    x(index) = arr(i,1);
    y(index) = arr(i,4);
    index = index + 1;
    scatter(arr(i,1),arr(i,4),60,'filled');
    hold on
    errorbar(arr(i,1),arr(i,4), arr(i,5),arr(i,6), arr(i,2),arr(i,3), 'LineStyle','none','CapSize', 0);
end
corr(x,y);
vals = [x,y];
Fits = bootstrp(1000,@polyfits,vals);
med = median(Fits);
med(1)
p_val = sum(Fits(:,1) > 0)/1000
plot(x,polyval(med,x),'LineWidth',2);
title('baseline')
xlim([-.8 .8])
ylim([-.1 1])
xlabel('confidence reports')
ylabel('conf-acc correlation')


function [Fits] = polyfits(vals)
    Fits = polyfit(vals(:,1),vals(:,2),1);
end