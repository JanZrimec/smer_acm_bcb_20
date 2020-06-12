function [F, Swithin, Stotal] = Anova_F_ratio_dist(in,y)
% nonparametric anova anderson 2001 - Euclidean distance
% in .. input distance between data m observations x n variables
% y  .. class labels on observations

labels = unique(y); % unique classes
C = length(labels);
[nm,nn] = size(in);

Swithin = 0;
Sbetween = 0; % Ctotal - Cwithin
Stotal = 0;

for i = 1:C
    tmp = in(y==labels(i),y==labels(i));
    tmp2 = triu(ones(sum(y==labels(i))),1);
    tmp3 = sum(power(tmp(logical(tmp2)),2));
    Swithin = Swithin + tmp3/sum(y==labels(i));
end
tmp2 = triu(ones(nm),1);
Stotal = sum(power(in(logical(tmp2)),2))/nm;
Sbetween = Stotal - Swithin;

F = (Sbetween/(C-1))/(Swithin/(nm-C));
end