function [F, Swithin, Stotal] = Anova_F_ratio(in,y)
% nonparametric anova anderson 2001 - Euclidean distance
% in .. input original data m observations x n variables
% y  .. class labels on observations

labels = unique(y); % unique classes
C = length(labels);
[nm,nn] = size(in);

Swithin = 0;
Sbetween = 0; % Ctotal - Cwithin
Stotal = 0;

for i = 1:C
    tmp = sum(power(pdist(in(y==labels(i),:)),2));
    Swithin = Swithin + tmp/sum(y==labels(i));
end
Stotal = sum(power(pdist(in),2))/nm;
Sbetween = Stotal - Swithin;

F = (Sbetween/(C-1))/(Swithin/(nm-C));
end