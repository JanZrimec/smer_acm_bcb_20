function [F, Swithin2, Stotal2] = Anova_F_ratio_dist_pairwise(in,y)
% nonparametric anova anderson 2001 - Euclidean distance
% in .. input distance between data m observations x n variables
% y  .. class labels on observations

% output equal to t-test if squared

labels = unique(y); % unique classes
C = length(labels);

Swithin = zeros(C,1);
F = zeros(C);

for i = 1:C
    tmp = in(y==labels(i),y==labels(i));
    tmp2 = triu(ones(sum(y==labels(i))),1);
    tmp3 = sum(power(tmp(logical(tmp2)),2));
    Swithin(i) = tmp3/sum(y==labels(i));
end
% tmp2 = triu(ones(nm),1);
% Stotal = sum(power(in(logical(tmp2)),2))/nm;

for i=1:C
    for j=i+1:C
        Swithin2 = Swithin(i)+Swithin(j);
        nm = sum(y==labels(i)|y==labels(j));
        
        tmp = in(y==labels(i)|y==labels(j),y==labels(i)|y==labels(j));
        tmp2 = triu(ones(sum(y==labels(i)|y==labels(j))),1);
        tmp3 = sum(power(tmp(logical(tmp2)),2));
        Stotal2 = tmp3/nm;
        
        Sbetween2 = Stotal2 - Swithin2;
        F(i,j) = (Sbetween2/Swithin2)*(nm-2);
    end
end
end