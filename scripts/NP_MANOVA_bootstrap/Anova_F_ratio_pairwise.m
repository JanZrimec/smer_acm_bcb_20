function [F, Swithin2, Stotal2] = Anova_F_ratio_pairwise(in,y)
% nonparametric anova anderson 2001 - Euclidean distance
% in .. input original data m observations x n variables
% y  .. class labels on observations

labels = unique(y); % unique classes
C = length(labels);

Swithin = zeros(C,1);
F = zeros(C);

for i = 1:C
    tmp = sum(power(pdist(in(y==labels(i),:)),2));
    Swithin(i) = tmp/sum(y==labels(i));
end

for i=1:C
    for j=i+1:C
        Swithin2 = Swithin(i)+Swithin(j);
        nm = sum(y==labels(i)|y==labels(j));
        
        tmp = sum(power(pdist(in(y==labels(i)|y==labels(j),:)),2));
        Stotal2 = tmp/nm;
        
        Sbetween2 = Stotal2 - Swithin2;
        F(i,j) = (Sbetween2/Swithin2)*(nm-2);
    end
end

end