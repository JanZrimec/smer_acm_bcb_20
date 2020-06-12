function out=pvalue_boots(in,inboot,num_boots)
% calculate p-values from bootstrap replicates of F statistic
% in ... F statistic
% inboot ... vector of bootstrap replicates of F statistic
% num_boots ... size of inboot

out = sum(inboot > in)/num_boots;

end
