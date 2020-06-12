function [out,outboot] = norm_to_boots(in,inboot)
% norm F statistic according to distribution of bootstrap replicates of F
% in ... F statistic
% inboot ... vector of bootstrap replicates of F statistic

% force column vector
inboot = inboot(:);

amin = min([in;inboot]);
amax = max([in;inboot]);
out = featureNormalize01(in,amin,amax);
outboot = featureNormalize01(inboot,amin,amax);

end