function [X_norm] = featureNormalize01(X,minx,maxx)
%FEATURENORMALIZE Normalizes the features in X 

[m,n] = size(X);
X_norm = (X-repmat(minx,m,n))./(repmat(maxx-minx,m,n));

end