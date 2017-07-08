function b = multivariate_kurtosis(data)
% Authors: Hua (Oliver) Xie & Jason E. Hill
% Date 17 JAN 2017
% multivariate_kurtosis computes the multivariate kurtosis 

nT = size(data,1);

C = cov(data,1);

d = abs((data/C)*data');

b = sum(diag(d).^2)/nT;
