%=========== Description: ===========
%
% Modification of maniMDS of Chen and Muller
%
%================ Input Arguments: ================
%
%Input X_reg: 	matrix of dimension NxN where N is the number of subjects
% and M is the number of points in the grid T 
%
%Input	T: 	grid of M points 
%
%
%Input Dis: pairwise L2 distances 
%
%================ Output Arguments: ================
%
% Output : 
%
% deltacv : best delta chosen with cross-validation
% 
% manidis : matrix containing pairwise geo estimation

function [deltacv,Manidis] = mani_p_iso(X,T)

N = length(X);
t_pooled = unique(cell2mat(T));

t = t_pooled;
M = length(t);
X_reg = reshape(cell2mat(X),M,N)';
X_cv = X_reg;
Dis = L2_distance(X_reg',X_reg',1)*sqrt(range(t)/(M-1));

Dis = (Dis+Dis')/2;

[Kcv,deltacv,hcv,mse] = par10cv_mod(X,T,5,X_cv,t_pooled,'epan',Dis,'p-isomap',4);


[Y,Index,FVE,Manidis] = maniMethods(Dis,'p-isomap',4,Kcv,deltacv);


end
