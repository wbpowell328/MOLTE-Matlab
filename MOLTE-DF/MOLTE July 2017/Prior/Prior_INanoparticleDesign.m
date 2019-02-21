
%  If nargin == 0, use default prior;
function [mu_0, covM, beta_W]=Prior_INanoparticleDesign(varargin)
noiselevel=0.01;
priorMatrix=load('Default_INanoparticleDesign.mat');
mu_0=priorMatrix.mu_0;
covM=priorMatrix.covM;
beta_W=ones(size(mu_0))./(noiselevel*range(mu_0))^2;
end

