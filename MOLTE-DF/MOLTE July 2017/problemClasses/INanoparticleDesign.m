function [mu,beta_W, numD]=INanoparticleDesign(noiseRation,varargin)
% truth-from-prior experiments
     mu_0=varargin{1};
     covM=varargin{2};
     beta_W=varargin{3};
     numD=2;
     mu=mvnrnd(mu_0, covM)';
end