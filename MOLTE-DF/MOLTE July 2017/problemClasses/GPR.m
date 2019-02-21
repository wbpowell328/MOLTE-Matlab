function [mu, beta_W, numD]=GPR(noiseRation,varargin)

     mu_0=varargin{1};
     covM=varargin{2};
     beta_W=varargin{3};
     mu=mvnrnd(mu_0, covM)';
     numD=1;
     
end