function  [mu, beta_W, numD]=PayloadD(noiseRation,varargin)
%In this general Payload Delivary setting, there is two unknown parameters coming from some uniform distributions. 
%The truth is the expected function values over the distribution of the
%two unknown parameters  which is approximately obtained by sampling num=100 times.

%if nargin==0, return the expected function value used for prior generation
%and noise level specification
%if nargin==1, return the function value of a sampled set of parameters
%from the uniform distributions as the truth.
numD=2;
M=100;
MIN_PHI_W = 0.1;      % minimum value of phi_w (volume fraction of water droplets inside oil droplets)
MAX_PHI_W = 0.9;          % maximum value of phi_w
MIN_PHI_O = 0.1;      % minimum value of phi_o (volume fraction of oil droplets in solution)
MAX_PHI_O = 0.9;           % maximum value of phi_o
% Constants
MATERIALS_SCALE = 0.3;
num = 100;

% Discretization
NUM_PHI_W = 10; 
NUM_PHI_O = M/NUM_PHI_W;
STEP_PHI_W = (MAX_PHI_W - MIN_PHI_W)./(NUM_PHI_W - 1); 
STEP_PHI_O = (MAX_PHI_O - MIN_PHI_O)./(NUM_PHI_O - 1); 

%return the expected function value used for prior generation
              %and noise level specification
% flat indexing
% iterate over phi_w
    for i=0:(NUM_PHI_W-1)
        phi_w = MIN_PHI_W + i*STEP_PHI_W;

        % iterate over phi_o
        for j=0:(NUM_PHI_O-1)
            phi_o = MIN_PHI_O + j*STEP_PHI_O;

            % calculate utility function and store it in 1D index i*NUM_PHI_O + j
            X(1, i*NUM_PHI_O + j+1) = phi_w;
            X(2, i*NUM_PHI_O + j+1) = phi_o;
        end
    end
    for s=1:num
      theta(s,:) =  PayloadD_sample_parameters();
        % iterate over alternatives
        for i=1:M
            % get x and y
            x = X(1, i);
            y = X(2, i);

            U(i) = PayloadD_U(x, y, theta(s,:));
        end
        MATERIALS_VAR = MATERIALS_SCALE*(max(U) - min(U));

           % Generate samples

        % Sample a truth
        samples(s,:) = U + randn(size(U))*MATERIALS_VAR;
    end
%      mu=mean(samples)';
     beta_W=1./var(samples)';

    % Sample a truth
    theta =  PayloadD_sample_parameters();
    for i=0:(NUM_PHI_W-1)
        phi_w = MIN_PHI_W + i*STEP_PHI_W;

        % iterate over phi_o
        for j=0:(NUM_PHI_O-1)
            phi_o = MIN_PHI_O + j*STEP_PHI_O;

            % calculate utility function and store it in 1D index i*NUM_PHI_O + j
            mu(i*NUM_PHI_O + j+1, 1) = PayloadD_U(phi_w, phi_o, theta);
        end
    end


end



