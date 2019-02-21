function [cumReward, profit] = MLE(varargin)

addpath('StepPolicies');

% MLE problem class with one parameter
% Input: steprule, numIterations, vector of tunable parameters for the 
% stepsizes, and number of sample paths 

% Output: 
% cumulative reward - vector of size 1xnumPaths, were each 
%  entry contains the sum of the MSE at each iteration in the sample path
% profit - vector of size 1xnumPaths, where each entry is the 
% negative of the final MSE (at the end of one sample path)

% read in inputs 
steprule = varargin{1};
numIterations = varargin{2};
tuneparam = varargin{3};  
numPaths = varargin{4};

namerule = func2str(steprule);

switch(namerule)
    % initialize parameters for adam
    case 'adam'
        mpast = 0;
        vpast = 0;
        gradF = 1;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            alphanought = tuneparam(1);
        else 
            alphanought = 0.00004;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            beta1 = tuneparam(2);
        else 
            beta1 = 0.95;
        end 
        if(tuneparam(3) ~= 0 && ~isnan(tuneparam(3)))
            beta2 = tuneparam(3);
        else 
            beta2 = 0.95;
        end 
    % initialize parameters for adagad
    case 'adagrad'
        gradF = 1;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            adagradstepsize = tuneparam(1);
        else 
            adagradstepsize = .0004;
        end 

    % initialize parameters for GHS
    case 'GHS'
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            GHSalpha = tuneparam(1);
        else 
            GHSalpha = .018;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            GHStheta = tuneparam(2);
        else 
            GHStheta = 500;
        end 

  % initialize parameters for Polynomial learning rates 
   case 'polylearning'        
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            Polyalpha = tuneparam(1);
        else 
            Polyalpha = .00003;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            Polybeta = tuneparam(2);
        else 
            Polybeta = 0.8;
        end 
  case 'kestens'
        % initialize parameters for Kestens 
        K = 0;
        prevgradF = 1;
        gradterm = 1;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            kestenalpha = tuneparam(1);
        else 
            kestenalpha = .0003;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            kestentheta = tuneparam(2);
        else 
            kestentheta = 10;
        end 

end 


% initialization of variables 
original = 10;
theta = original'; % values of initial parameters 
numparams = size(original, 2);
data = zeros(numIterations, numparams); %data matrix  
sigma2 = 10; % variance of noise
fvalues = zeros(1, numIterations);
err = zeros(1, numPaths); %% vector of final MSE for each sample path
profit = zeros(1, numPaths); %% -(1) * err
cumReward = zeros(1, numPaths); %% sum of the MSEs for each sample path
epsilon = 1e-8; % error term for stepsizes 


% for each sample path
for path = 1:numPaths 


    % for a particular sample path, keep track of the MSE at each 
    % iteration
    F = zeros(1, numIterations); 

    % generate data 
    xvect = zeros(1, numparams)';
    for i = 1:numIterations
            for j = 1:numparams 
                x = randi([1, 200]);
                xvect(j) = x;
            end 
            f = dot(theta, xvect);
            y = f + normrnd(0, sigma2^(.5));
            data(i,:) = xvect'; 
            fvalues(i) = y;
    end

     % use stochastic gradienet algorithm to find parameters
     % initialize our estimates of theta as zeros.
     est = zeros(1, numparams)'; 
     
%      estvect = zeros(numIterations, numparams); 
%      funct = zeros(numIterations); 
     
        for k = 1:numIterations 
            
            % choose the stepsize rule based on input 
            if namerule == string('GHS')
                a = steprule(GHSalpha, GHStheta, k);
            end 
            if namerule == string('adam')
                a = adam(j, alphanought, beta1, beta2, mpast, vpast, gradF,epsilon); 
            end 
            if namerule == string('polylearning')
                a = polylearning(Polyalpha, j, Polybeta);
            end 
            if namerule == string('adagrad')
                if (k == 1) 
                histgrad = zeros(1, 1);
                [a, g] = steprule(adagradstepsize, gradF,histgrad, 1, epsilon, 1);
                histgrad = g;
                else 
                [a, g] = steprule(adagradstepsize, gradF, histgrad, 1, epsilon, 1);
                histgrad = g;
                end  
            end 
           if namerule == string('kestens')
                if (k == 1) 
                     K = 0;
                end 
                    [a, newk] = steprule(kestenalpha, kestentheta, K, prevgradF,gradterm);
                    K = newk;    
           end 
            
            % for each of the parameters use the partial derivative
            % to update estimates of the parameters 
            for i = 1:numparams
                gradterm = computeGrad(fvalues(k), data(k, :), est, i);
                prevgradF = gradterm;
                est(i) = est(i) + a*gradterm;
                F(k) = immse(original',est); 
%                 estvect(path*k, i) = est(i);
            end

%            funct(path*k) = dot(data(k, :), est); 
            cumReward(1, path) = cumReward(1, path) + immse(original',est);
        end 
      
       err(path) = immse(original', est);
       profit(path) = immse(original', est); %% tack on a negative sign
end

end
    
    function gradterm = computeGrad(y, data, est,paranum)
    
        gradterm = (y - dot(est,data))* data(paranum);
    end 