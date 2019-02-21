function [finalcumulative, profit] = newsvendor3(varargin)
% newsvendor problem with demand exponentially distributed and price much
% larger than cost 

% Input:
% numTruth: number of times to run one sample path of newsvendor simulation 
% steprule: alpha 

% Output:
% y: vector of estimates of our optimal demand 
addpath('StepPolicies');

steprule = varargin{1};
numIterations = varargin{2};
tuneparam = varargin{3};
numPaths = varargin{4};
% numPaths = 5;

% initialization of newsvendor problem parameters 
c = 100; % cost of newspaper
p = 120; % price of newspaper
mu = 120; % mean of demand distribution 
% sigma = 20; % variance of demand distribution 
grad = @vanillanewsvendorgrad; % gradient 
epsilon = 1e-8;

namerule = func2str(steprule);

switch(namerule) 
    % initialize parameters for stepsize rules 
    case 'BAKF'
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            lambda = tuneparam(1);
        else 
            lambda = 1;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            BAKFalpha = tuneparam(2);
        else 
            BAKFalpha = 1;
        end 
        if(tuneparam(3) ~= 0 && ~isnan(tuneparam(3)))
            alpha = tuneparam(3);
        else 
            alpha = 1;
        end 
        nu = .1;
        v = 1; 
        beta = .05;
    case 'adam'

        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            alphanought = tuneparam(1);
        else 
            alphanought = .2;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            beta1 = tuneparam(2);
        else 
            beta1 = 0.95;
        end 
        if(tuneparam(3) ~= 0 && ~isnan(tuneparam(3)))
            beta2 = tuneparam(3);
        else 
            beta2 = 0.80;
        end 

        mpast = 15;
        vpast = 10;
        
    case 'adagrad'
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            adagradstepsize = tuneparam(1);
        else 
            adagradstepsize = 60;
        end 
    
    case 'GHS'
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            GHSalpha = tuneparam(1);
        else 
            GHSalpha = 30;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            GHStheta = tuneparam(2);
        else 
            GHStheta = 1;
        end 
   case 'polylearning'
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            Polyalpha = tuneparam(1);
        else 
            Polyalpha = 30;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            Polybeta = tuneparam(2);
        else 
            Polybeta = 0.9;
        end 
  case 'kestens'
        K = 0;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            kestenalpha = tuneparam(1);
        else 
            kestenalpha = 60;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            kestentheta = tuneparam(2);
        else 
            kestentheta = 50;
        end 

end 

N = numIterations;

profit = zeros(1, numPaths);
cumaward = zeros(1, numPaths);

for i = 1:numPaths
    
% initialization of variables 
x = 0; 
y = zeros(1, N); % vector of estimates of parameter 
steps = zeros(1, N); % store the stepsizes 
F = zeros(1, N); % profit vector for a single sample path
gradF = 1;
gradFvect = zeros(1, N);

  
% iterate through one sample path
    for j = 1:N

     w = exprnd(mu); % distribution of demand 

     F(j) = p*min(x, w) - c*x; % compute profit 
     cumaward(i) = cumaward(i) + F(j);
     prevgradF = gradF;
     gradF = grad(w, x, p, c); % get gradient 
     gradFvect(j) = gradF;
     xprev = x;

     % ADAM
     if namerule == string('adam')
        alpha = steprule(j, alphanought, beta1, beta2, mpast, vpast, gradF, epsilon);
     end 
     % GHS
     if namerule == string('GHS')
       alpha = steprule(GHSalpha, GHStheta, j); 
     end 
     % Polynomial learning rates
     if namerule == string('polylearning')
       alpha = steprule(Polyalpha, j, Polybeta);    
     end 
     % adagrad 
     if namerule == string('adagrad')
        if (j == 1) 
            histgrad = 0;
            [alpha, g] = steprule(adagradstepsize, gradF, histgrad, 1, epsilon, 1);
            histgrad = g;
        else 
        [alpha, g] = steprule(adagradstepsize, gradF, histgrad, 1, epsilon, 1);
        histgrad = g;
        end 
     end 

     % BAKF 
     if namerule == string('BAKF') 
         x = x + alpha*gradF;
        [alpha, beta, v, lambda] = steprule(j, x, xprev, nu, beta, v, lambda, BAKFalpha);  
     end 

     % kestens  
     if namerule == string('kestens')

        if (j == 1) 
            K = 0;
        end 
        [alpha, newk] = steprule(kestenalpha, kestentheta, K, prevgradF,gradF);
        K = newk;    
     end 
     
     if namerule ~= string('BAKF')
        x = x + alpha*gradF;
     end 
     
     y(j) = x;
     steps(j) = alpha;
     
    end
        finprofit = F(N);
        profit(i) = finprofit; 
end
    finalcumulative = cumaward; 
%      plot(F);
end 