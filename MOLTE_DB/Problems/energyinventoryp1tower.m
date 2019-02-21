

% given battery storage pricing data find optimal set of sell/buy prices 
% such that profit is maximized

% in this version the price process is given by 
% a jump diffusion martingale process where: 
% P^t+1 = P^t + .5(energycost - price) + epsilon + jumpdiffusionfactor 

function [maxtheta, profit] = energyinventoryp1tower(varargin)

addpath('StepPolicies');
numrestart = 5; % number of random restarts 

steprule = varargin{1};
numIterations = varargin{2};
tuneparam = varargin{3};  
numPaths = varargin{4};

epsilon = 1e-8;

namerule = func2str(steprule);

F = zeros(1, numIterations);
switch(namerule)
    % initialize parameters for adam
    case 'adam'
    mpast = 0;
    vpast = 0;
    gradF = 1;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            alphanought = tuneparam(1);
        else 
            alphanought = 0.001;
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
    case 'adagrad'
        % initialize parameters for adagad
        % adagradstepsize = 30;
        gradF = 1;
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            adagradstepsize = tuneparam(1);
        else 
            adagradstepsize = .0001;
        end 
    
    case 'GHS'
        % initialize parameters for GHS
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            GHSalpha = tuneparam(1);
        else 
            GHSalpha = .0001;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            GHStheta = tuneparam(2);
        else 
            GHStheta = 1;
        end 
   case 'polylearning'
        % initialize parameters for Polynomial learning rates 
        if(tuneparam(1) ~= 0 && ~isnan(tuneparam(1)))
            Polyalpha = tuneparam(1);
        else 
            Polyalpha = .0001;
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
            kestenalpha = .001;
        end 
        if(tuneparam(2) ~= 0 && ~isnan(tuneparam(2)))
            kestentheta = tuneparam(2);
        else 
            kestentheta = 10;
        end 

end 


maxprofit = -inf;
maxtheta = zeros(1,2)';
% prices = findprices();

gradvect = zeros(2,1); % initialize gradient vector  
prevgrad = [1, 1];

N = numIterations;

profit = zeros(1, numPaths);

for path = 1:numPaths

for k = 1:numrestart
theta = getRandom()'; % initialize theta_S, theta_B

thetas = zeros(1, N);
thetab = zeros(1, N);

for i = 1:N 
    % GHS 
    
    switch(namerule) 
        case 'GHS' 
            a = steprule(GHSalpha, GHStheta, i);    
        case 'adam'
            a = adam(i, alphanought, beta1, beta2, mpast, vpast, gradF,epsilon); 
        case 'polylearning' 
            a = polylearning(Polyalpha, i, Polybeta);
        case 'adagrad'
            if (i == 1) 
            histgrad = zeros(1, 1);
            [a, g] = steprule(adagradstepsize, gradF,histgrad, 1, epsilon, 1);
            histgrad = g;
            else 
            [a, g] = steprule(adagradstepsize, gradF, histgrad, 1, epsilon, 1);
            histgrad = g;
            end 
        case 'kestens' 
            if (i == 1) 
                K = 0;
            end 
            [a, newk] = steprule(kestenalpha, kestentheta, K, prevgradF,gradterm);
            K = newk;         
    end     
    
    gradvect = getgradF(theta);
    theta = theta + a*gradvect;
    
    
    % safeguard against sell/buy price dipping below zero
    if (theta(1) < 0)
        theta(1) = 0;
    end 
    if (theta(2) < 0)
        theta(2) = 0;
    end 
    if (theta(1) < theta(2))
        theta(2) = theta(1);
    end 
    thetas(i) = theta(1); 
    thetab(i) = theta(2); 
    thetaest = [theta(1), theta(2)]';
    F(i) = getF(thetaest);
    
end 
    thetatry = theta; 
    [~, profittry] = getF(thetatry);
    % update profit if new starting point produces better local max
    if (profittry >= maxprofit) 
        maxprofit = profittry; 
        maxtheta = thetatry;
        maxF = F;
    end 
end 

    MAX = maxprofit;
    disp(MAX);
%     plot(maxF);
    profit(path) = MAX;
end 
end 
% compute numerical gradient 
% partial derivatives are computed by perturbing the current sell and 
% buy prices by delta and computing the resulting profits with 
% perturbed sell and buy prices and then computing 
% a slopes
function gradient = getgradF(theta)
    delta = 5; 
    gradient = zeros(2,1);
    
    % perturb buy/sell prices and obtain new profits
    thetaSperturb = theta + [1 0]'* delta;
    thetaBperturb = theta + [0 1]'* delta;
    [Fs, ~] = getF(thetaSperturb); 
    [Fb, ~] = getF(thetaBperturb); 
    [F, ~] = getF(theta);

    % calculate numerical partial derivatives
    gradFs = (Fs-F)/delta;
    gradFb = (Fb-F)/delta;
    gradient(1) = gradFs; 
    gradient(2) = gradFb;
end 

% get profit 
function [F, finalprofit] = getF(theta) 
    
    T = 1000; %% simulate for # iterations equal to # prices
    sellprice = theta(1);
    buyprice = theta(2);
    R = zeros(1,1000);
    inventory = 50;
    
    profit = 0;
    F = 0;
    price = 0;
    
    for t = 1:T
        
        % martingale with jumpdiffusion 
        energycost = 25; 
        epsilon = normrnd(0,10); 
        jumpprob = .05; 
        jumpdiffusion = normrnd(100,50);

        x = rand(); 
        if (x <= jumpprob)
            price = price + .5 * (energycost - price) + epsilon + jumpdiffusion; 
        else 
            price = price + .5 * (energycost - price) + epsilon; 
        end 

        
        if(isnan(price)) 
            R(t) = inventory;
            continue;
        end 
        if(price < 0) 
            R(t) = inventory;
            continue;
        end 
        
        decision = getDecision(price, sellprice, buyprice,inventory);

        if (inventory > 100) 
            inventory = 100;
        end 
        if (inventory < 0) 
            inventory = 0;
        end 
        % sell
        if (decision == -1) 
           inventory = inventory - 1;
           profit = profit + .95*price;
           F = F + .95*price;
        % buy
        elseif (decision == 1) 
            inventory = inventory + .95;
            profit = profit - price;
            F = F - .95*price;
        else 
        end 
        R(t) = inventory;
    end 
    finalprofit = profit;
end 

function findecision = getDecision(price, sellprice, buyprice, R)

    if (price > sellprice && R > 0) 
        decision = -1;
    elseif (price < buyprice)
        decision = 1;
    else
        decision = 0;
    end 
    findecision = decision;
    
end 
% generate a random set of sell/buy prices 
function randtheta = getRandom() 

    thetas = randi([10, 100]);
    thetab = randi([10,100]); 
    while ((thetas < thetab) || (thetas <0) || (thetab< 0))
        thetas = randi([5, 100]);
        thetab = randi([5,100]); 
    end 
    randtheta(1) = thetas; 
    randtheta(2) = thetab;
end
