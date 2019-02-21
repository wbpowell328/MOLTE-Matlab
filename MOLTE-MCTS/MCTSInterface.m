close all;
clear;

spreadsheet='MCTSProblems.xls';
[ndata, ~, alldata] = xlsread(spreadsheet,'','','basic');
[numPPp, ~]=size(alldata);
numPP=numPPp-1;

[~,~,pathdata] = xlsread(spreadsheet, 'Path', '', 'basic');
path = pathdata{2,1};
addpath(genpath(path));

for problemi = 1:numPP
    problemInstance = alldata{problemi+1,1};
    problemType     = alldata{problemi+1,2};
    rolloutPolicy   = str2func(alldata{problemi+1,3});
    d_thr           = alldata{problemi+1,4};
    e_thr           = alldata{problemi+1,5};
    alpha           = alldata{problemi+1,6};
    budget          = alldata{problemi+1,7};
end

if strcmp(problemInstance, 'TSP')    
    % Load Java dependency
    dpath = {strcat(path, '/TSP/EqualByValueDoubleSet.jar')};
    javaclasspath(dpath);
    results = TSPf(problemType, rolloutPolicy, d_thr, e_thr, alpha, budget);
elseif strcmp(problemInstance, 'Selling')
    results = SellingProblemf(problemType, rolloutPolicy, d_thr, e_thr, alpha, budget);
elseif strcmp(problemInstance, 'Pricing')
    % Set param (weigher) to 1 if using KGCBLin, otherwise read in param 
    % (epsilon) for epsilon greedy algorithm 
    if strcmp(func2str(rolloutPolicy), 'KGCBLinPricing')
        param = 1;
    else
        param = alldata{problemi+1, 10};
    end
    
    % Load problem specific data
    numExp = alldata{problemi+1, 8};
    L = alldata{problemi+1, 9};
    resolution = alldata{problemi+1, 11};
    weight = alldata{problemi+1,12};
    
    results = PricingProblemf(problemType, rolloutPolicy, d_thr, e_thr, alpha, budget, numExp, L, param, resolution, weight);
else
    fprintf('Error: Problem %s does not exist!\n', problemInstance);
end