% Description: MOLTE_DB calls different stepsize policies and stochastic
% optimization problems and compares performances using a histogram of 
% profits. 

% INPUT:
% spreedsheet - the file name of the spreedsheet that contains different 
%             problem classes

spreadsheet= 'MOLTE-DB-spreadsheet.xlsx'; 
[ndata, text, alldata] = xlsread(spreadsheet,'','','basic');
[numPPp, numA]=size(ndata);
numPP=numPPp;
numPaths = 10; % number of times we run a sample path

addpath('StepPolicies');
addpath('Problems');
% loop through the problems ie. the rows of the spreadsheet
for problemi = 1:numPP 

problem = alldata{problemi+1, 1};

% number of iterations in one sample path 
numSim = alldata{problemi+1, 2}; 

% get problem type (ex. newsvendor)
problemtype = alldata{problemi+1, 3}; 

numberofpolicies=alldata{problemi+1,4};

pol=strings(numberofpolicies,1);


% read in the policies as strings 
for i=1:numberofpolicies
%     pol(i)= string(text(2,i+2));
    pol(i)= string(text(problemi+1,i+4));
end


% transform input string to function handle so we can call stepsize
% functions 
problemclass = str2func(problem);
policies=alldata(problemi+1,5:4+numberofpolicies);


% get default parameters inputted through spreadsheet, if they exist 
tuned = zeros(1,numberofpolicies);
paraProvided=zeros(1, numberofpolicies);
alpha=zeros(3,numberofpolicies); % store the value of tunable parameters for each policy

for jjj=1:numberofpolicies % find out the policies name and whether they need to be tuned
    charp=char(policies(jjj));
    temppp=length(charp);
    if charp(temppp)==')'   
        before=find(charp=='(');
        if charp(temppp-1)=='*'  
            tuned(jjj)=1;
         else
%             alpha(jjj)=str2double(charp(before+1:temppp-1));
            comma1 = find(charp == ',', 1, 'first');
            comma2 = find(charp == ',', 1, 'last');
            alpha(1, jjj) = str2double(charp(before+1:comma1-1));
            alpha(2, jjj) = str2double(charp(comma1+1:comma2-1));
            alpha(3, jjj) = str2double(charp(comma2+1:char(temppp)-1));
            paraProvided(jjj)=1;
         end
    end

end

% for running time purposes, we only do one sample path of energy inventory
% if contains(problem, 'energyinventory')
%     profitmatrix = zeros(numberofpolicies, 2);
% else 
    profitmatrix=zeros(numberofpolicies, numPaths); % init profit matrix 
% end 

% loop through the stepsize rules 
% get the profit and cumreward 
% for MLE, this profit is the MSE 
for i=1:numberofpolicies

    if contains(pol(i), 'adam')
        [cumreward, profit] = problemclass(@adam, numSim, alpha(:,i), numPaths);
        if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
    if contains(pol(i), 'GHS')
        [cumreward, profit] = problemclass(@GHS, numSim, alpha(:,i), numPaths);
        if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
    if contains(pol(i), 'polylearning') 
        [cumreward, profit] = problemclass(@polylearning, numSim, alpha(:,i), numPaths);
        if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
    if contains(pol(i), 'adagrad')
       [cumreward, profit] = problemclass(@adagrad, numSim, alpha(:,i), numPaths);
       if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
    if contains(pol(i), 'BAKF')
       [cumreward, profit] = problemclass(@BAKF, numSim, alpha(:,i), numPaths);
        if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
    if contains(pol(i), 'kestens')
        [cumreward, profit] = problemclass(@kestens, numSim, alpha(:,i), numPaths);
        if strcmp(problemtype, 'online') && ~contains(problem, 'energyinventory')
            profitmatrix(i,:)= cumreward;
        else
            profitmatrix(i,:)= profit;
        end 
    end
end

% if contains(problem, 'energyinventory')
%     numPaths = 2;
% end 

referenceprofit=0;

% our reference profit is how much profit is generated 
% using the first stepsize rule/policy 
for i=1:numPaths
    referenceprofit= referenceprofit+profitmatrix(1,i);
end

mean_referenceprofit=referenceprofit/numPaths;
compprofitmatrix=zeros(numberofpolicies-1,numPaths);

% compute "profit"or "MSE"
% positive equals better compared to policy in first column of spreadsheet
for i=1:(numberofpolicies-1)
    for j=1:numPaths
        compprofitmatrix(i,j)=profitmatrix(i+1,j)-mean_referenceprofit;
    end
end

trans= zeros(numPaths,numberofpolicies-1);

for i=1:numPaths
    for j=1:(numberofpolicies-1)
        trans(i,j)= compprofitmatrix(j,i);
    end
end

% create histograms 
color=[1,0,0;0,0.95,0.95;0.6,0.3,0.9;0,1,0;1,1,0;1,0.25,0.75;0,0.5,0.75;0.5,0.95,0.75;];
firstPolicy=policies{1};
remainingPolicies=policies(2:length(policies));
figure
% hist(trans(:,1:(numberofpolicies-1)),200);

hist(trans(:,1:(numberofpolicies-1)), 200);

% this is to make sure the y limits are correct for energy inventory
% histogram since we only obtain one profit for each policy
% if contains(problem, 'energyinventory')
%    ylim([0 1]);
% end
colormap (color);
legend(remainingPolicies, 'Location', 'best');

end 
