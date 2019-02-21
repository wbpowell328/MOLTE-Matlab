%{
This function generates a prior (using MLE if Prior=='MLE') of the problemClass and returns the external noise level
beta_W(=1/sigma_W^2). numD is the dimension of the text function.
%}
function [mu_0, covM, beta_W] = priorTruthGen(problemClass, noiseRatio, bm,paras, Prior)
% bm: belief model, whether independent or correlated

if strcmp(Prior, 'MLE')  
        %Use MLE to fit a prior distribution
%         [mu, beta_W, numD]=truthGen(problemClass);   
        [mu_0, covM, beta_W] = MLEpriorGen( problemClass,noiseRatio, bm);
%        [mu_0, covM] = EmPriorGen( mu, beta_W, bm, numD);
else if strcmp(Prior, 'uninformative')
         [mu, beta_W]=truthGen(problemClass,noiseRatio);   
         M=length(mu);
         mu_0=zeros(M,1);
         covM=zeros(M,M);
         for i=1:M
             covM(i,i)=10^10;
         end
        
    else     % either use a default prior or use the prior given by the user
        [mu_0, covM, beta_W]= priorGen(problemClass,noiseRatio, Prior, paras); 
    end

end
end


