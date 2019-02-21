%{
Tune the parameter of the policy within the problemClass using the offline
objective funciton: the expected reward of the final recommended alternative.
%}
function alpha_tuned=offtune( policy, problemClass,noiseRatio, M, N, bm,paras, Prior)
% INPUT:
% policy: the policy we are interested in tuning the parameter alpha
% problemClass: the name of the problemClass
% M: number of alternatives
% N: measurement budget
% bm: belief model, whether independent or correlated
% paras: the parameters for the problem class
% Prior: whether to use a default prior/ a given prior, or fit a prior
% using MLE

%number of (truth, prior, sample paths);
numSamples=200; % a reasonable value for approximating the expected total rewards is 500 - 1000


RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

%generate 'numSamples' random (truth, prior, sample paths);
mu_0=zeros(M, numSamples);
covM=zeros(M,M,numSamples);
mu=zeros(M, numSamples);
beta_W=zeros(M, numSamples);
samples=zeros(M,N, numSamples);


for i=1:numSamples  % generate numSamples (truth, prior, samples)

    [mu_01, covM1,beta_W1] = priorTruthGen(problemClass,noiseRatio, bm, paras, Prior); %generate a prior
    %generate one possible truth      
    [mu1,beta_W1]=truthGen(problemClass,noiseRatio, mu_01, covM1,beta_W1);     
    mu(:,i)=mu1; beta_W(:,i)= beta_W1;  
    mu_0(:,i)=mu_01;
    covM(:,:,i)=covM1;
    [samples(:,:,i) ] = sampleGen( mu(:,i), beta_W(:,i), N );  %generate sample paths;
end

total_ereward=zeros(11,1);
for power=-5:5
    alpha=10^power;
    for i=1: numSamples      
           [mu_est,recommendedArm]=policyT(policy,mu_0(:,i),beta_W(:,i),covM(:,:,i),samples(:,:,i), mu(:,i), alpha);
           total_ereward(power+6)=total_ereward(power+6)+mu(recommendedArm,i)/numSamples; % offline objective: only care about the final recommendedArm
    end
end
[reward_max, alpha_max]=max(total_ereward);
l_range=10^-5;
if alpha_max>1
   l_range=10^(alpha_max-6-1);
end
r_range=10^(alpha_max-6+1);
total_ereward=zeros(11,1);
for t=1:3 %repeat search for 3 times
    for index=1:11
  
        alpha=l_range+(r_range-l_range)/10*(index-1);  %test the range in 10 percent increments to find the best alpha
         for i=1: numSamples
      
                 [mu_est,recommendedArm]=policyT(policy,mu_0(:,i),beta_W(:,i),covM(:,:,i),samples(:,:,i), mu(:,i), alpha);
                 total_ereward(index)=total_ereward(index)+mu(recommendedArm,i)/numSamples; % offline objective: only care about the final recommendedArm
         end
    end
    [reward_max, alpha_max]=max(total_ereward);
    if alpha_max>1
       l_range=l_range+(r_range-l_range)/10*(alpha_max-2);
    end
    r_range=l_range+(r_range-l_range)/10*alpha_max;
    total_ereward=zeros(11,1);
    alpha_tuned=l_range+(r_range-l_range)/10*(alpha_max-1);
end



end

