%{
notation for the following:
M is the number of alternatives.
N is the number of time-steps
M x N stands for a matrix with M rows and N columns

INPUT:
mu_0:   prior for the mean (M x 1)
beta_W: measurement precision (1/lambda(x)) (M x 1)
covM:   initial covariance matrix (M,M)
samples: pre-generated observation realizations (M x N)
alpha: tunable parameter
tune: whether tune the parameter or not: if tune=0, use the default value

OUTPUT:
mu_est:     Final estimates for the means (M x 1)
count:      The number of each alternative being measured (M x 1)
recommendedArm:  The index of the arm recommended by this policy after the
measurement budget exhausted
%}

function [ mu_est, count,recommendedArm ] = OLKG( mu_0,beta_W,covM,samples,  alpha, tune)
[M,N]=size(samples);
count=zeros(M,1);  % count the times each alternative is measured
choice=zeros(1, N);
reward_sum=0;
ereward_sum=0;

mu_est=mu_0;
beta_est=1./diag(covM,0);
 
 for i=1:N
     sigmatilde=sqrt(1./beta_est-1./(beta_est+beta_W));
     %calculate zeta
     aux=repmat(mu_est, 1, M);
     aux1=aux(diag(ones(M,1))~=1);
     aux2= reshape(aux1, M-1,M);
     zeta=-abs((mu_est-max(aux2)')./ sigmatilde);
     f=zeta.*normcdf(zeta)+normpdf(zeta);
     KGfactor=  sigmatilde.*f;
     [max_est, x]=max(mu_est+(N-i+1)*KGfactor);
     count(x)=count(x)+1;
     W=samples(x, count(x));
   
     %update belief
     mu_est(x)=(beta_est(x)*mu_est(x)+beta_W(x)*W)/(beta_est(x)+beta_W(x));
     beta_est(x)=beta_est(x)+beta_W(x);

 end
 [value, recommendedArm] = max(mu_est); %the recommended arm is the last surviving arm. 
end

