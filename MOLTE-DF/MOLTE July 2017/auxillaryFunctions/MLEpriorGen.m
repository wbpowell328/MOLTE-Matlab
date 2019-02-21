function [ mu_0, covM, beta_W] = MLEpriorGen(problemClass, noiseRatio, bm)
% Use Latin hypecube designs to draw the initial points and use MLE to construct a uniform prior distribution.
% For correlated beliefs, we use 
% a exponential decay prior covariance matrix. We use MLE to fit this 
% covariance structure using the initial sample points. 

ox=8; 
tx=12;

[mu, beta_W, numD]=truthGen(problemClass, noiseRatio);   
M=length(mu);
covM=zeros(M,M);
if strcmp('independent', lower(bm))==1
    index=ceil(lhsdesign(ox,1).*M);
    sample=zeros(ox,1);
    for i=1:ox
        %generate the truth values
        mu=truthGen(problemClass, noiseRatio);   
        sample(i)=mu(index(i));
    end
    mu_0=mean(sample)*ones(M,1);
    sigma_i=var(sample);
    for i=1:M
       covM(i,i)=sigma_i;
    end
end

if strcmp('correlated', lower(bm))==1
      if numD==1
        index=ceil(lhsdesign(ox,1).*M);
        sample=zeros(ox,1);
        for i=1:ox
            mu=truthGen(problemClass, noiseRatio);   
            sample(i)=mu(index(i));
        end
        [alpha betas ave] = MLE( index , [], sample );
        mu_0=ave*ones(M,1);
        %construct covariance matrix
        for i=1:M
            for j=1:M
                covM(i,j)=alpha*exp(-betas(1)*(i-j)^2 );
            end
        end
      end
    
    
   if numD==2

        lhX=ceil(lhsdesign(tx,2)*sqrt(M));
        x=lhX(:,1);
        y=lhX(:,2);
        index=x+(y-1)*sqrt(M);
        index_reshape=reshape(index,[tx,1]);

        sample=zeros(tx,1);
        for i=1:tx
            mu=truthGen(problemClass,noiseRatio);   
            sample(i)=mu(index_reshape(i));
        end
        [alpha betas ave] = MLE( x , y, sample );
        mu_0=ave*ones(M,1);
        
        %construct covariance matrix
        xx=1:sqrt(M); yy=1:sqrt(M);
        [X,Y]=meshgrid(xx,yy);
        x_reshape=reshape(X,[M,1]);
        y_reshape=reshape(Y,[M,1]);
        for i=1:M
            for j=1:M
                covM(i,j)=alpha*exp(-betas(1)*(x_reshape(i)-x_reshape(j))^2 -betas(2)*(y_reshape(i)-y_reshape(j))^2 );
            end
        end
    end
    
end

end

