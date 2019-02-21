function [step] = adam(n, alpha, beta1, beta2, mpast, vpast, gradt, epsilon)
% adam stepsize rule with bias correction mechanism 
% http://cs231n.github.io/neural-networks-3/

% Input 
% n - iteration number 
% alpha - learning rate 
% beta1 - decay factor for running average of the gradient
% beta2 - decay factor for running average of the gradient^2
% mpast - running average of gradient - think of as mean
% vpast - running average of gradient^2 - think of as measure of variance
% gradt - gradient 
% prevent division from zero error. (usually use10^-8)

m = beta1*mpast + (1-beta1)* gradt; % first moment estimation 
v = beta2*vpast + (1-beta2)* gradt*gradt; % second moment

mhat = m /(1-realpow(beta1,n));  
vhat = v / (1-realpow(beta2,n));

step = (alpha /(sqrt(vhat) + epsilon)) * mhat; 

end

