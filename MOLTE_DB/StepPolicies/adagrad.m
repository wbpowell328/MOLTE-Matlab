function [alpha, G] = adagrad(varargin) 
% let our stochastic gradient algorithm be
% given by theta = theta + n* adjusted_grad
% where adjusted_grad= (grad/(epsilon + sqrt(sumofhistoricalgrad)) 

% the user must create a diagonal gradient nxn matrix where n 
% is the number of parameters

% Input: 
% 
% grad: function handle to stochatic gradient
% stepSize: scalar step size 
% epsilon: error term 
% i: index of gradient to use 
% 
% 
% Output: 
% alpha where alpha * grad is the "adjusted gradient"
% G : matrix where the diagonals give gradients 


% adagrad(stepsize, gradient,histgrad, whichparameter, epsilon, numparams)

%  if (varargin{6} == 1) 
     stepSize = varargin{1};
     grad = varargin{2};
     histgrad = varargin{3};
     i = varargin{4};
     epsilon = varargin{5};
     numparams = varargin{6};
     histgrad = histgrad(i, i) + grad*grad;
     G = histgrad;
     alpha = stepSize / (sqrt(histgrad) + epsilon);
     
%  elseif (varargin{6} > 1)
%      stepSize = varargin{1};
%      grad = varargin{2};
%      histgrad = varargin{3};
%      i = varargin{4};
%      epsilon = varargin{5};
%      numparams = varargin{6};
%      G = updateG(histgrad, grad, i);
%      alpha = stepSize / (sqrt(G(i,i)) + epsilon);     
%  else 
%      error('adagrad requires 4 or 5 parameters');
%  end 

end

% % method to update the gradient according to the adagrad rule 
% function G = updateG(G, grad, i) 
%     G(i,i) = G(i, i) + grad*grad;    
% end 

