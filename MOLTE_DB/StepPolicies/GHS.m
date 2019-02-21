function [alpha] = GHS(alphanought, a, n)
% Generalized Harmonic Stepsize 
% returns stepsize of function with tuning parameter 
% a, initial alpha, and iteration number

% Output: 
% alpha - stepsize at iteration n

% Input: 
% a - tuning parameter 
% alphanought - stepsize/learning rate 
% n - iteration number 
if (n == 0) 
    error('N must be greater than zero');
end 
if (a == 0 || alphanought == 0) 
    error('theta and alphanought cannot equal zero');
end 
    
      alpha = alphanought * (1/(a + n - 1));
%     alpha = alphanought * (a/(a + n - 1));
end 




