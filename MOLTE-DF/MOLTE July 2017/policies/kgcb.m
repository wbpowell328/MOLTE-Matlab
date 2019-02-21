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
function [ mu_est, count, recommendedArm ] = kgcb( mu_0,beta_W,covM,samples, alpha, tune)

[M, N]=size(samples);
 % M number of available choices, N number of measurements

mu_est=mu_0;

count=zeros(M,1);  % count the times each alternative is measured

for k=1:N %try the kgcb for N number of times

    %Plogy is the log values of KG for alternatives
    Plogy=[];
    Py=[];
    
    for iter1=1:M
        a=mu_est';
        b=covM(iter1,:)/sqrt(1/beta_W(iter1) + covM(iter1,iter1));
        
        [KG,LogKG] = LogEmaxAffine(a,b);
        
        [Plogy]=[Plogy, LogKG];
        [Py]=[Py,KG];
    end
    
    [maxh,x]=max(Plogy);
    
    %max_value is the best estimated value of the KG 
    %x is the argument that produces max_value

    %observe the outcome of the decision
    count(x)=count(x)+1;
    W_k=samples(x, count(x));

    e_x=zeros(M,1);
    e_x(x)=1;
    
    %updating equations for Normal-Normal model with covariance
    addscalar = (W_k - mu_est(x))/(1/beta_W(x) + covM(x,x));
    mu_est=mu_est + addscalar*covM*e_x;
    covM = covM - (covM*e_x*e_x'*covM)/((1/beta_W(x)) + covM(x,x));

end
    [value, recommendedArm] = max(mu_est);  
end




% logy = LogEmaxAffine(a,b)
% Calculates log(Exp[max_x a_x + b_x Z]-max_x a_x), where Z is a standard
% normal random variable and a,b are 1xM input vectors.
function [y,logy, a,b,c] = LogEmaxAffine(a,b)
    if (any(isnan(a)) || any(isnan(b)))
        warning('a or b is NaN');
    end
    assert(all(isreal(a)));
    assert(all(isreal(b)));

    a = a';
    b = b';

    % Check that a and b are column vectors of the right size
    if (any(size(a) ~= size(b)))
        error('LogEmaxAffine: a and b must be column vectors of the same size');
    end
    
    [a,b] = AffineBreakpointsPrep(a,b);
    
    [c, keep] = AffineBreakpoints(a,b); 
    a = a(keep);
    b = b(keep);
    c = c([1,keep+1]);
    M = length(keep);
    assert(all(isreal(c)));

    % I need logbdiff=log(diff(b)).  I thought that the following code would be
    % more numerically stable, able for example to distinguish cases like 
    % logb = [-25 -.3] vs. logb = [-35 -.3], but it doesn't seem to be able to.
    % Indeed, in the debugging output that I have below, the difference was 0.
    %{
    logb = log(abs(b)); % If b is 0, this is -Inf.
    sgnb = sign(b); % If b is 0, this is 0.
    logbdiff = zeros(size(c(2:M)));
    for i=1:length(b)-1
	[logbdiff(i),logbdiffsgn] = LogPlusExpSigned(logb(i),sgnb(i),logb(i+1),-sgnb(i+1));
	%assert(logbdiffsgn>=0);  % The b are distinct, so bdiff(i) can't be 0.
    end
    disp(sprintf('log(b)=%s log(diff(b))=%g logbdiff=%g difference=%g',mat2str(log(b)),log(diff(b)),logbdiff,log(diff(b))-logbdiff));
    %}
    logbdiff = log(diff(b))';  
    
    if M==1
        logy=log(a);
    elseif M>=2
        logy = LogSumExp(logbdiff+LogEI(-abs(c(2:M))));
    end

    logy=real(logy);
    y=exp(logy);
end
% Prepares vectors for passing to AffineEmaxBreakpoints, changing their
% order and removing elements with duplicate slope.

function [a,b] = AffineBreakpointsPrep(a,b)
    % Make sure a and b are column vectors.
    rows = size(a); if (rows == 1), a=a'; end
    rows = size(b); if (rows == 1), b=b'; end
    
    % 11/29/2008 PF: Experimental preprocessing step, which I hope will remove
    % a large number of the entries.
    [b1, i1] = min(b); % [a1,b1] is best at z=-infinity
    [a2, i2] = max(a); % [a2,b2] is best at z=0
    [b3, i3] = max(b); % [a3,b3] is best at z=+infinity
    a1 = a(i1);
    b2 = b(i2);
    a3 = a(i3);
    cleft = (a - a1)./(b1 - b); % intersection with leftmost line. 
    cright = (a - a3)./(b3 - b); % intersection with rightmost line.
    c2left = (a2 - a1)./(b1 - b2); % intersection with leftmost line. 
    c2right = (a2 - a3)./(b3 - b2); % intersection with rightmost line.
    keep = find(b==b1 | b==b3 | cleft <= c2left | cright >= c2right);
    %disp(sprintf('Preprocessing cut %d of %d entries', length(a)-length(keep), length(a)));
    a = a(keep);
    b = b(keep);
    clear keep cleft cright
   
        
    % Form a matrix for which ba(x,1) is the slope b(x) and ba(x,2) is the
    % y-intercept a(x).  Sort this matrix in ascending order of slope, 
    % breaking ties in slope with the y-intercept.  
    ba = [b, a];
    ba = sortrows(ba,[1,2]);
    a = ba(:,2);
    b = ba(:,1);
    
    % Then, from each pair of indices with the b component equal, remove
    % the one with smaller a component.  This code works because the sort
    % above enforced the condition: if b(i) == b(i+1), then a(i) <= a(i+1).
    keep = [find(diff(b)); length(b)];
    % This previous line is equivalent to:
    % keep = [];
    % for i=[1:length(b)-1]
    %    if b(i)~=b(i+1)
    %        keep = [keep, i];
    %    end
    %end 
    %keep = [keep, length(b)];  % We always keep the last one.
    
    % Note that the elements of keep are in ascending order.
    % This makes it so that b(keep) is still sorted in ascending order.
    a = a(keep);
    b = b(keep);
end

% Inputs are two M-vectors, a and b.
% Requires that the b vector is sorted in increasing order.
% Also requires that the elements of b all be unique.
% This function is used in AffineEmax, and the preparation of generic
% vectors a and b to satisfy the input requirements of this function are
% shown there.
%
% The output is an (M+1)-vector c and a vector A ("A" is for accept).  Think of
% A as a set which is a subset of {1,...,M}.  This output has the property
% that, for any i in {1,...,M} and any real number z,
%   i \in argmax_j a_j + b_j z
% iff
%   i \in A and z \in [c(j+1),c(i+1)],
%   where j = sup {0,1,...,i-1} \cap A.
%
% A note about indexing:
% Since Matlab does not allow indexing from 0, but instead requires
% indexing from 1, what is called c_i in the paper is written in matlab as
% c(1+i).  This is because in the paper we reference c_0.  For the vectors
% a and b, however, we don't need to reference a_0 or b_0, so we reference
% a_i and b_i by a(i) and b(i) respectively, rather than a(i+1) or b(i+1).
% 
function [c,A] = AffineBreakpoints(a,b)
    % Preallocate for speed.  Instead of resizing the array A whenever we add
    % to it or delete from it, we keep it the maximal size, and keep a length
    % indicator Alen telling us how many of its entries are good.  When the
    % function ends, we remove the unused elements from A before passing
    % it.
    M = length(a);
    c = zeros(1,M+1);
    A = zeros(1,M);
    
    % Step 0
    i=0;
    c(1+i) = -inf;
    c(1+i+1) = +inf;
    A(1) = 1;
    Alen = 1;
    
    for i=[1:M-1]
        c(1+i+1) = +inf;
        while(1)
            j = A(Alen); % jindex = Alen
            c(1+j) = (a(j) - a(i+1))/(b(i+1)-b(j));
	    % The if statement below replaces these lines from version 2 of the
	    % function.
	    %    kindex = jindex-1 = Alen-1
            %    if kindex > 0 && c(1+j)<=c(1+A(kindex))
            if Alen > 1 && c(1+j)<c(1+A(Alen-1))
		Alen = Alen-1; % Remove last element j
                % continue in while(1) loop
            else
                break % quit while(1) loop
            end
        end
	A(Alen+1) = i+1;
	Alen = Alen + 1;
    end
    A = A(1:Alen);
end

% Returns the log of Exp[(s+Z)^+], where s is a constant and Z is a standard
% normal random variable.  For large negative arguments Exp[(s+Z)^+] function
% is close to 0.  For large positive arguments, the function is close to the
% argument.  For s large enough, s>-10, we use the formula
% Exp[(s+Z)^+] = s*normcdf(s) + normpdf(s).  For smaller s we use an asymptotic
% approximation based on Mill's ratio.  EI stands for "expected improvement",
% since Exp[(s+Z)^+] would be the log of the expected improvement by measuring
% an alternative with excess predictive mean s over the best other measured
% alternative, and predictive variance 0.
function logy = LogEI(s)

% Use the asymptotic approximation for these large negative s.  The
% approximation is derived via:
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-|s|normcdf(-|s|)/normpdf(s)]
% and noting that normcdf(-|s|)/normpdf(s) is the Mill's ratio at |s|, which is
% asymptotically approximated by |s|/(s^2+1) [Gordon 1941, also documented in
% Frazier,Powell,Dayanik 2009 on page 14].  This gives,
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-s^2/(s^2+1)] = normpdf(s)/(s^2+1).

i=find(s<-10);
if (length(i)>0)
    logy(i) = LogNormPDF(s(i)) - log(s(i).^2 + 1);
end

% Use straightforward routines for s in the more numerically stable region.
i=find(s>=-10);
if (length(i)>0)
    logy(i) = log(s(i).*normcdf(s(i))+normpdf(s(i)));
end

assert(all(isreal(logy)));
end

% logy = LogNormPDF(z)
% Returns the log of the normal pdf evaluated at z.  z can be a vector or a scalar.
function logy = LogNormPDF(z)
	const = -.5*log(2*pi); % log of 1/sqrt(2pi).
	logy = const - z.^2/2;
end

% function y=LogSumExp(x)
% Computes log(sum(exp(x))) for a vector x, but in a numerically careful way.
function y=LogSumExp(x)
xmax = max(x);
y = xmax + log(sum(exp(x-xmax)));
end
