function [currNode, Omega] = generateObservations_Selling(currNode, price, type, var)

if strcmp(type, 'D')
    Omega = price;
    currNode.pmf = containers.Map(price, 1);
else
    % This check isn't absolutely necessary
    if ~strcmp(type, 'S')
        fprintf('Invalid type %s! Using S (stochastic).', type);
    end
    Omega = unique(price + normrnd(0, var^2, 10, 1));
    currNode.pmf = containers.Map(Omega, ones(length(Omega),1)/length(Omega));
end

end