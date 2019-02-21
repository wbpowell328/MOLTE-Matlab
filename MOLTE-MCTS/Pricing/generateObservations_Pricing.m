function [currNode, Omega] = generateObservations_Pricing(currNode, theta_truth, noise, x)
    Omega = unique([1 x x^2]*theta_truth + randn(1)*sqrt(noise));
    currNode.pmf = containers.Map(Omega, 1/length(Omega)*ones(length(Omega), 1));
end