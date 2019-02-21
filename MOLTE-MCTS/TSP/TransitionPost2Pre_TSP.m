function newNode = TransitionPost2Pre_TSP(oldNode, newNode, realizedCost)
newNode.visited = [oldNode.visited, oldNode.action];
newNode.city = oldNode.action;
newNode.accCost = oldNode.accCost + realizedCost;
end