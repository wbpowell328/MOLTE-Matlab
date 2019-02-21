function grad = vanillanewsvendorgrad(demand, x, price, cost)
% returns gradient for the vanilla newsvendor problem
% notice that the demand belongs to some probability distribution
% F = price*min(demand, x) - cx 


 if (x <= demand) 
     gradF = price-cost;  
 else
     gradF = -cost; 
 end 

 grad = gradF;
end
