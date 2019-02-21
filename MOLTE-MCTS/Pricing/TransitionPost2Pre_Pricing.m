 function currNode = TransitionPost2Pre_Pricing(prevNode, currNode, w, var)
 
    x = [1; prevNode.action; prevNode.action^2];
    e       = w - prevNode.theta'*x;
    gamma   = 1 + x'*prevNode.B*x;
    
    currNode.theta   = prevNode.theta + 1/gamma*prevNode.B*x*e;
    currNode.B       = prevNode.B - 1/gamma*prevNode.B*(x*x')*prevNode.B;
    currNode.C       = currNode.B*var;
    currNode.accR    = prevNode.accR + w;
    
 end