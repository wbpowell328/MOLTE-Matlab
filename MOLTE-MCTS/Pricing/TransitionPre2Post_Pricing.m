function currNode = TransitionPre2Post_Pricing(prevNode, currNode)
    currNode.theta = prevNode.theta;
    currNode.B = prevNode.B;
    currNode.C = prevNode.C;
    currNode.accR = prevNode.accR;
end