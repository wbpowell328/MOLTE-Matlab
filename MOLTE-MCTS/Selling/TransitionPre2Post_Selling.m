function currNode = TransitionPre2Post_Selling(prevNode, currNode, x)
    currNode.inv = prevNode.inv - x;
    currNode.acc_profit = prevNode.acc_profit;
end