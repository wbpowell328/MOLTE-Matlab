function currNode = TransitionPost2Pre_Selling(prevNode, currNode, w)
%     Tree(new_ind).inv = max([Tree(curr_ind).inv - w 0]);
    currNode.inv = prevNode.inv;
%     Tree(new_ind).V_x = Tree(curr_ind).V_x + w*Tree(curr_ind).action;
    if currNode.acc_profit == 0
        currNode.acc_profit = prevNode.acc_profit + w*prevNode.action;
    end
end