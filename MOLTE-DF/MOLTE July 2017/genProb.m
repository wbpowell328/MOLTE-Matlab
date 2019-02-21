%chose the problem class that one wants to generate the prob. of
%optimality/winning for
numberFolders=1
                                                                       
for i=1:numberFolders
    aa=load(['./', int2str(i),'/objectiveFunction.mat']);
    offline_avgcost=aa.offline_avgcost;
    %get the name of all policies
    policies=aa.policies;
    %chose the policies in comparison
    policy=[1,2,3,4,8,9]; 
    %calculate the mean and std of the final fit
    avg=mean(offline_avgcost(:,policy));
    stdd=std(offline_avgcost(:,policy));
    figure
    hold on
    cmap = colormap(cool);
    %chose color
    index=[1,8,24,40,56,64];
    
    %plot the oppurtunity cost the std of each policy
    for ii=1:length(policy)       
        bar(ii,avg(ii),'FaceColor',cmap(index(ii),:));
    end
    legend(policies{policy},'Location','NorthWest')
    errorbar(avg,stdd,'r.');
    ylabel('Opportunity cost','FontSize',18)
    set(gca,'xticklabel',[],'FontSize',15)
    hold off
%         print ('-painters', '-dpdf', '-r600', ['./',int2str(i), '/OC_hist.pdf']);
    print('-dtiff',['./',int2str(i),'/OC_hist.tif']);

    %calculate the prob. of not finding the actual best alternatives
    probError=sum(offline_avgcost(:,policy)==0)/length(offline_avgcost(:,1));
    
    %calculate the prob. of each policy being the best policy 
    probBest=zeros(1,length(policy));
    for j=1:length(policy)
        probBest(j)=sum(min(offline_avgcost(:,policy)')==offline_avgcost(:,policy(j))')/length(offline_avgcost(:,1));
    end
    figure
    colormap(cool)
    bar([probError;probBest]);
    legend( policies{policy},'Location','NorthWest');
    name={'Probability of optimality','Probability of winning'}
    set(gca,'xticklabel',name,'FontSize',15)
    
%          print ('-painters', '-dpdf', '-r600', ['./',int2str(i), '/Bar_hist.pdf']);
     print('-dtiff',['./',int2str(i),'/Bar_hist.tif']);
end
    