
%filedir is the location of the desired choice.mat
filedir='./1/'
numT=2; % draw the final fit of the first numT trails

a=load([filedir,'FinalFit.mat']);

mkdir([filedir,'FinalFit']);
mu_est=a.choice;
mu_ALL=a.policies;
aa=load([filedir,'choice.mat']);
policies=aa.policies;
color=[1,0,0;0,0.95,0.95;0.6,0.3,0.9;0,1,0;1,1,0;1,0.25,0.75;0,0.5,0.75;0.5,0.95,0.75;];
for j=1:numT

    for i=1:length(policies)
      
        plot(mu_ALL(j,:));
        hold on
        plot(mu_est(j,:,i),'r');
        ylabel('number of times each alternative is chosen','FontSize',14);
        xlabel('number of alternatives','FontSize',14);   
        axis tight  
        hold off
        print('-dtiff',[filedir,'FinalFit/',int2str(j),'_',char(policies(i)),'.tif']);
    end

end

