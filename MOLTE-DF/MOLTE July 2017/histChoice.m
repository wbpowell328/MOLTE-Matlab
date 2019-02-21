
%filedir is the location of the desired choice.mat
filedir='./1/';
numT=2; % draw the frenqucy of the chosen arms of the first 5 trails

a=load([filedir,'choice.mat']);

mkdir([filedir,'choices']);
count=a.choice;
policies=a.policies;
[a1,a2,a3]=size(count);
for j=1:a3
    for i=1:length(policies)
        bar(count(i,:,j), 'b','EdgeColor','b');
        axis tight 
        ylabel('number of times each alternative is chosen','FontSize',14);
        xlabel('number of alternatives','FontSize',14);
        print('-dtiff',[filedir,'choices/', int2str(j),'_',char(policies(i)),'.tif']);
      
    end
end

