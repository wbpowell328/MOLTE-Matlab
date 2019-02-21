
% INPUT:
%spreedsheet; the file name of the spreedsheet that contains different 
%             problem classes

%This simulator will run each problem class for numP times
%In each run, it will run numTruth trails, e.g. numTruth (truth, prior,
%sample paths),computes the value of the objective function for each sample path and
%then averages the  numTruth  trials as the expected final reward or the expected cumulative rewards.

% run each problem class for numP times
numP=200;
%number of trails for each run  
numTruth=1;


spreadsheet='run.xls';
[ndata, text, alldata] = xlsread(spreadsheet,'','','basic');
[numPPp, numA]=size(ndata);
numPP=numPPp-1;


addpath('policies')
addpath('problemClasses')
addpath('auxillaryFunctions')
addpath('Prior')

matlabpool 8
for problemi=1:numPP
    
    %read in the characteristcs of all problem classes
    bm=char(text(problemi+1,5));
    if length(char(text(problemi+1,6)))==7
        offline=1;
    else
        offline=0;
    end
    NR=alldata{problemi+1, 4}; % ratio of number of alternatives
    Prior=alldata{problemi+1, 2}; % how can we get the prior
    problemClass=char(text(problemi+1,1));
    noiseRatio=alldata{problemi+1, 3}; % Noise level of problem classes
    policies=alldata(problemi+1,8:7+alldata{problemi+1,7});
    tuned=zeros(1,alldata{problemi+1,7}); %record whether each policy need to be tuned
    paraProvided=zeros(1,alldata{problemi+1,7});
    alpha=zeros(1,alldata{problemi+1,7}); % store the value of tunable parameters for each policy
    for jjj=1:alldata{problemi+1,7} % find out the policies name and whether they need to be tuned
        charp=char(policies(jjj));
        temppp=length(charp);
       if charp(temppp)==')'
           before=find(charp=='(');

          if charp(temppp-1)=='*'  
               tuned(jjj)=1;
          else
              alpha(jjj)=str2double(charp(before+1:temppp-1));
              paraProvided(jjj)=1;
          end

          policies{jjj}=charp(1:before-1);
       end

    end
    
    mkdir(int2str(problemi));                 
    paras={};
    % check whether the problem class have parameters
    if length(find( problemClass =='('))~=0 %#ok<ISMT> % this class has parameters 
         before=find(problemClass=='(');
         after=find(problemClass==')');
         paras=char(problemClass( before+1:after-1)); % parameters for GPR
         problemClass=problemClass(1:before-1);
    end

   %generate priors for each problem class
   [mu_0, covM, beta_W] = priorTruthGen(problemClass,noiseRatio, bm, paras, Prior); 

    M=length(mu_0);
    N=ceil(NR*M);
 %tune policies and store the tuned values in alpha.txt
    a=fopen(['./', int2str(problemi), '/alpha.txt'],'w');
    if offline==1
         fprintf(a,['offline', '\n']);
    else        
         fprintf(a,['online', '\n']);
    end      
    fclose(a);
    for i=1: length(policies)
        if tuned(i)==1   % need to tune the parameter
            if offline==1 % offline objective
                alpha(i)=offtune(char(policies(i)), problemClass, noiseRatio, M, N, bm,  paras, Prior);
            else   % online objective
                alpha(i)=ontune(char(policies(i)), problemClass, noiseRatio, M, N, bm, paras, Prior);
            end
            a=fopen(['./', int2str(problemi), '/alpha.txt'],'a');
            fprintf(a,[char(policies(i)), '=', num2str(alpha(i), 4), '\n']);
            fclose(a);
        end
        if(paraProvided(i)==1)
            tuned(i)=1;
        end
    end
    display('After tuning')
    if problemi ==1
        file_handle=fopen('latex.txt','w');
        fclose(file_handle);
    end

    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
    offline_avgcost=zeros(numP,length(policies));% opportunity cost, offline objective
    online_avgesum=zeros(numP,length(policies)); % expected totel reward, online objective
  
    parfor i=1:numP 
     %for i=1:numP       
        % for every problem class generate numTruth truths
        o_cost_off=zeros(numTruth,length(policies));
        mu_est=zeros(numTruth, M, length(policies));
        mu_ALL=zeros(numTruth, M);
        choices=zeros(length(policies), M, numTruth);
        ereward_sum=zeros(numTruth, length(policies));
        recommendedArm=zeros(numTruth, length(policies));

        for j=1:numTruth 
               %generate one possible truth      
               [mu,beta_W_p]=truthGen(problemClass,noiseRatio, mu_0, covM, beta_W);     
               mu_ALL(j,:)=mu;
               [samples ] = sampleGen( mu, beta_W_p, N );
               [ mu_est(j,:,:), choices(:,:,j), ereward_sum(j,:), recommendedArm(j,:) ] = policyRun(mu_0,beta_W_p,covM, samples, policies, mu, tuned, alpha);
               o_cost_off(j,:)=repmat(max(mu),length(policies), 1)-mu(recommendedArm(j,:)); %offline objective function, expected value of the 'best'
               o_cost_off(j,:)=o_cost_off(j,:)/range(mu);
               ereward_sum(j,:)=ereward_sum(j,:)/range(mu);
        end

        %evaluate policy using objective functions
        %offline OC:
        offline_avgcost(i,:)=sum(o_cost_off,1)/numTruth;
        
        %online average expected total rewards
        online_avgesum(i,:)=sum(ereward_sum,1)/numTruth/N;
        if i==1
          savefile1 = 'choice.mat';
          saveFile(savefile1, choices,policies, problemi) ; 
          savefile2='FinalFit';
          saveFile(savefile2, mu_est,mu_ALL, problemi) ; 
        end

    end


    %online comparison with the first policy
    KGesum=repmat(online_avgesum(:,1), 1,length(policies) );
    OC=(KGesum-online_avgesum);

    %offline comparison with the first policy
    KGesum1=repmat(offline_avgcost(:,1), 1,length(policies) );
    OC1=(offline_avgcost-KGesum1);

    savefile = 'objectiveFunction.mat';
    save(fullfile('./',int2str(problemi),savefile), 'offline_avgcost', 'online_avgesum','policies');


    color=[1,0,0;0,0.95,0.95;0.6,0.3,0.9;0,1,0;1,1,0;1,0.25,0.75;0,0.5,0.75;0.5,0.95,0.75;];
    
    %figure display
    firstPolicy=policies{1};
    remainingPolicies=policies(2:length(policies));
    if offline==1  %offline objective function, compared with the first policy
        figure
        hist(OC1(:,2:length(policies)),200);
        colormap (color);
        ylabel('number of experiments','FontSize',15);
        xlabel([firstPolicy,'-\pi'],'FontSize',15);
        legend(remainingPolicies);
        print ('-painters', '-dpdf', '-r600', ['./',int2str(problemi), '/offline_hist.pdf']);
        print('-dtiff',['./',int2str(problemi),'/offline_hist.tif']);

        % write numbers in latex table format
        file_handle=fopen('latex.txt','a');
        if strcmp(bm,'independent')
            bmha='I';
        else
            bmha='C';
        end

        fprintf(file_handle,[problemClass ' & ' bmha ' & ' '%.2f'],NR);
        for li=2:length(policies)
            temp=mean(OC1(:,li));
            temp3=std(OC1(:,li))/sqrt(numP);
            temp2=sum(OC1(:,li)<0)/length(OC1(:,li));
            fprintf(file_handle,[' & %.3f & %.3f & %.2f'],temp,temp3, temp2);
        end
        fprintf(file_handle,'\\\\ \n');
        fclose(file_handle);

    else            %online objective function,compared with the first policy
        figure
        hist(OC(:,2:length(policies)),200);
        colormap (color);
        ylabel('number of experiments','FontSize',15);
        xlabel([firstPolicy,'-\pi'],'FontSize',15);
        legend(remainingPolicies);
        print ('-painters', '-dpdf', '-r600', ['./',int2str(problemi), '/online_hist.pdf']);
        print('-dtiff',['./',int2str(problemi),'/online_hist.tif']);

        
        % write numbers in latex table format
        file_handle=fopen('latex.txt','a');
        if strcmp(bm,'independent')
            bmha='I';
        else
            bmha='C';
        end
        fprintf(file_handle,[problemClass ' & ' bmha ' & ' '%.2f'],NR);
        

        for li=2:length(policies)
            temp=mean(OC(:,li));
            temp3=std(OC(:,li))/sqrt(numP);
            temp2=sum(OC(:,li)<0)/length(OC(:,li));
            fprintf(file_handle,[' & %.3f & %.3f & %.2f'],temp,temp3, temp2);
        end
        fprintf(file_handle,'\\\\ \n');
        fclose(file_handle);
        
    end
    
 
    
end
matlabpool close;

