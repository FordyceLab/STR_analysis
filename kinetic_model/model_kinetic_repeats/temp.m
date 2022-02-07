%%
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
koffData=nanmean(koffDataMaxSelected,2);
KDData=nanmean(KDDataMaxSelected,2);
%koffData=koffDataPho4Selected;
%KDData=KDDataPho4Selected;
weights=ones(numInData,1);
weightsAll=ones(numInData,1);
motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;
weights=zeros(numInData,1);
flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
for i=1:numInData
   if ismember(seqIdentifier(i,2),flanksToTrainOn)&&seqIdentifier(i,1)>0
       weights(i)=1;
   end
end
%weights(14)=1; %motif 1 with random flank, 14 18 22
%weights(3)=1;%motif 2
weights(motif1inds)=1;
weights(motif2inds)=1;%
weights(flankonlyinds)=1; %Don't fit on flank only;
%weights(26)=0;
%KDData(weights==1)=[1000;1500;2000];
%weights(seqIdentifier(:,1)==2)=0;%Dont fit on second motif
% p_core_rel p_flank_rel are used as parameters, where kon,micro,core =
% p_core_rel*koffmicro_flank

%weights=ones(numInData,1);
weights(sum(seqIdentifier,2)==0)=0; %Don't fit on the negative control seqences, that has no core or repat flank

numCores=regionsAdded(1);
numFlanks=regionsAdded(2);
paramIndex=[2 2+numCores 2+numCores+numFlanks 2+numCores+numFlanks+numCores]; %index of the start of [koffmicro_core koffmicro_flank p_core_rel p_flank_rel] index 1: konMax, index2: koff,M
numParams=1+2*regionsAdded(1)+2*regionsAdded(2);
numGuesses=10;

paramIndexEnd=paramIndex(2:end)-1;
paramIndexEnd(4)=numParams;

lb=0*ones(1,numParams);
ub=inf*ones(1,numParams);


lb(1)=0;%10^-5;%10^-2;



lb(paramIndex(3):paramIndexEnd(3))=0;%0.01; %p_core_rel 
lb(paramIndex(4):paramIndexEnd(4))=0;%0.01; %p_flank_rel

ub(1)=inf;%inf;%10^-5;%10^-2; %kon,max



%ub(paramIndex(1):paramIndexEnd(1))=10^2; %KDmicro_cores
%ub(paramIndex(2):paramIndexEnd(2))=10^2; %KDmicro_flank 
ub(paramIndex(3):paramIndexEnd(3))=inf; %p_core_rel 
ub(paramIndex(4):paramIndexEnd(4))=inf; %p_flank_rel


[params,residual]=fitCoreRepeatModelInfiniteTesting(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub);
[params2,residual2]=fitCoreRepeatModelInfiniteTesting(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub); % Run agian to see if same minima is found

[KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params2,weightsAll,seqIdentifier,paramIndex);
[~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);

figure();
plot(params(paramsFitted==1),'b.');
test=gca;
test.YScale='log';
xlabel('Used parameter index');
ylabel('params fit 1')
title(['residual 1: ' num2str(residual) 'residual 2: ' num2str(residual2)]);


figure();
plot(params2(paramsFitted==1)./params(paramsFitted==1),'b.');
test=gca;
test.YScale='log';
xlabel('Used parameter index');
ylabel('(params fit 2)/(params fit 1)')
title(['residual 1: ' num2str(residual) 'residual 2: ' num2str(residual2)]);
%[KDModelTemp,koffModelTemp] = getManyCoreRepeatModelKDmicro(params,1,seqIdentifierTemp,paramIndex);

f1=figure();
hold on;
title('Training data')
plot(KDModel(weights~=0),KDData(weights~=0,:),'b.');



xlabel('Model K_D (nM)');
ylabel('Experimental K_D (nM)');
xL=xlim();
yL=ylim();
%maxL=max([yL xL]);
%minL=min([yL xL]);
%xlim([minL maxL]);
%ylim([minL maxL]);

f2=figure();
hold on;
title('Training data')
plot(koffModel(weights~=0),koffData(weights~=0,:),'r.');
xlabel('Model koff (s^-^1)');
ylabel('Experimental koff (s^-^1)');


figure(f1);

plot(KDModel(motif1inds),KDData(motif1inds,:),'co');
plot(KDModel(motif2inds),KDData(motif2inds,:),'mo');
plot(KDModel(flankonlyinds),KDData(flankonlyinds,:),'kx');
xL=xlim();
yL=ylim();

maxL=max([yL xL]);
minL=min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);
plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
figure(f2);

plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
xL=xlim();
yL=ylim();
maxL=max([yL xL]);
minL=min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);
plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

