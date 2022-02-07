% This script analyses kinetic measurement data of motif and repeats DNA
% binding Max and Pho4 from Connor
%impleementing constrains in fitting, so that all fitted p_rel and
%koffmicro are only allwoed be withing some factor from each oterh
%211012_normcscore:Do'nt usen oramlized squared error, normalzied to each
%KDData, normalize to one valie for KD and one for koff
addpath(genpath('Functions'));


%% Hej
seqPath='../data/library1.csv';
dataPath='../data/library1_Kd_koff.csv';
concensusMotif='GTCACGTGAC';


%Load data
seqAndNames=importdata(seqPath);
numInSeq=numel(seqAndNames);
for i=1:numInSeq
    seqAndNames{i}=strsplit(seqAndNames{i},',');
end

data=importdata(dataPath,',');
%% Set parameters
startIndexes=1:3;
flank1Indexes=4:16;
motifBaseIndexes=17:26;
flank2Indexes=27:39;
endIndexes=40:54;

numRegions=2; %core and flanks, two states

setNoMotifZero=1;%If binding probabilies for no motif is to be fixed at 0
setRandFlankZero=1;%If binding probabilies for random flanks is to be fixed at 0
lockFlankControls=1;%Since flanks + no motif are not identical to the sequences with motif, lock the flank control sequences to the same parameter as for the most simlar flank+motif

%Get the sequence and name of each sequence, with the same indexing
numInData=size(data.data,1);
sequences=cell(1,numInData);
newNames=cell(1,numInData);
motifNames=cell(1,numInData);
flankNames=cell(1,numInData);

for i=1:numInSeq
    for j=1:numInData
        if strcmp(data.textdata{j+1,2},seqAndNames{i}{2})
            sequences{j}=seqAndNames{i}{3};
            newNames{j}=seqAndNames{i}{2};
            thestr=strsplit(newNames{j},' + ');
            motifNames{j}=thestr{1};
            flankNames{j}=thestr{2};
            break;
        end
    end
end
koffDataMaxAll=data.data(:,1:4:12);
KDDataMaxAll=data.data(:,3:4:12); % I think units are nM
koffDataPho4All=data.data(:,13:4:end);
KDDataPho4All=data.data(:,15:4:end);
sequencesOrig=sequences;

%Chance the flank only sequences so they mathces the cloesest mathces for a
%sequence that has the motif, done manually
if lockFlankControls
    sequences{26}(flank1Indexes)=sequences{17}(flank1Indexes);
    sequences{26}(flank2Indexes)=sequences{17}(flank2Indexes);
    
    sequences{29}(flank1Indexes)=sequences{15}(flank1Indexes);
    sequences{29}(flank2Indexes)=sequences{15}(flank2Indexes);
    
end


%Check if each sequence has consencus motif, and if flanks are reverse
%complement of each other

flankRevCompl=zeros(1,numInData);
hasConcensusMotif=zeros(1,numInData);

for j=1:numInData
  %  if strcmp(sequences{j}(flank1Indexes),seqrcomplement(sequences{j}(flank2Indexes)))
  %      flankRevCompl(j)=1;
  %      
  %  end
    
    if strcmp(sequences{j}(motifBaseIndexes),concensusMotif)
        hasConcensusMotif(j)=1;
    end
end

% Create unique identified for each core and flanks

seqIdentifier=zeros(numInData,numRegions);
uniqueMotifs={};
uniqueFlank1={};
uniqueFlank2={};
motifHasBindProb=[];
flankHasBindProb=[];
regionsAdded=zeros(1,numRegions);
numMotifs=0;
numFlanks=0;
for i=1:numInData
    %Check for this core motif
    if (setNoMotifZero==1&&~strcmp(motifNames{i},'No motif'))||setNoMotifZero==0
        for j=1:regionsAdded(1)
            if strcmp(sequences{i}(motifBaseIndexes),uniqueMotifs{j})
                seqIdentifier(i,1)=j;
                break;
            end
        end
        if seqIdentifier(i,1)==0 %First time seing this sequence
            regionsAdded(1)=regionsAdded(1)+1;
            uniqueMotifs{regionsAdded(1)}=sequences{i}(motifBaseIndexes);
            seqIdentifier(i,1)=regionsAdded(1);
        end
    end
    %Check for this flanks
    if (setNoMotifZero==1&&~strcmp(flankNames{i}(1:6),'random'))||setRandFlankZero==0
        for j=1:regionsAdded(2)
            if strcmp(sequences{i}(flank1Indexes),uniqueFlank1{j})&&strcmp(sequences{i}(flank2Indexes),uniqueFlank2{j})
                seqIdentifier(i,2)=j;
                break;
            end
        end
        if seqIdentifier(i,2)==0 %First time seing this sequence
            regionsAdded(2)=regionsAdded(2)+1;
            uniqueFlank1{regionsAdded(2)}=sequences{i}(flank1Indexes);
            uniqueFlank2{regionsAdded(2)}=sequences{i}(flank2Indexes);
            seqIdentifier(i,2)=regionsAdded(2);
        end
    end
    
end
koffDataMaxSelected=koffDataMaxAll;
KDDataMaxSelected=KDDataMaxAll; % I think units are nM
koffDataPho4Selected=koffDataPho4All;
KDDataPho4Selected=KDDataPho4All;

%trainingDataset=1:numInData; %Maybe only train on motif and mutatted motifs and the flanks who exists for both, 2*2*7 datapoints, 1+2*2+2*7 variables


%% testing theoretical experssion of mean first passage times/ eigenvalues
times=[0:0.005:800];
%times=[0:0.00005:1];
params=[1 10^6 0.05 10 0.1 0.4];% KD = koffSlope/konSlope but not equal %%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
%params=[1 10^6 0.05 0.02 0.1 0.4];
%params=[1 10^6 0.05 100 0.1 0.4];
%params=[1 10^6 0.05 0.02 0.4 0.4];
%koffMeanFirstPassage/konMeanFirstPassge
%params=[1 0.6*10^6 0.05 10 0.1/0.6 0];%KD = koffSlope/konSlope =koffMeanFirstPassage/konMeanFirstPassge
%params=[1 10^6 0.05 0.1 0.1 0.4];
%params=[1 10^6 0.05 100 0.05 0.9];
%params=[1 0.1*10^6 0.05 100 0.5 0];
solutionConc=1;
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient
slopePoints=5;

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

funcValsOn=solveODEs(Aon,bon,x0on,times);
funcValsOff=solveODEs(Aoff,boff,x0off,times);
%funcValsOff=solveODEs(Aoff,boff,[0 params(5)/(params(5)+params(6)) params(6)/(params(5)+params(6))]',times);

onCurve=sum(funcValsOn(2:3,:));

theSum=sum(x0off);

KD = solutionConc*(1-theSum)/theSum;


%survivalFuncOff=sum(funcValsOff)/theSum;
survivalFuncOff=sum(funcValsOff(2:3,:));
survivalFuncOff=survivalFuncOff/survivalFuncOff(1);
%Formula for expecation value using the cdf (1-survival function)
offTime=trapz(times,survivalFuncOff);
koffIntegral=1/offTime;

p=polyfit(times(1:slopePoints),survivalFuncOff(1:slopePoints),1);
koffSlope=-p(1);
p=polyfit(times(1:slopePoints),onCurve(1:slopePoints),1);
konSlope=p(1)/solutionConc;

%Try with mean first passagetimes

%tiem from dissocatied to core OR flank

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage=1/offMeanFirstPassage;

koffMeanFirstPassageInitial=1/offMeanFirstPassageInitial;

figure();
Eqn='exp(-b*x)';
options = fitoptions(Eqn);
options.Lower=0;
options.Upper=10^6;
thefit=fit(times',survivalFuncOff',Eqn,options);
plot(thefit,times,survivalFuncOff);
%% Single fit, no bootstrapping
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
koffData=nanmean(koffDataPho4Selected,2);
KDData=nanmean(KDDataPho4Selected,2);
%koffData=koffDataPho4Selected;
%KDData=KDDataPho4Selected;
weights=ones(numInData,1);
weightsAll=ones(numInData,1);
motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;



weights=zeros(numInData,1);
flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
combToTrainOn=zeros(numInData,1);
for i=1:numInData
   if ismember(seqIdentifier(i,2),flanksToTrainOn)&&seqIdentifier(i,1)>0
       weights(i)=1;
       combToTrainOn(i)=1;
   end
end
motif1andflankinds=(seqIdentifier(:,1)==1)+combToTrainOn==2;
motif2andflankinds=(seqIdentifier(:,1)==2)+combToTrainOn==2;

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
numGuesses=1000;
koffmicrofac=10;
prelfac=10;

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


[params,residual]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
[params2,residual2]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac); % Run agian to see if same minima is found

[KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params2,weightsAll,seqIdentifier,paramIndex);
[~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);

singleFit.params=params;
singleFit.param2=params2;
singleFit.residual=residual;
singleFit.residual2=residual2;
singleFit.KDModel=KDModel;
singleFit.koffModel=koffModel;



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
plot(KDModel(motif1andflankinds),KDData(motif1andflankinds,:),'cd');
plot(KDModel(motif2andflankinds),KDData(motif2andflankinds,:),'md');
test=gca;
test.XScale='log';
test.YScale='log';

xL=xlim();
yL=ylim();

maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);


plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
savefig('KD_data_vs_model');

figure(f2);

plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
plot(koffModel(motif1andflankinds),koffData(motif1andflankinds,:),'cd');
plot(koffModel(motif2andflankinds),koffData(motif2andflankinds,:),'md');
xL=xlim();
yL=ylim();
maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);
plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

savefig('koff_data_vs_model')
%% Bootstrapping
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
%koffData=nanmean(koffDataMaxSelected,2);
%KDData=nanmean(KDDataMaxSelected,2);
koffmicrofac=10;
prelfac=10;
numReps=200;%size(koffDataMaxSelected,2);
numFits=2;
numGuesses=500;
reps=cell(1,numReps);

paramsReps=[];
koffModelReps=[];
KDModelReps=[];

koffDataAll=koffDataPho4Selected;
KDDataAll=KDDataPho4Selected;


koffError=zeros(numInData,1);
KDError=zeros(numInData,1);

for i=1:numInData
   koffError(i)=nanstd(koffDataAll(i,:))/sqrt(sum(~isnan(koffDataAll(i,:)))); 
   KDError(i)=nanstd(KDDataAll(i,:))/sqrt(sum(~isnan(KDDataAll(i,:)))); 
    
end

f1=figure();
f2=figure();
f3=figure();
f4=figure();

f5=figure();
for i=1:numReps
    disp(i);
    koffData=nanmean(koffDataAll,2)+randn(numInData,1).*koffError;
    KDData=nanmean(KDDataAll,2)+randn(numInData,1).*KDError;
    
    
    
    weights=ones(numInData,1);
    weightsAll=ones(numInData,1);
    motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
    motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
    flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;
    weights=zeros(numInData,1);
    flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
    for k=1:numInData
        if ismember(seqIdentifier(k,2),flanksToTrainOn)&&seqIdentifier(k,1)>0
            weights(k)=1;
        end
    end
    %weights(14)=1; %motif 1 with random flank, 14 18 22
    %weights(3)=1;%motif 2
    weights(motif1inds)=1;
    weights(motif2inds)=1;%
    weights(flankonlyinds)=1; %Don't fit on flank only;
    %weights(29)=0;
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
    
    
    paramIndexEnd=paramIndex(2:end)-1;
    paramIndexEnd(4)=numParams;
    
    lb=0*ones(1,numParams);
    ub=inf*ones(1,numParams);
    
    
    lb(1)=0;%0;%10^-5;%10^-2;
    
    
    
    lb(paramIndex(3):paramIndexEnd(3))=0;%0.01; %p_core_rel
    lb(paramIndex(4):paramIndexEnd(4))=0;%0.01; %p_flank_rel
    
    ub(1)=inf;%inf;%inf;%10^-5;%10^-2; %kon,max
    
    
    
    %ub(paramIndex(1):paramIndexEnd(1))=10^2; %KDmicro_cores
    %ub(paramIndex(2):paramIndexEnd(2))=10^2; %KDmicro_flank
    ub(paramIndex(3):paramIndexEnd(3))=inf; %p_core_rel
    ub(paramIndex(4):paramIndexEnd(4))=inf; %p_flank_rel
    paramsVec=[];
    residualVec=[];
    paramsAllVec=cell(1,numFits);
    residualAllVec=[];
    KDModelVec=[];
    koffModelVec=[];
    for j=1:numFits
        disp(['Fit num: ' num2str(j)]);
        [paramsCurr,residualCurr,paramsAllCurr,residualAllCurr]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
        paramsVec(j,:)=paramsCurr;
        residualVec(j,:)=residualCurr;
        paramsAllVec{j}=paramsAllCurr;
        residualAllVec(j,:)=residualAllCurr;
        params=paramsVec(j,:);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        KDModelVec(j,:)=KDModel';
        koffModelVec(j,:)=koffModel';
    end
    params=paramsVec(1,:);
    
    [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
    [~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);
    
% % % % % %     figure(f1);
% % % % % %     hold off;
% % % % % %     plot(params(paramsFitted==1),'b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('params fit 1')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     
% % % % % %     
% % % % % %     figure(f2);
% % % % % %     hold off;
% % % % % %     plot((paramsVec(:,paramsFitted==1)./repmat(params(paramsFitted==1),[numFits 1]))','b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('(params fit 2)/(params fit 1)')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     %[KDModelTemp,koffModelTemp] = getManyCoreRepeatModelKDmicro(params,1,seqIdentifierTemp,paramIndex);
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(KDModel(weights~=0),KDData(weights~=0,:),'b.');
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     xlabel('Model K_D (nM)');
% % % % % %     ylabel('Experimental K_D (nM)');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     %maxL=max([yL xL]);
% % % % % %     %minL=min([yL xL]);
% % % % % %     %xlim([minL maxL]);
% % % % % %     %ylim([minL maxL]);
% % % % % %     
% % % % % %     figure(f4);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(koffModel(weights~=0),koffData(weights~=0,:),'r.');
% % % % % %     xlabel('Model koff (s^-^1)');
% % % % % %     ylabel('Experimental koff (s^-^1)');
% % % % % %     
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold on;
% % % % % %     plot(KDModel(motif1inds),KDData(motif1inds,:),'co');
% % % % % %     plot(KDModel(motif2inds),KDData(motif2inds,:),'mo');
% % % % % %     plot(KDModel(flankonlyinds),KDData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
% % % % % %     figure(f4);
% % % % % %     hold on;
% % % % % %     plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
% % % % % %     plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
% % % % % %     plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
    
    reps{i}.paramsVec=paramsVec;
    reps{i}.residualVec=residualVec;
    reps{i}.paramsAllVec=paramsAllVec;
    reps{i}.residualAllVec=residualAllVec;
    reps{i}.KDModelVec=KDModelVec;
    reps{i}.koffModelVec=koffModelVec;
    
    paramsReps(i,:)=params;
    KDModelReps(:,i)=KDModel;
    koffModelReps(:,i)=koffModel;
    
    figure(f5);
    nCompleted=i;
    hold off;
    plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
    hold on;
    test=gca;
    test.YScale='log';
    xlabel('paramIndex');
    ylabel('paramValue')
    title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    
end
nCompleted=numReps;
figure();
plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
hold on;
test=gca;
test.YScale='log';
xlabel('paramIndex');
ylabel('paramValue')
title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    


figure();
errorbar(mean(log10(paramsReps(1:nCompleted,paramsFitted==1))),std(log10(paramsReps(1:nCompleted,paramsFitted==1))),'bx');

%% Combined plots single and boot

perc=0.68;

paramsCenters=singleFit.param2(paramsFitted==1);
paramsBoots=paramsReps(1:nCompleted,paramsFitted==1);
paramsCentersBoot=mean(paramsBoots,1);

numCurrParams=numel(paramsCenters);

startInd=round(((1-perc)/2)*numReps)+1;

endInd=numReps-startInd+1;

upperBounds=zeros(1,numCurrParams);
lowerBounds=zeros(1,numCurrParams);

for i=1:numCurrParams
   sortedVals=sort(paramsBoots(:,i),'ascend');
   upperBounds(i)=sortedVals(endInd);
   lowerBounds(i)=sortedVals(startInd);
end

figure();
errorbar(1:numCurrParams,paramsCenters,paramsCentersBoot-lowerBounds,upperBounds-paramsCentersBoot,'ro');
test=gca;
test.YScale='log';
xlabel('Param index');
ylabel('Param value');
title(['confidence interval : ' num2str(perc)]);
%savefig('fitted_params')

%% Single fit, no bootstrapping, dont fit on motif1 only data
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
koffData=nanmean(koffDataPho4Selected,2);
KDData=nanmean(KDDataPho4Selected,2);
%koffData=koffDataPho4Selected;
%KDData=KDDataPho4Selected;
weights=ones(numInData,1);
weightsAll=ones(numInData,1);
motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;



weights=zeros(numInData,1);
flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
combToTrainOn=zeros(numInData,1);
for i=1:numInData
   if ismember(seqIdentifier(i,2),flanksToTrainOn)&&seqIdentifier(i,1)>0
       weights(i)=1;
       combToTrainOn(i)=1;
   end
end
motif1andflankinds=(seqIdentifier(:,1)==1)+combToTrainOn==2;
motif2andflankinds=(seqIdentifier(:,1)==2)+combToTrainOn==2;

%weights(14)=1; %motif 1 with random flank, 14 18 22
%weights(3)=1;%motif 2
weights(motif1inds)=0;
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
numGuesses=1000;
koffmicrofac=10;
prelfac=10;

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


[params,residual]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
[params2,residual2]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac); % Run agian to see if same minima is found

[KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params2,weightsAll,seqIdentifier,paramIndex);
[~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);

singleFitTest1.params=params;
singleFitTest1.param2=params2;
singleFitTest1.residual=residual;
singleFitTest1.residual2=residual2;
singleFitTest1.KDModel=KDModel;
singleFitTest1.koffModel=koffModel;



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

%savefig('KD_data_vs_model')

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
plot(KDModel(motif1andflankinds),KDData(motif1andflankinds,:),'cd');
plot(KDModel(motif2andflankinds),KDData(motif2andflankinds,:),'md');
test=gca;
test.XScale='log';
test.YScale='log';

xL=xlim();
yL=ylim();

maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);


plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

figure(f2);

plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
plot(koffModel(motif1andflankinds),koffData(motif1andflankinds,:),'cd');
plot(koffModel(motif2andflankinds),koffData(motif2andflankinds,:),'md');
xL=xlim();
yL=ylim();
maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);
plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

%savefig('koff_data_vs_model')

figure();
hold on;
bar([mean(KDData(motif1inds)) mean(KDModel(motif1inds))]);
ylabel('KD (nM)')
xticks([1 2]);
test=gca;
test.XTickLabel{1}='motif1 data';
test.XTickLabel{2}='motif1 model';

savefig('testing_KD_motif1');

figure();
hold on;
bar([mean(koffData(motif1inds)) mean(koffModel(motif1inds))]);
ylabel('koff (s^-^1)')
xticks([1 2]);
test=gca;
test.XTickLabel{1}='motif1 data';
test.XTickLabel{2}='motif1 model';
savefig('testing_koff_motif1');

%% Bootstrapping, smaller training dataset, don't fit on motif1 only
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
%koffData=nanmean(koffDataMaxSelected,2);
%KDData=nanmean(KDDataMaxSelected,2);
koffmicrofac=10;
prelfac=10;
numReps=200;%size(koffDataMaxSelected,2);
numFits=2;
numGuesses=500;
repsTest1=cell(1,numReps);

paramsRepsTest1=[];
koffModelRepsTest1=[];
KDModelRepsTest1=[];

koffDataAll=koffDataPho4Selected;
KDDataAll=KDDataPho4Selected;


koffError=zeros(numInData,1);
KDError=zeros(numInData,1);

for i=1:numInData
   koffError(i)=nanstd(koffDataAll(i,:))/sqrt(sum(~isnan(koffDataAll(i,:)))); 
   KDError(i)=nanstd(KDDataAll(i,:))/sqrt(sum(~isnan(KDDataAll(i,:)))); 
    
end

f1=figure();
f2=figure();
f3=figure();
f4=figure();

f5=figure();
for i=1:numReps
    disp(i);
    koffData=nanmean(koffDataAll,2)+randn(numInData,1).*koffError;
    KDData=nanmean(KDDataAll,2)+randn(numInData,1).*KDError;
    
    
    
    weights=ones(numInData,1);
    weightsAll=ones(numInData,1);
    motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
    motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
    flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;
    weights=zeros(numInData,1);
    flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
    for k=1:numInData
        if ismember(seqIdentifier(k,2),flanksToTrainOn)&&seqIdentifier(k,1)>0
            weights(k)=1;
        end
    end
    %weights(14)=1; %motif 1 with random flank, 14 18 22
    %weights(3)=1;%motif 2
    weights(motif1inds)=0;
    weights(motif2inds)=1;%
    weights(flankonlyinds)=1; %Don't fit on flank only;
    %weights(29)=0;
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
    
    
    paramIndexEnd=paramIndex(2:end)-1;
    paramIndexEnd(4)=numParams;
    
    lb=0*ones(1,numParams);
    ub=inf*ones(1,numParams);
    
    
    lb(1)=0;%0;%10^-5;%10^-2;
    
    
    
    lb(paramIndex(3):paramIndexEnd(3))=0;%0.01; %p_core_rel
    lb(paramIndex(4):paramIndexEnd(4))=0;%0.01; %p_flank_rel
    
    ub(1)=inf;%inf;%inf;%10^-5;%10^-2; %kon,max
    
    
    
    %ub(paramIndex(1):paramIndexEnd(1))=10^2; %KDmicro_cores
    %ub(paramIndex(2):paramIndexEnd(2))=10^2; %KDmicro_flank
    ub(paramIndex(3):paramIndexEnd(3))=inf; %p_core_rel
    ub(paramIndex(4):paramIndexEnd(4))=inf; %p_flank_rel
    paramsVec=[];
    residualVec=[];
    paramsAllVec=cell(1,numFits);
    residualAllVec=[];
    KDModelVec=[];
    koffModelVec=[];
    for j=1:numFits
        disp(['Fit num: ' num2str(j)]);
        [paramsCurr,residualCurr,paramsAllCurr,residualAllCurr]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
        paramsVec(j,:)=paramsCurr;
        residualVec(j,:)=residualCurr;
        paramsAllVec{j}=paramsAllCurr;
        residualAllVec(j,:)=residualAllCurr;
        params=paramsVec(j,:);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        KDModelVec(j,:)=KDModel';
        koffModelVec(j,:)=koffModel';
    end
    params=paramsVec(1,:);
    
    [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
    [~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);
    
% % % % % %     figure(f1);
% % % % % %     hold off;
% % % % % %     plot(params(paramsFitted==1),'b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('params fit 1')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     
% % % % % %     
% % % % % %     figure(f2);
% % % % % %     hold off;
% % % % % %     plot((paramsVec(:,paramsFitted==1)./repmat(params(paramsFitted==1),[numFits 1]))','b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('(params fit 2)/(params fit 1)')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     %[KDModelTemp,koffModelTemp] = getManyCoreRepeatModelKDmicro(params,1,seqIdentifierTemp,paramIndex);
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(KDModel(weights~=0),KDData(weights~=0,:),'b.');
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     xlabel('Model K_D (nM)');
% % % % % %     ylabel('Experimental K_D (nM)');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     %maxL=max([yL xL]);
% % % % % %     %minL=min([yL xL]);
% % % % % %     %xlim([minL maxL]);
% % % % % %     %ylim([minL maxL]);
% % % % % %     
% % % % % %     figure(f4);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(koffModel(weights~=0),koffData(weights~=0,:),'r.');
% % % % % %     xlabel('Model koff (s^-^1)');
% % % % % %     ylabel('Experimental koff (s^-^1)');
% % % % % %     
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold on;
% % % % % %     plot(KDModel(motif1inds),KDData(motif1inds,:),'co');
% % % % % %     plot(KDModel(motif2inds),KDData(motif2inds,:),'mo');
% % % % % %     plot(KDModel(flankonlyinds),KDData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
% % % % % %     figure(f4);
% % % % % %     hold on;
% % % % % %     plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
% % % % % %     plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
% % % % % %     plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
    
    repsTest1{i}.paramsVec=paramsVec;
    repsTest1{i}.residualVec=residualVec;
    repsTest1{i}.paramsAllVec=paramsAllVec;
    repsTest1{i}.residualAllVec=residualAllVec;
    repsTest1{i}.KDModelVec=KDModelVec;
    repsTest1{i}.koffModelVec=koffModelVec;
    
    paramsRepsTest1(i,:)=params;
    KDModelRepsTest1(:,i)=KDModel;
    koffModelRepsTest1(:,i)=koffModel;
    
    figure(f5);
    nCompleted=i;
    hold off;
    plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
    hold on;
    test=gca;
    test.YScale='log';
    xlabel('paramIndex');
    ylabel('paramValue')
    title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    
end
nCompleted=numReps;
figure();
plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
hold on;
test=gca;
test.YScale='log';
xlabel('paramIndex');
ylabel('paramValue')
title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    


figure();
errorbar(mean(log10(paramsReps(1:nCompleted,paramsFitted==1))),std(log10(paramsReps(1:nCompleted,paramsFitted==1))),'bx');

%% Combined plots single and boot

perc=0.68;

paramsCenters=singleFit.params(paramsFitted==1);
paramsBoots=paramsReps(1:nCompleted,paramsFitted==1);

numCurrParams=numel(paramsCenters);

startInd=round(((1-perc)/2)*numReps)+1;

endInd=numReps-startInd+1;

upperBounds=zeros(1,numCurrParams);
lowerBounds=zeros(1,numCurrParams);

for i=1:numCurrParams
   sortedVals=sort(paramsBoots(:,i),'ascend');
   upperBounds(i)=sortedVals(endInd);
   lowerBounds(i)=sortedVals(startInd);
end

figure();
errorbar(1:numCurrParams,paramsCenters,paramsCenters-lowerBounds,upperBounds-paramsCenters,'ro');
test=gca;
test.YScale='log';
xlabel('Param index');
ylabel('Param value');
title(['confidence interval : ' num2str(perc)]);
savefig('fitted_params_lesstraining')


%% Single fit, no bootstrapping, dont fit on motif2 only data
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
koffData=nanmean(koffDataPho4Selected,2);
KDData=nanmean(KDDataPho4Selected,2);
%koffData=koffDataPho4Selected;
%KDData=KDDataPho4Selected;
weights=ones(numInData,1);
weightsAll=ones(numInData,1);
motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;



weights=zeros(numInData,1);
flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
combToTrainOn=zeros(numInData,1);
for i=1:numInData
   if ismember(seqIdentifier(i,2),flanksToTrainOn)&&seqIdentifier(i,1)>0
       weights(i)=1;
       combToTrainOn(i)=1;
   end
end
motif1andflankinds=(seqIdentifier(:,1)==1)+combToTrainOn==2;
motif2andflankinds=(seqIdentifier(:,1)==2)+combToTrainOn==2;

%weights(14)=1; %motif 1 with random flank, 14 18 22
%weights(3)=1;%motif 2
weights(motif1inds)=1;
weights(motif2inds)=0;%
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
numGuesses=1000;
koffmicrofac=10;
prelfac=10;

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


[params,residual]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
[params2,residual2]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac); % Run agian to see if same minima is found

[KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params2,weightsAll,seqIdentifier,paramIndex);
[~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);

singleFitTest2.params=params;
singleFitTest2.param2=params2;
singleFitTest2.residual=residual;
singleFitTest2.residual2=residual2;
singleFitTest2.KDModel=KDModel;
singleFitTest2.koffModel=koffModel;



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

%savefig('KD_data_vs_model')

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
plot(KDModel(motif1andflankinds),KDData(motif1andflankinds,:),'cd');
plot(KDModel(motif2andflankinds),KDData(motif2andflankinds,:),'md');
test=gca;
test.XScale='log';
test.YScale='log';

xL=xlim();
yL=ylim();

maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);


plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

figure(f2);

plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
plot(koffModel(motif1andflankinds),koffData(motif1andflankinds,:),'cd');
plot(koffModel(motif2andflankinds),koffData(motif2andflankinds,:),'md');
xL=xlim();
yL=ylim();
maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);
plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

%savefig('koff_data_vs_model')

figure();
hold on;
bar([mean(KDData(motif2inds)) mean(KDModel(motif2inds))]);
ylabel('KD (nM)')
xticks([1 2]);
test=gca;
test.XTickLabel{1}='motif2 data';
test.XTickLabel{2}='motif2 model';

savefig('testing_KD_motif2');

figure();
hold on;
bar([mean(koffData(motif2inds)) mean(koffModel(motif2inds))]);
ylabel('koff (s^-^1)')
xticks([1 2]);
test=gca;
test.XTickLabel{1}='motif2 data';
test.XTickLabel{2}='motif2 model';
savefig('testing_koff_motif2');

%% Bootstrapping, smaller training dataset, don't fit on motif2 only
%koffData=nanmean(koffDataPho4Selected,2);
%koffData(:)=NaN;
%KDData=nanmean(KDDataPho4Selected,2);
%koffData=nanmean(koffDataMaxSelected,2);
%KDData=nanmean(KDDataMaxSelected,2);
koffmicrofac=10;
prelfac=10;
numReps=200;%size(koffDataMaxSelected,2);
numFits=2;
numGuesses=500;
repsTest2=cell(1,numReps);

paramsRepsTest2=[];
koffModelRepsTest2=[];
KDModelRepsTest2=[];

koffDataAll=koffDataPho4Selected;
KDDataAll=KDDataPho4Selected;


koffError=zeros(numInData,1);
KDError=zeros(numInData,1);

for i=1:numInData
   koffError(i)=nanstd(koffDataAll(i,:))/sqrt(sum(~isnan(koffDataAll(i,:)))); 
   KDError(i)=nanstd(KDDataAll(i,:))/sqrt(sum(~isnan(KDDataAll(i,:)))); 
    
end

f1=figure();
f2=figure();
f3=figure();
f4=figure();

f5=figure();
for i=1:numReps
    disp(i);
    koffData=nanmean(koffDataAll,2)+randn(numInData,1).*koffError;
    KDData=nanmean(KDDataAll,2)+randn(numInData,1).*KDError;
    
    
    
    weights=ones(numInData,1);
    weightsAll=ones(numInData,1);
    motif1inds=(seqIdentifier(:,1)==1)+(seqIdentifier(:,2)==0)==2;
    motif2inds=(seqIdentifier(:,1)==2)+(seqIdentifier(:,2)==0)==2;
    flankonlyinds=(seqIdentifier(:,1)==0)+(seqIdentifier(:,2)~=0)==2;
    weights=zeros(numInData,1);
    flanksToTrainOn=[3 8];%1:16;%[1 3 9 11];
    for k=1:numInData
        if ismember(seqIdentifier(k,2),flanksToTrainOn)&&seqIdentifier(k,1)>0
            weights(k)=1;
        end
    end
    %weights(14)=1; %motif 1 with random flank, 14 18 22
    %weights(3)=1;%motif 2
    weights(motif1inds)=1;
    weights(motif2inds)=0;%
    weights(flankonlyinds)=1; %Don't fit on flank only;
    %weights(29)=0;
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
    
    
    paramIndexEnd=paramIndex(2:end)-1;
    paramIndexEnd(4)=numParams;
    
    lb=0*ones(1,numParams);
    ub=inf*ones(1,numParams);
    
    
    lb(1)=0;%0;%10^-5;%10^-2;
    
    
    
    lb(paramIndex(3):paramIndexEnd(3))=0;%0.01; %p_core_rel
    lb(paramIndex(4):paramIndexEnd(4))=0;%0.01; %p_flank_rel
    
    ub(1)=inf;%inf;%inf;%10^-5;%10^-2; %kon,max
    
    
    
    %ub(paramIndex(1):paramIndexEnd(1))=10^2; %KDmicro_cores
    %ub(paramIndex(2):paramIndexEnd(2))=10^2; %KDmicro_flank
    ub(paramIndex(3):paramIndexEnd(3))=inf; %p_core_rel
    ub(paramIndex(4):paramIndexEnd(4))=inf; %p_flank_rel
    paramsVec=[];
    residualVec=[];
    paramsAllVec=cell(1,numFits);
    residualAllVec=[];
    KDModelVec=[];
    koffModelVec=[];
    for j=1:numFits
        disp(['Fit num: ' num2str(j)]);
        [paramsCurr,residualCurr,paramsAllCurr,residualAllCurr]=fitCoreRepeatModelInfiniteTestingConNewNorm(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub,koffmicrofac,prelfac);
        paramsVec(j,:)=paramsCurr;
        residualVec(j,:)=residualCurr;
        paramsAllVec{j}=paramsAllCurr;
        residualAllVec(j,:)=residualAllCurr;
        params=paramsVec(j,:);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
        KDModelVec(j,:)=KDModel';
        koffModelVec(j,:)=koffModel';
    end
    params=paramsVec(1,:);
    
    [KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weightsAll,seqIdentifier,paramIndex);
    [~,~,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);
    
% % % % % %     figure(f1);
% % % % % %     hold off;
% % % % % %     plot(params(paramsFitted==1),'b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('params fit 1')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     
% % % % % %     
% % % % % %     figure(f2);
% % % % % %     hold off;
% % % % % %     plot((paramsVec(:,paramsFitted==1)./repmat(params(paramsFitted==1),[numFits 1]))','b.');
% % % % % %     test=gca;
% % % % % %     test.YScale='log';
% % % % % %     xlabel('Used parameter index');
% % % % % %     ylabel('(params fit 2)/(params fit 1)')
% % % % % %     title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2))]);
% % % % % %     %[KDModelTemp,koffModelTemp] = getManyCoreRepeatModelKDmicro(params,1,seqIdentifierTemp,paramIndex);
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(KDModel(weights~=0),KDData(weights~=0,:),'b.');
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     xlabel('Model K_D (nM)');
% % % % % %     ylabel('Experimental K_D (nM)');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     %maxL=max([yL xL]);
% % % % % %     %minL=min([yL xL]);
% % % % % %     %xlim([minL maxL]);
% % % % % %     %ylim([minL maxL]);
% % % % % %     
% % % % % %     figure(f4);
% % % % % %     hold off;
% % % % % %     title('Training data')
% % % % % %     plot(koffModel(weights~=0),koffData(weights~=0,:),'r.');
% % % % % %     xlabel('Model koff (s^-^1)');
% % % % % %     ylabel('Experimental koff (s^-^1)');
% % % % % %     
% % % % % %     
% % % % % %     figure(f3);
% % % % % %     hold on;
% % % % % %     plot(KDModel(motif1inds),KDData(motif1inds,:),'co');
% % % % % %     plot(KDModel(motif2inds),KDData(motif2inds,:),'mo');
% % % % % %     plot(KDModel(flankonlyinds),KDData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
% % % % % %     figure(f4);
% % % % % %     hold on;
% % % % % %     plot(koffModel(motif1inds),koffData(motif1inds,:),'co');
% % % % % %     plot(koffModel(motif2inds),koffData(motif2inds,:),'mo');
% % % % % %     plot(koffModel(flankonlyinds),koffData(flankonlyinds,:),'kx');
% % % % % %     xL=xlim();
% % % % % %     yL=ylim();
% % % % % %     maxL=max([yL xL]);
% % % % % %     minL=min([yL xL]);
% % % % % %     xlim([minL maxL]);
% % % % % %     ylim([minL maxL]);
% % % % % %     plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
    
    repsTest2{i}.paramsVec=paramsVec;
    repsTest2{i}.residualVec=residualVec;
    repsTest2{i}.paramsAllVec=paramsAllVec;
    repsTest2{i}.residualAllVec=residualAllVec;
    repsTest2{i}.KDModelVec=KDModelVec;
    repsTest2{i}.koffModelVec=koffModelVec;
    
    paramsRepsTest2(i,:)=params;
    KDModelRepsTest2(:,i)=KDModel;
    koffModelRepsTest2(:,i)=koffModel;
    
    figure(f5);
    nCompleted=i;
    hold off;
    plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
    hold on;
    test=gca;
    test.YScale='log';
    xlabel('paramIndex');
    ylabel('paramValue')
    title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    
end
nCompleted=numReps;
figure();
plot(paramsReps(1:nCompleted,paramsFitted==1)','b.');
hold on;
test=gca;
test.YScale='log';
xlabel('paramIndex');
ylabel('paramValue')
title(['residual 1: ' num2str(residualVec(1)) 'residual 2: ' num2str(residualVec(2)) ' bootNumber: ' num2str(i)]);
    
    


figure();
errorbar(mean(log10(paramsReps(1:nCompleted,paramsFitted==1))),std(log10(paramsReps(1:nCompleted,paramsFitted==1))),'bx');

%% trainging vs testing

perc=0.68;
koffData=nanmean(koffDataPho4Selected,2);
KDData=nanmean(KDDataPho4Selected,2);

ind1=find(motif1inds==1,1);
ind2=find(motif2inds==1,1);

%KD
modelCenters(1)=singleFitTest1.KDModel(ind1);
modelCenters(2)=singleFitTest2.KDModel(ind2);

modelVals1=zeros(1,numReps);
modelVals2=zeros(1,numReps);
for i=1:numReps
    modelVals1(i)=repsTest1{i}.KDModelVec(1,ind1);
    modelVals2(i)=repsTest2{i}.KDModelVec(1,ind2);
end 


sortedVals1=sort(modelVals1,'ascend');
sortedVals2=sort(modelVals2,'ascend');

startInd=round(((1-perc)/2)*numReps)+1;
endInd=numReps-startInd+1;

lowerBound1=sortedVals1(startInd);
upperBound1=sortedVals1(endInd);


lowerBound2=sortedVals2(startInd);
upperBound2=sortedVals2(endInd);

yvals=[mean(KDData(motif1inds)) modelCenters(1) mean(KDData(motif2inds)) modelCenters(2)];
%yvals=[mean(KDData(motif1inds)) mean(modelVals1) mean(KDData(motif2inds)) mean(modelVals2)];
yErrors=[std(KDData(motif1inds))/sqrt(numel(KDData(motif1inds))) std(modelVals1) std(KDData(motif2inds))/sqrt(numel(KDData(motif2inds))) std(modelVals2)];
yneg=yErrors;
ypos=yErrors;
yneg(2)=yvals(2)-lowerBound1;
ypos(2)=upperBound1-yvals(2);

yneg(4)=yvals(4)-lowerBound2;
ypos(4)=upperBound2-yvals(4);

figure();
hold on;
errorbar([1 2 4 5],yvals,yneg,ypos,'ro');
%plot(2*ones(1,numReps),modelVals1,'b.');
%plot(5*ones(1,numReps),modelVals2,'b.');
ylabel('KD (nM)')
xticks([1 2 4 5]);
test=gca;
test.XTickLabel{1}='motif1 data';
test.XTickLabel{2}='motif1 model';

test.XTickLabel{3}='motif2 data';
test.XTickLabel{4}='motif2 model';
test.YScale='log';
savefig('testing_KD_motif1_motif2')

%% koff
modelCenters(1)=singleFitTest1.koffModel(ind1);
modelCenters(2)=singleFitTest2.koffModel(ind2);

modelVals1=zeros(1,numReps);
modelVals2=zeros(1,numReps);
for i=1:numReps
    modelVals1(i)=repsTest1{i}.koffModelVec(1,ind1);
    modelVals2(i)=repsTest2{i}.koffModelVec(1,ind2);
end 


sortedVals1=sort(modelVals1,'ascend');
sortedVals2=sort(modelVals2,'ascend');

startInd=round(((1-perc)/2)*numReps)+1;
endInd=numReps-startInd+1;

lowerBound1=sortedVals1(startInd);
upperBound1=sortedVals1(endInd);


lowerBound2=sortedVals2(startInd);
upperBound2=sortedVals2(endInd);

yvals=[mean(koffData(motif1inds)) modelCenters(1) mean(koffData(motif2inds)) modelCenters(2)];
%yvals=[mean(koffData(motif1inds)) mean(modelVals1) mean(koffData(motif2inds)) mean(modelVals2)];
yErrors=[std(koffData(motif1inds))/sqrt(numel(koffData(motif1inds))) std(modelVals1) std(koffData(motif2inds))/sqrt(numel(koffData(motif2inds))) std(modelVals2)];
yneg=yErrors;
ypos=yErrors;
yneg(2)=yvals(2)-lowerBound1;
ypos(2)=upperBound1-yvals(2);

yneg(4)=yvals(4)-lowerBound2;
ypos(4)=upperBound2-yvals(4);

figure();
hold on;
errorbar([1 2 4 5],yvals,yneg,ypos,'ro');
%plot(2*ones(1,numReps),modelVals1,'b.');
%plot(5*ones(1,numReps),modelVals2,'b.');
ylabel('koff (s^-^1)')
xticks([1 2 4 5]);
test=gca;
test.XTickLabel{1}='motif1 data';
test.XTickLabel{2}='motif1 model';

test.XTickLabel{3}='motif2 data';
test.XTickLabel{4}='motif2 model';
test.YScale='log';

savefig('testing_koff_motif1_motif2')

%% params
paramsCenters=singleFit.param2(paramsFitted==1);
paramsBoots=zeros(numReps,sum(paramsFitted==1));
for i=1:numReps
    paramsBoots(i,:)=reps{i}.paramsVec(1,paramsFitted==1);
end

paramsCentersBoot=mean(paramsBoots,1);

numCurrParams=numel(paramsCenters);

startInd=round(((1-perc)/2)*numReps)+1;

endInd=numReps-startInd+1;

upperBounds=zeros(1,numCurrParams);
lowerBounds=zeros(1,numCurrParams);

for i=1:numCurrParams
   sortedVals=sort(paramsBoots(:,i),'ascend');
   upperBounds(i)=sortedVals(endInd);
   lowerBounds(i)=sortedVals(startInd);
end

figure();
errorbar(1:numCurrParams,paramsCenters,paramsCenters-lowerBounds,upperBounds-paramsCenters,'ro');
test=gca;
test.YScale='log';
xlabel('Param index');
ylabel('Param value');
title(['confidence interval : ' num2str(perc)]);
savefig('fitted_params_end')
%% write params to csv

cHeader = {'parameter_name' 'least_squares_fit' 'lowerbound_68percent_CI' 'upperbound_68percent_CI'}; %dummy header
%textHeader = strjoin(cHeader,','); %cHeader in text with commas
%write header to file

%fid = fopen('yourfile.csv','w'); 
%fprintf(fid,'%s\n',textHeader)
%fclose(fid)

%write data to end of file

numFittedParams=numel(paramsCenters);

paramsNames={'k_on_max' 'k_off_micro_motif' 'k_off_micro_mutated_motif','k_off_micro_flank1','k_off_micro_flank2','f_motif','f_mutated_motif','f_flank1','f_flank2'};

mydata=cell(numFittedParams,4);
for i=1:numFittedParams
    mydata{i,1}=paramsNames{i};
    mydata{i,2}=paramsCenters(i);
    mydata{i,3}=lowerBounds(i);
    mydata{i,4}=upperBounds(i);
end

T=cell2table(mydata,'VariableNames',cHeader);


writetable(T,'Max_fit.csv');


%% trainging eval
KDModel=singleFit.KDModel;
koffModel=singleFit.koffModel;

%KD
figure();

vals=KDModel;

hold on;

xlabel('Model K_D (nM)');
ylabel('Experimental K_D (nM)');
KDBoots=zeros(numReps,numInData);

for i=1:numReps
    KDBoots(i,:)=reps{i}.KDModelVec(1,:);
end 

lowerBounds=zeros(1,numInData);
upperBounds=zeros(1,numInData);

xPos=zeros(1,numInData);
xNeg=zeros(1,numInData);

startInd=round(((1-perc)/2)*numReps)+1;

endInd=numReps-startInd+1;

for j=1:numInData
    
    sortedVals=sort(KDBoots(:,j),'ascend');
    lowerBounds(j)=sortedVals(startInd);
    upperBounds(j)=sortedVals(endInd);
   
    
end
 xPos=upperBounds-vals';
 xNeg=vals'-lowerBounds;
%errorbar(KDModel(motif1inds),KDData(motif1inds,:),KDError(motif1inds),'co');
%errorbar(KDModel(motif2inds),KDData(motif2inds,:),KDError(motif2inds),'mo');
%errorbar(KDModel(flankonlyinds),KDData(flankonlyinds,:),KDError(flankonlyinds),'kx');
%errorbar(KDModel(motif1andflankinds),KDData(motif1andflankinds,:),KDError(motif1andflankinds),'cd');
%errorbar(KDModel(motif2andflankinds),KDData(motif2andflankinds,:),KDError(motif2andflankinds),'md');

errorbar(KDModel(motif1inds),KDData(motif1inds,:),KDError(motif1inds),KDError(motif1inds),xNeg(motif1inds),xPos(motif1inds),'co');
errorbar(KDModel(motif2inds),KDData(motif2inds,:),KDError(motif2inds),KDError(motif2inds),xNeg(motif2inds),xPos(motif2inds),'mo');
errorbar(KDModel(flankonlyinds),KDData(flankonlyinds,:),KDError(flankonlyinds),KDError(flankonlyinds),xNeg(flankonlyinds),xPos(flankonlyinds),'kx');
errorbar(KDModel(motif1andflankinds),KDData(motif1andflankinds,:),KDError(motif1andflankinds),KDError(motif1andflankinds),xNeg(motif1andflankinds),xPos(motif1andflankinds),'cd');
errorbar(KDModel(motif2andflankinds),KDData(motif2andflankinds,:),KDError(motif2andflankinds),KDError(motif2andflankinds),xNeg(motif2andflankinds),xPos(motif2andflankinds),'md');


test=gca;
test.XScale='log';
test.YScale='log';

xL=xlim();
yL=ylim();

maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);


plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);
savefig('KD_data_vs_model_end');

%% koff
figure();

vals=koffModel;

hold on;

xlabel('Model koff (s^-^1)');
ylabel('Experimental koff (s^-^1)');
koffBoots=zeros(numReps,numInData);

for i=1:numReps
    koffBoots(i,:)=reps{i}.koffModelVec(1,:);
end 

lowerBounds=zeros(1,numInData);
upperBounds=zeros(1,numInData);

xPos=zeros(1,numInData);
xNeg=zeros(1,numInData);

startInd=round(((1-perc)/2)*numReps)+1;

endInd=numReps-startInd+1;

for j=1:numInData
    
    sortedVals=sort(koffBoots(:,j),'ascend');
    lowerBounds(j)=sortedVals(startInd);
    upperBounds(j)=sortedVals(endInd);
   
    
end
 xPos=upperBounds-vals';
 xNeg=vals'-lowerBounds;
%errorbar(koffModel(motif1inds),koffData(motif1inds,:),koffError(motif1inds),'co');
%errorbar(koffModel(motif2inds),koffData(motif2inds,:),koffError(motif2inds),'mo');
%errorbar(koffModel(flankonlyinds),koffData(flankonlyinds,:),koffError(flankonlyinds),'kx');
%errorbar(koffModel(motif1andflankinds),koffData(motif1andflankinds,:),koffError(motif1andflankinds),'cd');
%errorbar(koffModel(motif2andflankinds),koffData(motif2andflankinds,:),koffError(motif2andflankinds),'md');

errorbar(koffModel(motif1inds),koffData(motif1inds,:),koffError(motif1inds),koffError(motif1inds),xNeg(motif1inds),xPos(motif1inds),'co');
errorbar(koffModel(motif2inds),koffData(motif2inds,:),koffError(motif2inds),koffError(motif2inds),xNeg(motif2inds),xPos(motif2inds),'mo');
errorbar(koffModel(flankonlyinds),koffData(flankonlyinds,:),koffError(flankonlyinds),koffError(flankonlyinds),xNeg(flankonlyinds),xPos(flankonlyinds),'kx');
errorbar(koffModel(motif1andflankinds),koffData(motif1andflankinds,:),koffError(motif1andflankinds),koffError(motif1andflankinds),xNeg(motif1andflankinds),xPos(motif1andflankinds),'cd');
errorbar(koffModel(motif2andflankinds),koffData(motif2andflankinds,:),koffError(motif2andflankinds),koffError(motif2andflankinds),xNeg(motif2andflankinds),xPos(motif2andflankinds),'md');


test=gca;
test.XScale='log';
test.YScale='log';

xL=xlim();
yL=ylim();

maxL=1.1*max([yL xL]);
minL=0.9*min([yL xL]);
xlim([minL maxL]);
ylim([minL maxL]);


plot([minL maxL],[minL maxL],'Color',[0.7 0.7 0.7]);

savefig('koff_data_vs_model_end')
%% Calculat passagetimes

solutionConc =1; %nM of DNA in solution
paramsUsed=singleFit.param2(paramsFitted==1);
paramsTest =paramsUsed([1 2 4 6 8]);
switchrate=paramsTest(3);
%paramsTest(3)=100;
p_flank_rel_vec=logspace(-4,3,100);
numpoints=numel(p_flank_rel_vec);
passagetimesCore=zeros(1,numpoints);
passagetimesCoreSwitch=zeros(1,numpoints);
for i=1:numel(p_flank_rel_vec)
    paramsNow=paramsTest;
    paramsNow(5)=p_flank_rel_vec(i);
    paramsNowSwitch=[paramsNow switchrate];
    [passagetimesCore(i)] = getPassagetimeCoreInfiniteTesting(paramsNow,solutionConc);
    [passagetimesCoreSwitch(i)] = getPassagetimeCoreInfiniteTestingSwitchAllow(paramsNowSwitch,solutionConc);
end

figure();
hold on;
plot(p_flank_rel_vec,passagetimesCore,'b');
plot(p_flank_rel_vec,passagetimesCoreSwitch,'r');
test=gca;
test.XScale='log';
ylabel('first passage time to core (s)');
xlabel('f\_flank');
savefig('firstpassage_example')

koffmicro_flank_vec=logspace(-4,3,100);

passagetimesCore2=zeros(1,numpoints);
passagetimesCoreSwitch2=zeros(1,numpoints);
for i=1:numel(p_flank_rel_vec)
    paramsNow=paramsTest;
    paramsNow(3)=koffmicro_flank_vec(i);
    paramsNowSwitch=[paramsNow koffmicro_flank_vec(i)];
    [passagetimesCore2(i)] = getPassagetimeCoreInfiniteTesting(paramsNow,solutionConc);
    [passagetimesCoreSwitch2(i)] = getPassagetimeCoreInfiniteTestingSwitchAllow(paramsNowSwitch,solutionConc);
end

figure();
hold on;
plot(p_flank_rel_vec,passagetimesCore2,'b');
plot(p_flank_rel_vec,passagetimesCoreSwitch2,'r');
test=gca;
test.XScale='log';
ylabel('first passage time to core (s)');
xlabel('koff\_mikro');


%% total site

solutionConc =1; %nM of DNA in solution
paramsUsed=singleFit.param2(paramsFitted==1);
paramsTest =paramsUsed([1 2 4 6 8]);
%paramsTest(1)=1;
switchrate=paramsTest(3);
%paramsTest(3)=100;
p_flank_rel_vec=logspace(-4,3,100);
numpoints=numel(p_flank_rel_vec);
passagetimesTotal=zeros(1,numpoints);

for i=1:numel(p_flank_rel_vec)
    paramsNow=paramsTest;
    paramsNow(5)=p_flank_rel_vec(i);
    paramsNowSwitch=[paramsNow switchrate];
    [passagetimesTotal(i)] = getPassagetimeTotalInfiniteTesting(paramsNow,solutionConc);
end

figure();
hold on;
plot(p_flank_rel_vec,passagetimesTotal,'b');
test=gca;
test.XScale='log';
ylabel('first passage time to motif+flanks (s)');
xlabel('f\_flank');
savefig('firstpassage_totalsite_example')


%% total site MFPT 2D

solutionConc =1; %nM of DNA in solution
paramsUsed=singleFit.param2(paramsFitted==1);
paramsTest =paramsUsed([1 2 4 6 8]);
%paramsTest(1)=1;
switchrate=paramsTest(3);
%paramsTest(3)=100;
p_flank_rel_vec=logspace(-2,2,100);
p_core_rel_vec=logspace(-2,2,100);%logspace(log10(0.7),log10(0.9),50);
numpointsF=numel(p_flank_rel_vec);
numpointsC=numel(p_core_rel_vec);
passagetimesTotal2D=zeros(numpointsC,numpointsF);

for i=1:numpointsF
    for j=1:numpointsC
    paramsNow=paramsTest;
    paramsNow(5)=p_flank_rel_vec(i);
    paramsNow(4)=p_core_rel_vec(j);
    paramsNowSwitch=[paramsNow switchrate];
    [passagetimesTotal2D(j,i)] = getPassagetimeTotalInfiniteTesting(paramsNow,solutionConc);
    end
end

figure();
hold on;
surf(p_flank_rel_vec,p_core_rel_vec,passagetimesTotal2D/(1/(paramsTest(1)*solutionConc)));
test=gca;
test.XScale='log';
test.YScale='log';
test.ColorScale='log';
cc=colorbar;
colormap hot
cc.Ticks=[1:10 20 30 40 50 60 70 80];
test.CLim=[1 51];
ylabel(cc,'first passage time to motif+flanks 1/(konmax*DNA)');
xlabel('f\_flank');
ylabel('f\_motif');
shading interp;
savefig('firstpassage_totalsite_norm2D_example')
%% total site normalied time


figure();
hold on;
plot(p_flank_rel_vec,passagetimesTotal/(1/(paramsTest(1)*solutionConc)),'b');
test=gca;
test.XScale='log';
ylabel('first passage time to motif+flanks 1/(konmax*DNA)');
xlabel('f\_flank');
savefig('firstpassage_totalsite_norm_example')

mydata=[p_flank_rel_vec' passagetimesTotal' passagetimesTotal'/(1/(paramsTest(1)*solutionConc))];

%cHeader={'f_flank','nonnormalized_time_(s)' 'normalized_time_(1/(k_on,max*[DNA]))'};
cHeader={'f_flank','nonnormalized_time' 'normalized_time'};
T=table(p_flank_rel_vec',passagetimesTotal',passagetimesTotal'/(1/(paramsTest(1)*solutionConc)),'VariableNames',cHeader);


writetable(T,'passage_time_curves.csv');
%% additative f_flank and f_motif to add sites together (what is assuemd in the mode fitting)
times=[0:1:800];
%times=[0:0.00005:1];
%paramsTest =[konmax koffmicro_core koffmicro_flank f_core _rel f_flank_rel]
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
paramsUsed=singleFit.param2(paramsFitted==1);
paramsTest =paramsUsed([1 2 4 6 8]);
%paramsTest =paramsUsed([1 2 2 6 6]);
solutionConc=100;
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient

%convert paramters 

f_flank_rel_vec=logspace(-4,3,100);
numpoints=numel(f_flank_rel_vec);
occupancyTotal=zeros(3,numpoints); %Testing state, core, flank
occupancyTotalMotifOnly=zeros(3,numpoints);
occupancyTotalFlankOnly=zeros(3,numpoints);
konTotal=zeros(1,numpoints);
konTotalMotifOnly=zeros(1,numpoints);
konTotalFlankOnly=zeros(1,numpoints);
for i=1:numpoints
    
%both motif and flank
%f_flank_curr=paramsTest(5);
f_flank_curr=f_flank_rel_vec(i);
p_core = paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
p_flank= f_flank_curr/(1+paramsTest(4)+f_flank_curr);
params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];
bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotal(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_1=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_1=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_1=1/offMeanFirstPassageInitial;

konTotal(i)=konMeanFirstPassage_1;

%motif only
f_flank_curr=0;
p_core = paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
p_flank= f_flank_curr/(1+paramsTest(4)+f_flank_curr);
params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotalMotifOnly(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_2=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_2=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_2=1/offMeanFirstPassageInitial;
konTotalMotifOnly(i)=konMeanFirstPassage_2;



%flamnk only only
%f_flank_curr=paramsTest(5);
f_flank_curr=f_flank_rel_vec(i);
p_core = 0;%paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
p_flank= f_flank_curr/(1+0+f_flank_curr);
params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotalFlankOnly(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_3=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_3=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_3=1/offMeanFirstPassageInitial;

konTotalFlankOnly(i)=konMeanFirstPassage_3;

end

figure();
hold on;
plot(f_flank_rel_vec,occupancyTotal(2,:),'r');
plot(f_flank_rel_vec,occupancyTotal(3,:),'c');
plot(f_flank_rel_vec,(occupancyTotal(2,:)+occupancyTotal(3,:)),'b');
test=gca;
test.XScale='log';


ylabel('occupancy');
xlabel('f\_flank');

legend('on motif','on flank','on motif + flank');

figure();
hold on;
plot(f_flank_rel_vec,occupancyTotalMotifOnly(2,:),'k--');
plot(f_flank_rel_vec,occupancyTotalFlankOnly(3,:),'k-.');
plot(f_flank_rel_vec,occupancyTotalMotifOnly(2,:)+occupancyTotalFlankOnly(3,:),'k');
plot(f_flank_rel_vec,(occupancyTotal(2,:)+occupancyTotal(3,:)),'b');
test=gca;
test.XScale='log';


ylabel('occupancy');
xlabel('f\_flank');

legend('motif only DNA','flanks only DNA','on motif +flanks, independent sites','on motif + flanks, DNA with both sites');

%funcValsOn=solveODEs(Aon,bon,x0on,times);
%funcValsOff=solveODEs(Aoff,boff,x0off,times);

%% (1-p) way of adding sites together
times=[0:1:800];
%times=[0:0.00005:1];
%paramsTest =[konmax koffmicro_core koffmicro_flank f_core _rel f_flank_rel]
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
paramsUsed=singleFit.param2(paramsFitted==1);
paramsTest =paramsUsed([1 2 4 6 8]);
%paramsTest =paramsUsed([1 2 2 6 6]);
solutionConc=100;
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient

%convert paramters 

f_flank_rel_vec=logspace(-4,3,100);
numpoints=numel(p_flank_rel_vec);
occupancyTotal=zeros(3,numpoints); %Testing state, core, flank
occupancyTotalMotifOnly=zeros(3,numpoints);
occupancyTotalFlankOnly=zeros(3,numpoints);
konTotal=zeros(1,numpoints);
konTotalMotifOnly=zeros(1,numpoints);
konTotalFlankOnly=zeros(1,numpoints);
for i=1:numpoints
    
%both motif and flank
%f_flank_curr=paramsTest(5);
f_flank_curr=f_flank_rel_vec(i);
%p_core = paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
%p_flank= f_flank_curr/(1+paramsTest(4)+f_flank_curr);

p_core_single = paramsTest(4)/(1+paramsTest(4));
p_flank_single = f_flank_curr/(1+f_flank_curr);

ptot=1-(1-p_core_single)*(1-p_flank_single);

%Approximatom, probably not true, need to run simulations to check what
%this actually is
p_core=p_flank_single*((1-p_core_single))+p_flank_single*p_core_single/2;
p_flank=p_core_single*((1-p_flank_single))+p_flank_single*p_core_single/2;

params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];
bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotal(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_1=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_1=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_1=1/offMeanFirstPassageInitial;

konTotal(i)=konMeanFirstPassage_1;

%motif only
f_flank_curr=0;
p_core = paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
p_flank= f_flank_curr/(1+paramsTest(4)+f_flank_curr);
params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotalMotifOnly(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_2=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_2=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_2=1/offMeanFirstPassageInitial;
konTotalMotifOnly(i)=konMeanFirstPassage_2;



%flamnk only only
%f_flank_curr=paramsTest(5);
f_flank_curr=f_flank_rel_vec(i);
p_core = 0;%paramsTest(4)/(1+paramsTest(4)+f_flank_curr);
p_flank= f_flank_curr/(1+0+f_flank_curr);
params =[paramsTest(1) 10^6 paramsTest(2) paramsTest(3) p_core p_flank];

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];

x0off=Aon\(-bon);

occupancyTotalFlankOnly(:,i)=x0off;

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage_3=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;

ssProb3=1/(1+(params(6)*params(3)/(params(4)*params(5))));
ssProb4=1-ssProb3;
offMeanFirstPassage2=ssProb3*offMeanFirstPassageCore+ssProb4*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage_3=1/offMeanFirstPassage;

koffMeanFirstPassageInitial_3=1/offMeanFirstPassageInitial;

konTotalFlankOnly(i)=konMeanFirstPassage_3;

end

figure();
hold on;
plot(f_flank_rel_vec,occupancyTotal(2,:),'r');
plot(f_flank_rel_vec,occupancyTotal(3,:),'c');
plot(f_flank_rel_vec,(occupancyTotal(2,:)+occupancyTotal(3,:)),'b');
test=gca;
test.XScale='log';


ylabel('occupancy');
xlabel('f\_flank');

legend('on motif','on flank','on motif + flank');

figure();
hold on;
plot(f_flank_rel_vec,occupancyTotalMotifOnly(2,:),'k--');
plot(f_flank_rel_vec,occupancyTotalFlankOnly(3,:),'k-.');
plot(f_flank_rel_vec,occupancyTotalMotifOnly(2,:)+occupancyTotalFlankOnly(3,:),'k');
plot(f_flank_rel_vec,(occupancyTotal(2,:)+occupancyTotal(3,:)),'b');
test=gca;
test.XScale='log';


ylabel('occupancy');
xlabel('f\_flank');

legend('motif only DNA','flanks only DNA','on motif +flanks, independent sites','on motif + flanks, DNA with both sites');

%funcValsOn=solveODEs(Aon,bon,x0on,times);
%funcValsOff=solveODEs(Aoff,boff,x0off,times);

