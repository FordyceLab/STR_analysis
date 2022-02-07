function [KDModel,koffModel,paramsFitted] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex)

%params =  [%konmax %all_koffmicro_core all_koffmicro_flank all_p_core_rel all_p_flank_rel]

[numInData]=size(seqIdentifier,1);
paramsFitted=zeros(1,numel(params));
KDModel=NaN*zeros(numInData,1);
koffModel=NaN*zeros(numInData,1);
paramsCurr=zeros(1,5);
for i=1:numInData
    if weights(i)>0
       
        paramsCurr(1)=params(1); %konmax
        paramsFitted(1)=1;
        coreIndex = seqIdentifier(i,1);
        if coreIndex==0
            paramsCurr(2)=1;%Shouldn matter what you set here
            paramsCurr(4)=0;
        else
            paramsCurr(2)=params(paramIndex(1)-1+coreIndex);% koffmicro_core
            paramsCurr(4)=params(paramIndex(3)-1+coreIndex);%p_core_rel
            paramsFitted(paramIndex(1)-1+coreIndex)=1;
            paramsFitted(paramIndex(3)-1+coreIndex)=1;
        end
        
        flankIndex = seqIdentifier(i,2);
        if flankIndex==0
            paramsCurr(3)=1;%Shouldn matter what you set here
            paramsCurr(5)=0;
        else
            paramsCurr(3)=params(paramIndex(2)-1+flankIndex);% koffmicro_core
            paramsCurr(5)=params(paramIndex(4)-1+flankIndex);%p_core_rel
            paramsFitted(paramIndex(2)-1+flankIndex)=1;
            paramsFitted(paramIndex(4)-1+flankIndex)=1;
        end
        if paramsCurr(5)>0||paramsCurr(4)>0
            [KDModel(i),koffModel(i)] = getRatesCoreRepeatFastInfiniteTesting(paramsCurr);
        else
            KDModel(i)=NaN;
            koffModel(i)=NaN;
        end
    end
end
end

