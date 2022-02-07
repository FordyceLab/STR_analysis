function [squaredDiffTot] = diffSquaredCoreRepeatModelInfiniteTesting(params,koffData,KDData,weights,seqIdentifier,paramIndex)
%params = [%konmax %all_koffmicro_core all_koffmicro_flank all_p_core_rel all_p_flank_rel]
[numInData,numReps]=size(koffData);

[KDModel,koffModel] = getManyCoreRepeatModelInfiniteTesting(params,weights,seqIdentifier,paramIndex);
%if sum(isnan(KDModel))>0||sum(isnan(KoffModel))>0
%    squaredDiffTot=inf;
%else
squaredDiffsKD=zeros(numInData,numReps);
squaredDiffskoff=zeros(numInData,numReps);
for i=1:numInData
    for j=1:numReps
        if weights(i)>0
            if isfinite(KDModel(i))&&isfinite(koffModel(i))
            squaredDiffsKD(i,j)=(((KDModel(i)-KDData(i,j))/KDData(i,j))^2)*weights(i);
            squaredDiffskoff(i,j)=(((koffModel(i)-koffData(i,j))/koffData(i,j))^2)*weights(i);
            else
              squaredDiffsKD(i,j)=10^8;  
              squaredDiffskoff(i,j)=10^8;
            end    
          
        end
    end
end
    squaredDiffTot=nansum(nansum(squaredDiffsKD))+nansum(nansum(squaredDiffskoff));
end




