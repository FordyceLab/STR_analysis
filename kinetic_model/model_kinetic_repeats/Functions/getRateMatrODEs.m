function [RateMatr] = getRateMatrODEs(pvals,totRateConst)
%Returns the rate matrix fo an dissocation experiment where the protein concentration in solution doesnt change 
dim=numel(totRateConst);
%dim = number of states in model -1

RateMatr=zeros(dim,dim);

pvalsnew=[0 pvals];

RateMatr(1,1)=-totRateConst(1);
RateMatr(1,2)=pvalsnew(2)*totRateConst(2);

for i=2:dim-1
    RateMatr(i,i-1)=(1-pvalsnew(i-1))*totRateConst(i-1); 
    RateMatr(i,i)=-totRateConst(i);
    RateMatr(i,i+1)=pvalsnew(i+1)*totRateConst(i+1);

end
%onRateMatr(dim-1,1:dim-1)=-totRateConst(dim);

%onRateMatr(dim-1,dim-2)=onRateMatr(dim-1,dim-2)+(1-pvalsnew(dim-2))*totRateConst(dim-2);

%onRateMatr(dim-1,dim-1)=onRateMatr(dim-1,dim-1)-totRateConst(dim-1);

%onRateMatr(dim-1,dim)=totRateConst(dim);

RateMatr(dim,dim-1)=(1-pvalsnew(dim-1))*totRateConst(dim-1);
RateMatr(dim,dim)=-totRateConst(dim);



end

