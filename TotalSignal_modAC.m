function s=TotalSignal_modAC(lambda,z,Conformers,T,dW)
% s=TotalSignal_modAC(lambda,z,Conformers,T,dW)
% dW: width of the suppression function (MHz)
dW = 2*pi*(dW/1000); %MHz to G(rad/s)
r=Conformers.Distance;
D=0.326983340102601;                  %Dipolar Interaction Constant, G(rad/s)·nm³=rad/ns·nm³
Z=repmat(z',[size(lambda,1) 1 size(lambda,3)]); %[nOffsets nTheta nConformers]
R=repmat(reshape(r,[1 1 size(r)]),[size(lambda,1) size(lambda,2) 1]); %[nOffsets nTheta nConformers]
s=zeros(size(lambda,1),size(T,1));
wdd = D.*(3.*Z.^2-1)./R.^3; %dipolar coupling frequency
G = exp(-((wdd.^2)./(dW.^2))); %suppression window
if isfield(Conformers,'Probability')
    P = Conformers.Probability/sum(Conformers.Probability);
    P = repmat(reshape(P,[1 1 size(P)]),[size(lambda,1) size(lambda,2) 1]);
    for k=1:size(T,1)
        p=sum(P.*lambda.*G.*(cos(wdd.*T(k))-1),3);
        s(:,k)=1+trapz(z,p,2);
    end
else
    for k=1:size(T,1)
        p=mean(lambda.*G.*(cos(wdd.*T(k))-1),3);
        s(:,k)=1+trapz(z,p,2);
    end
end