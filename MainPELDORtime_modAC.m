function Result=MainPELDORtime_modAC(EP,Conformers,t,dW)
while max(t)<500  
    t=t*10;             %set t in ns
end
    
%The experimental parameter (EP) can be either directly insert with the needed
%structur (see documentation) or is given in a script with the related name
%as a char.
if ischar(EP)==1
[EP]=ExpParameters(EP);
end
try
    EP=SetEP(EP);
catch
    disp('hallo')
end
% The magnetic parameters (MP) like g-tensor and hyperfine tensor are given in
% this function.
MP=MagneticParameters('Two Nitroxides');
% MP.R1.Sigma =5e-04;
% MP.R2.Sigma =5e-04;

%The excitationprofiles for the pulses
profiles=Excitation(EP,MP,'Derived from the pulse shape');

%Spherical Coordinates for the powder averaging
z=(0:0.01:1)'; theta=acos(z); phi=(0:0.02*pi:2*pi)';  
m1=(-1:1:1)'; m2=(-1:1:1)';%quantum numbers m_s

%Grid Variables
[V.THETA,V.PHI,V.M1,V.M2]=ndgrid(theta,phi,m1,m2);
%Refocused Echo Magnitude
v0=GridRefocusedEcho(profiles,EP.Settings,V,MP);
V0=repmat(trapz(z,trapz(phi,sum(sum(v0,5),4),3),2),[1 size(z,1)]);
%calculation of lambda (orientation function, pump pulse efficiency function) for each conformer
%Conformer loop
Nconf=size(Conformers.Distance,1);
L=zeros(size(EP.Settings.B0,1),numel(z),Nconf);

for j=1:Nconf
    ea1=Conformers.EulerAngles.R1(j,:); 
    ea2=Conformers.EulerAngles.R2(j,:);
    lambda=GridPELDOR(profiles,EP.Settings,V,ea1,ea2,MP);
    L(:,:,j)=trapz(phi,sum(sum(lambda,5),4),3)./V0; %integrate over the phi angle
%     if round(j/200)==j/200
%     plot(z',L(:,:,j)); pause(0.01);
%     end
    %disp(j); drawnow;
end

Result.LAMBDAS=L;
x=size(t);
%calclation of the resulting PELDOR time trace,
if x(1)>x(2)
 s=TotalSignal_modAC(L,z,Conformers,t,dW)';
else
 s=TotalSignal_modAC(L,z,Conformers,t',dW)';   
end

%Output variables
Result.Conformers=Conformers; Result.MagnPar=MP; Result.ExpPar=EP;
Result.time=t; Result.Sexp=s;

