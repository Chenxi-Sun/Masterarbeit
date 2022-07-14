bp1=1;
NConf=100;
h=33;
r=5.95;
ss=normgauss(NConf,0.01);
Dynamik.radius=ones(NConf,1).*r+ss(:,9);
ss=normgauss(NConf,0.000001);
Dynamik.hpitch=ones(NConf,1).*h+ss(:,7);
% if Bewegung.bend==0
Dynamik.alpha=0;
% else
%  ss=normgauss(NConf,Bewegung.bend);
%  Dynamik.alpha=zeros(NConf,1)+ss(:,8);
% end

p = haltonset(1, 'Leap', 10,'Skip',20);  %%%generate points in space (Halton sequence)
  Dynamik.Drehung = net(p, NConf)*2*pi;
  
  
 %Berechnet die Konformere f?r die DNA(1,9)
bp1=1;
bp2=9;
%SLframe berechnet das Achsensystem der Spinlabel
% [SLA,X1,Y1,Z1]=SLframealle(1,Dynamik,bp1,Bewegung.Bew9,IMUSpinlabel,1,2,3,'DNA');

StrangA.dh=1.18; StrangB.dh=2.19;
StrangA.phi=-0.554; StrangB.phi=-4.477;
StrangA.c=4.85;StrangB.c=4.85;

A=StrangA.c./sqrt(Dynamik.hpitch.^2+(2.*pi.*Dynamik.radius).^2);
fi=2.*pi.*A*bp1+StrangA.phi; 

X=repmat(Dynamik.radius,size(bp1)).*cos(fi);       
Y=repmat(Dynamik.radius,size(bp1)).*sin(fi);       
Z=Dynamik.hpitch.*A*bp1 + StrangA.dh;    

s(:,1,:)=X;  s(:,2,:)=Y;  s(:,3,:)=Z;



[SLB,X2,Y2,Z2]=SLframealle(-1,Dynamik,bp2,Bewegung.Bew9,IMUSpinlabel,2,3,1,'DNA');
R=SLA-SLB;

%aus dem Achsensystem werden die Eulerwinkel bestimmt
[o1,o2,distance]=RuCoordv1(R,X1,Y1,Z1,X2,Y2,Z2);
%und in eine Struktur f?r MainPELDOR umgewandelt
Conformers.EulerAngles.R1=o1;
Conformers.EulerAngles.R2=o2;
Conformers.Distance=distance./10;
Conformers.Distance=Conformers.Distance*1.02;
