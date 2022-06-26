% function Result = Spinlabelgeometry(MD1,MD2,EP,zeit)
%%Spinlabel.1
A = importdata(MD1);

[m,n]=size(A);
N  = [A(:,1) A(:,2) A(:,3)];
C2 = [A(:,4) A(:,5) A(:,6)];
C4 = [A(:,7) A(:,8) A(:,9)];
C5 = [A(:,10) A(:,11) A(:,12)];

NC5 = C5-N;
NC2 = C2-N;

%Bereche die Achse
xaxis1 = zeros(m,3);        
yaxis1 = zeros(m,3);
zaxis1 = zeros(m,3);

for ra=1:m
    xaxis1(ra,:) = NC5(ra,:)/norm(NC5(ra,:)); 
    NC2(ra,:) = NC2(ra,:)/norm(NC2(ra,:));
    zaxis1(ra,:) = cross(xaxis1(ra,:),NC2(ra,:));
    zaxis1(ra,:) = zaxis1(ra,:)/norm(zaxis1(ra,:));
    yaxis1(ra,:) = cross(zaxis1(ra,:),xaxis1(ra,:));  
%Rechne den Mittelpunkt der 5-Ring als Nullpunkt 
    M1(ra,:) = C5(ra,:)+0.5*(C4(ra,:)-N(ra,:)); %Mittelpunkt zwischen C4 und C5
    Mittelpunkt(ra,:) = M1(ra,:)+0.688*norm(C4(ra,:)-C5(ra,:))*xaxis1(ra,:); %Mittelpunkt der 5-Ring
    %Drehe die Achse
    theta = 27/360*2*pi; 
    [xaxis1_dre(ra,:) yaxis1_dre(ra,:) zaxis1_dre(ra,:)] = drehnung(xaxis1(ra,:),yaxis1(ra,:),zaxis1(ra,:),theta);
%Rechne M in der NO-Verbindung
    Ns(ra,:) = Mittelpunkt(ra,:)+xaxis1_dre(ra,:)*(3.959*norm(C4(ra,:)-C5(ra,:))); %Position des N in der NO-Verbindung
    Os(ra,:) = Mittelpunkt(ra,:)+xaxis1_dre(ra,:)*(3.959*norm(C4(ra,:)-C5(ra,:))+1.4); %Position des O in der NO-Verbindung
    M(ra,:) = 0.5*(Ns(ra,:)+Os(ra,:)); %Position des Mittelpunkt der NO-Verbindung
end


%%%%%%%%%%%%%%%%%%%%%%
%%Spinlabel.2
A2 = importdata(MD2);

[m2,n2]=size(A2);
N_2  = [A2(:,1) A2(:,2) A2(:,3)];
C2_2 = [A2(:,4) A2(:,5) A2(:,6)];
C4_2 = [A2(:,7) A2(:,8) A2(:,9)];
C5_2 = [A2(:,10) A2(:,11) A2(:,12)];

NC5_2   = C5_2-N_2;
NC2_2   = C2_2-N_2;

%Bereche die Achse
xaxis1_2 = zeros(m,3);
yaxis1_2 = zeros(m,3);
zaxis1_2 = zeros(m,3);

for ra=1:m
    xaxis1_2(ra,:) = NC5_2(ra,:)/norm(NC5_2(ra,:)); 
    NC2(ra,:) = NC2_2(ra,:)/norm(NC2_2(ra,:));
    zaxis1_2(ra,:) = cross(xaxis1_2(ra,:),NC2_2(ra,:));
    zaxis1_2(ra,:) = zaxis1_2(ra,:)/norm(zaxis1_2(ra,:));
    yaxis1_2(ra,:) = cross(zaxis1_2(ra,:),xaxis1_2(ra,:));    
%Rechne den Mittelpunkt der 5-Ring als Nullpunkt 
    M1_2(ra,:) = C5_2(ra,:)+0.5*(C4_2(ra,:)-N_2(ra,:)); %Mittelpunkt zwischen C4 und C5
    Mittelpunkt_2(ra,:) = M1_2(ra,:)+0.688*norm(C4_2(ra,:)-C5_2(ra,:))*xaxis1_2(ra,:); %Mittelpunkt der 5-Ring
    %Drehe die Achse
    
    theta_2 = 27/360*2*pi; 
    [xaxis1_2_dre(ra,:) yaxis1_2_dre(ra,:) zaxis1_2_dre(ra,:)] = drehnung(xaxis1_2(ra,:),yaxis1_2(ra,:),zaxis1_2(ra,:),theta_2);
%Rechne M in der NO-Verbindung
    Ns_2(ra,:) = Mittelpunkt_2(ra,:)+xaxis1_2_dre(ra,:)*(3.959*norm(C4_2(ra,:)-C5_2(ra,:))); %Position des N in der NO-Verbindung
    Os_2(ra,:) = Mittelpunkt_2(ra,:)+xaxis1_2_dre(ra,:)*(3.959*norm(C4_2(ra,:)-C5_2(ra,:))+1.4); %Position des O in der NO-Verbindung
    M2(ra,:) = 0.5*(Ns_2(ra,:)+Os_2(ra,:)); %Position des Mittelpunkt der NO-Verbindung
end

%R-Vektor



R = M-M2;
[o1,o2,r]=RuCoordv1(R,xaxis1_dre,yaxis1_dre,zaxis1_dre,xaxis1_2_dre,yaxis1_2_dre,zaxis1_2_dre);

Conformers.M=M;
Conformers.M2=M2;
% Conformers.R=R;
Conformers.EulerAngles.R1=[o1(:,:)];
% Conformers.EulerAngles.R1(:,2)=Conformers.EulerAngles.R1(:,2)+ones(length(Conformers.EulerAngles.R1(:,2)),1)*deg2rad(10);
Conformers.EulerAngles.R2=[o2(:,:)];
% Conformers.EulerAngles.R2(:,2)=Conformers.EulerAngles.R2(:,2)+ones(length(Conformers.EulerAngles.R2(:,2)),1)*deg2rad(10);
Conformers.Distance = r/10;
zeiten = zeit*1000;;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the program
% Result.dista = [vecnorm(N-N_2,2,2) vecnorm(C2-C2_2,2,2) vecnorm(C4-C4_2,2,2) vecnorm(C5-C5_2,2,2)];
% end 
