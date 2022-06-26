% function Result = Spinlabelgeometry(MD1,MD2,EP,zeit)
%%Spinlabel.1
MD1=load('Z:\Students\ChSun\MD_simulations_DNA\BP1_DNA_BSC1_odd.txt')
A = MD1;

[m,n]=size(A);
N  = [A(:,1) A(:,2) A(:,3)];
C2 = [A(:,4) A(:,5) A(:,6)];
C4 = [A(:,7) A(:,8) A(:,9)];
C5 = [A(:,10) A(:,11) A(:,12)];

%C1a points,C1b 
% NC2=C2-N;
NC2=(C2(1,:)-N(1,:))/1.4;
NC5=(C5(1,:)-N(1,:))/2.4;
z=cross(NC5,NC2);

NC1=NC2*cos(242.1/360*2*pi)+NC5*sin(242.1/360*2*pi);
NC1_rot=NC2*cos(278/360*2*pi)+NC5*sin(278/360*2*pi);
NC1_perp=(dot(NC1,NC1_rot)/norm(NC1_rot)^2)*NC1_rot;

C12_unit=NC1-NC1_perp;
C12_unit=-C12_unit;
C12_unit=C12_unit/norm(C12_unit);
C1=N(1,:)+NC1*1.5;
C2=N(1,:)+NC1*1.5+C12_unit*10.7;
C12=C2-C1;

%global z_axis



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
end


%%%%%%%%%%%%%%%%%%%%%%
%%Spinlabel.2
MD2=load('Z:\Students\ChSun\MD_simulations_DNA\BP9_DNA_BSC1_odd.txt')
A2 = MD2;

[m2,n2]=size(A2);
N_2  = [A2(:,1) A2(:,2) A2(:,3)];
C2_2 = [A2(:,4) A2(:,5) A2(:,6)];
C4_2 = [A2(:,7) A2(:,8) A2(:,9)];
C5_2 = [A2(:,10) A2(:,11) A2(:,12)];
NC5_2   = C5_2-N_2;
NC2_2   = C2_2-N_2;
%C1 points 

C1b=(N_2-C4_2)/norm(N_2-C4_2)*1.5+N_2;

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
end

%R-Vektor
alpha=75.5/360*2*pi;               %first rotation angle alpha
beta=6.8/360*2*pi;               %second rotation angle beta

for ra=1:m
    X_axis1_Spin1(ra,:)=C1b(ra,:)-C1a(ra,:);
    Y_axis_Spin1(ra,:)=cross(X_axis1_Spin1(ra,:),zaxis1(ra,:))
    Y_axis_Spin1(ra,:)=Y_axis_Spin1(ra,:)/norm(Y_axis_Spin1(ra,:));
    Z_axis_Spin1(ra,:)=cross(Y_axis_Spin1(ra,:),X_axis1_Spin1(ra,:));
    Z_axis_Spin1(ra,:)=Z_axis_Spin1(ra,:)/norm(Z_axis_Spin1(ra,:));
    %%first rotation
    X_axis2_Spin1(ra,:)=X_axis1_Spin1(ra,:).*cos(alpha)+Y_axis_Spin1(ra,:).*sin(alpha);
    %%second rotation
    x_axis3_Spin1(ra,:)=X_axis2_Spin1(ra,:).*cos(beta)+Z_axis_Spin1(ra,:).*sin(beta);
    %%find Spin1
    Spin_1(ra,:)=C1a(ra,:)+x_axis3_Spin1(ra,:).*11.4;  %length of C1 and spinlabel is 11.4A      
end

for ra=1:m
    X_axis1_Spin2(ra,:)=C1a(ra,:)-C1b(ra,:);
    Y_axis_Spin2(ra,:)=cross(X_axis1_Spin2(ra,:),zaxis1_2(ra,:))
    Y_axis_Spin2(ra,:)=Y_axis_Spin2(ra,:)/norm(Y_axis_Spin2(ra,:));
    Z_axis_Spin2(ra,:)=cross(Y_axis_Spin2(ra,:),X_axis1_Spin2(ra,:));
    Z_axis_Spin2(ra,:)=Z_axis_Spin2(ra,:)/norm(Z_axis_Spin2(ra,:));
    %%first rotation
    X_axis2_Spin2(ra,:)=X_axis1_Spin2(ra,:).*cos(alpha)-Y_axis_Spin2(ra,:).*sin(alpha);
    %%second rotation
    x_axis3_Spin2(ra,:)=X_axis2_Spin2(ra,:).*cos(beta)-Z_axis_Spin2(ra,:).*sin(beta);
    %%find Spin1
    Spin_2(ra,:)=C1b(ra,:)+x_axis3_Spin2(ra,:).*11.4;  %length of C1 and spinlabel is 11.4A      
end

R=Spin_1-Spin_2;

% R = M-M2;
[o1,o2,r]=RuCoordv1(R,xaxis1,yaxis1,zaxis1,xaxis1_2,yaxis1_2,zaxis1_2);

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