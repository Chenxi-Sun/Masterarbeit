
% function [Result,R_mean] = AMPELDOR_CdotDNA(nd,EP,zeit)
% function [Result,R_mean,sigma] = AMPELDOR_CdotDNA(nd,EP,zeit)
% sigma_y=14;
function [Result] = AMPELDOR_CdotDNA(sigma_y,nd,EP,zeit)
% function [Result] = AMPELDOR_CdotDNA(nd,sigma_r,EP,zeit)
% function [Result] = AMPELDOR_CdotDNA(nd,sigma_h,EP,zeit)
% function [Result] = AMPELDOR_CdotDNA(alpha,beta,EP,zeit)
% for i=5:14
nd=nd-4;
str='B';
% nd=i-4;
sigma_z=5;
% sigma_y=0;

%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

%%initial position for C from pymol

% %parameter for expression for C1-points at DNA from pymol
r0=5.68;
h0=33.75; %h0 for C
b=0.6288;
c1_x=-1.691;
c1_y=-0.1296;
c2_x=0.7453;
c2_y=2.324;
d=3.375;
e1=-69.2;
e2=-69.17;

%helix.A
x1=r0*sin(b.*n_bp+c1_x);
y1=r0*sin(b.*n_bp+c1_y);
z1=d.*n_bp+e1;  

%helix.B
x2=r0*sin(b.*n_bp+c2_x);
y2=r0*sin(b.*n_bp+c2_y);
z2=d.*n_bp+e2;

%end-end length
L0=z2(20)-z1(1);

%contour length C
C=sqrt((2*pi*r0)^2*L0^2/h0^2+L0^2); %countor length (AM paper)

%n_turns
% n_turns=b.*19/2/pi;  %the range of b.*n_bp+c1 = how much 2pi = how much turns
n_turns=(z1(20)-z1(1))/h0; %the total height of a helix=n_turns*pitch height

%z_axis
% Z=[x1(11)-x1(1) y1(11)-y1(1) z1(11)-z1(1)];
% z_axis=Z/norm(Z);
z_axis=[0 0 1];

%%rotationsangle 1st SL
alpha_1=76.20/360*2*pi;
beta_1=6.65/360*2*pi;
%%rotationsangle 2nd SL 5-14 bp;
alpha_2=[76.1446 76.1787 76.1766 76.2030 76.1871 76.1902 76.1983 76.1741 76.1634 76.1681]./360.*2.*pi;
beta_2=[6.6477 6.6232 6.6401 6.6957 6.7435 6.7840 6.7882 6.7712 6.7300 6.6836]./360.*2.*pi;

for k=1:500
switch (str)
    case 'B'
% parameter
% deltar=normrnd(0,sigma_r);   %standard deviation of radius=0.65 
deltar=normrnd(0,0.65);   %standard deviation of radius=0.65 
% deltar=0;
r=r0+deltar;
deltaL=-deltar*20/3.2;  %AM paper
L=L0+deltaL; 

new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r)^2+h0^2); %countor length=n_turns*sqrt((2pir)^2+h^2)
new_b=new_turns*2*pi/19+0.0005; %the range of b.*n_bp+c1 = how much 2pi = how much turns
new_c1_x=b+c1_x-new_b;  %keep the first position for C1a(1) (x1,y1,z1) same
new_c1_y=b+c1_y-new_b;
new_d=new_turns*h0/19;  
new_e1=d+e1-new_d;
new_c2_x=b+c2_x-new_b;
new_c2_y=b+c2_y-new_b;
new_e2=20*d-20*new_d+e2+deltaL; %the first position for C1b(20) is (x2,y2,z2+deltaL)

% new helix.A
new_x1=r*sin(new_b.*n_bp+new_c1_x);
new_y1=r*sin(new_b*n_bp+new_c1_y);
new_z1=new_d.*n_bp+new_e1;  

% new helix.B
new_x2=r*sin(new_b.*n_bp+new_c2_x);
new_y2=r*sin(new_b.*n_bp+new_c2_y);
new_z2=new_d.*n_bp+new_e2;


% %%%%% atan2(norm(cross(vector1,vector2)), dot(vector1,vector2)) %angle cal
% %%%%% vector2_proj=(dot(vector2,z_axis)/norm(z_axis))*z_axis; %projection
% 
% %



%%
   case 'A'
% Model A
% % %parameter
deltah=normrnd(0,3.85);   %standard deviation of pitch height= 3.85 
% deltah=normrnd(0,sigma_h);   %standard deviation of pitch height= 3.85 
% deltah=0;
h=h0+deltah;
deltaL=deltah*20/18.2;  %AM paper
L=L0+deltaL; 
r=r0;
new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r0)^2+h^2); %countor length=n_turns*sqrt((2pir)^2+h^2)
new_b=new_turns*2*pi/19+0.0005; %the range of b.*n_bp+c1 = how much 2pi = how much turns; -0.0138: when deltar=0,keep new_b = b
new_c1_x=b+c1_x-new_b;  %keep the first position for C1a(1) (x1,y1,z1) same
new_c1_y=b+c1_y-new_b;
new_d=new_turns*h/19;  
new_e1=d+e1-new_d;
new_c2_x=b+c2_x-new_b;
new_c2_y=b+c2_y-new_b;
new_e2=20*d-20*new_d+e2+deltaL; %the first position for C1b(20) is (x2,y2,z2+deltaL)

% new helix.A
new_x1=r*sin(new_b.*n_bp+new_c1_x);
new_y1=r*sin(new_b*n_bp+new_c1_y);
new_z1=new_d.*n_bp+new_e1;  

% new helix.B
new_x2=r*sin(new_b.*n_bp+new_c2_x);
new_y2=r*sin(new_b.*n_bp+new_c2_y);
new_z2=new_d.*n_bp+new_e2;
% % 
    otherwise
end 

%%%%new coordinate calculation
%new C1a/C1b position
C1a=[new_x1;new_y1;new_z1]';
C1b=[new_x2;new_y2;new_z2]';

%vector calculation 
x_axis1_Spin1=(C1b(3,:)-C1a(3,:))/norm(C1b(3,:)-C1a(3,:)); 
z_proj=(dot(z_axis,x_axis1_Spin1)/norm(x_axis1_Spin1)^2)*x_axis1_Spin1;
z_axis_Spin1=z_axis-z_proj;
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis_Spin1=cross(x_axis1_Spin1,z_axis);
y_axis_Spin1=y_axis_Spin1/norm(y_axis_Spin1);

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)+y_axis_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%%Kugelverteilung
% theta1=rand(1,1)*pi;
% phi1=rand(1,1)*2*pi;
% M(k,:)=Spin_1+[sin(theta1)*cos(phi1) sin(theta1)*sin(phi1) cos(theta1)];
% M(k,:)=Spin_1;

%C Label position rechnen
Z_Spin1(k,:)=cross(x_axis3_Spin1,x_axis1_Spin1);
Z_Spin1(k,:)=Z_Spin1(k,:)/norm(Z_Spin1(k,:));
y_Spin1(k,:)=cross(x_axis1_Spin1,Z_Spin1(k,:));
y_Spin1(k,:)=y_Spin1(k,:)/norm(y_Spin1(k,:));
N1_1(k,:)=C1a(3,:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin1+1.5*sin(54.2/360*2*pi)*y_Spin1(k,:);
C2_1(k,:)=C1a(3,:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin1+2.5*sin(24.3/360*2*pi)*y_Spin1(k,:);
C4_1(k,:)=C1a(3,:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin1+4.2*sin(52.6/360*2*pi)*y_Spin1(k,:);
C5_1(k,:)=C1a(3,:)+3.8*cos(72.2/360*2*pi)*x_axis1_Spin1+3.8*sin(72.2/360*2*pi)*y_Spin1(k,:);
C45_1(k,:)=(C4_1(k,:)+C5_1(k,:))/2;
N1C2_1(k,:)=(N1_1(k,:)+C2_1(k,:))/2;
X_Spin1(k,:)=C5_1(k,:)-N1_1(k,:);
X_Spin1(k,:)=X_Spin1(k,:)/norm(X_Spin1(k,:));
Y_Spin1(k,:)=cross(Z_Spin1(k,:),X_Spin1(k,:));
Y_Spin1(k,:)=Y_Spin1(k,:)/norm(Y_Spin1(k,:));

% % %%rotation around N-O axis 
%%first rotation around y-axis with 6 grad
rotate_theta1_1(k,:)=normrnd(0,sigma_y); %range -5 to 5grad
[X1(k,:),Y1(k,:),Z1(k,:)] = AxisAngleRotate(X_Spin1(k,:),Y_Spin1(k,:),Z_Spin1(k,:),Y_Spin1(k,:),rotate_theta1_1(k,:));
%second around z-axis with 5 grad
rotate_theta2_1(k,:)=normrnd(0,sigma_z); %range -5 to 5grad
[X2(k,:),Y2(k,:),Z2(k,:)] = AxisAngleRotate(X1(k,:),Y1(k,:),Z1(k,:),Z1(k,:),rotate_theta2_1(k,:));
%second coordinatesystem for Cdot Spin label
[xaxis1_dre(k,:),yaxis1_dre(k,:),zaxis1_dre(k,:)] = AxisAngleRotate(X2(k,:),Y2(k,:),Z2(k,:),Z2(k,:),27);
%find electron center
M(k,:)=N1C2_1(k,:)+X2(k,:)*(2.3+0.688*1.5)+xaxis1_dre(k,:)*(3.959*1.5+0.7);


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

x_axis1_Spin2=(C1a(7+(nd-1),:)-C1b(7+(nd-1),:))/norm(C1a(7+(nd-1),:)-C1b(7+(nd-1),:));
z_proj_Spin2=(dot(z_axis,x_axis1_Spin2)/norm(x_axis1_Spin2)^2)*x_axis1_Spin2;
z_axis_Spin2=z_axis-z_proj_Spin2;
z_axis_Spin2=z_axis_Spin2/norm(z_axis_Spin2);
y_axis_Spin2=cross(x_axis1_Spin2,z_axis);
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);

%%first rotation
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))-y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))-z_axis_Spin2.*sin(beta_2(nd));
%%find Spin2
Spin_2=C1b(7+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

%%Kugelverteilung
% theta2=rand(1,1)*pi;
% phi2=rand(1,1)*2*pi;
% M2(k,:)=Spin_2+[sin(theta2)*cos(phi2) sin(theta2)*sin(phi2) cos(theta2)];
% M2(k,:)=Spin_2;

%C Label position rechnen
Z_Spin2(k,:)=cross(x_axis3_Spin2,x_axis1_Spin2);
Z_Spin2(k,:)=Z_Spin2(k,:)/norm(Z_Spin2(k,:));
y_Spin2(k,:)=cross(x_axis1_Spin2,Z_Spin2(k,:));
y_Spin2(k,:)=y_Spin2(k,:)/norm(y_Spin2(k,:));
N1_2(k,:)=C1b(7+(nd-1),:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin2+1.5*sin(54.2/360*2*pi)*y_Spin2(k,:);
C2_2(k,:)=C1b(7+(nd-1),:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin2+2.5*sin(24.3/360*2*pi)*y_Spin2(k,:);
C4_2(k,:)=C1b(7+(nd-1),:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin2+4.2*sin(52.6/360*2*pi)*y_Spin2(k,:);
C5_2(k,:)=C1b(7+(nd-1),:)+3.8*cos(72.9/360*2*pi)*x_axis1_Spin2+3.8*sin(72.9/360*2*pi)*y_Spin2(k,:);
C45_2(k,:)=(C4_2(k,:)+C5_2(k,:))/2;
N1C2_2(k,:)=(N1_2(k,:)+C2_2(k,:))/2;
X_Spin2(k,:)=C5_2(k,:)-N1_2(k,:);
X_Spin2(k,:)=X_Spin2(k,:)/norm(X_Spin2(k,:));
Y_Spin2(k,:)=cross(Z_Spin2(k,:),X_Spin2(k,:));
Y_Spin2(k,:)=Y_Spin2(k,:)/norm(Y_Spin2(k,:));

% % % %%rotation around N-O axis 
%%first rotation around y-axis with 6 grad
rotate_theta1_2(k,:)=normrnd(0,sigma_y); %range -5 to 5grad
[X1_2(k,:),Y1_2(k,:),Z1_2(k,:)] = AxisAngleRotate(X_Spin2(k,:),Y_Spin2(k,:),Z_Spin2(k,:),Y_Spin2(k,:),rotate_theta1_2(k,:));
%second around z-axis with 5 grad
rotate_theta2_2(k,:)=normrnd(0,sigma_z); %range -5 to 5gradg
[X2_2(k,:),Y2_2(k,:),Z2_2(k,:)] = AxisAngleRotate(X1_2(k,:),Y1_2(k,:),Z1_2(k,:),Z1_2(k,:),rotate_theta2_2(k,:));
%second coordinatesystem for Cdot Spin label
[xaxis2_dre(k,:),yaxis2_dre(k,:),zaxis2_dre(k,:)] = AxisAngleRotate(X2_2(k,:),Y2_2(k,:),Z2_2(k,:),Z2_2(k,:),27);
%find electron center
M2(k,:)=N1C2_2(k,:)+X2_2(k,:)*(2.3+0.688*1.5)+xaxis2_dre(k,:)*(3.959*1.5+0.7);

end 
R=M2-M;

[o1,o2,r]=RuCoordv1(R,xaxis1_dre,yaxis1_dre,zaxis1_dre,xaxis2_dre,yaxis2_dre,zaxis2_dre);

Conformers.M=M;
Conformers.M2=M2;
% Conformers.R=R;
Conformers.EulerAngles.R1=[o1(:,:)];
% Conformers.EulerAngles.R1(:,2)=Conformers.EulerAngles.R1(:,2)+ones(length(Conformers.EulerAngles.R1(:,2)),1)*deg2rad(10);
Conformers.EulerAngles.R2=[o2(:,:)];
Conformers.Distance = r/10;
% pd=fitdist(Conformers.Distance,'Normal');
% R_mean=pd.mu;
% sigma=pd.sigma;
zeiten = zeit*1000;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the program
% mean_trend=mean(trend);
h=histfit(Conformers.Distance);
h(1).FaceColor = [.3 .75 .93];
h(2).Color = [.0 .0 1];
end 
