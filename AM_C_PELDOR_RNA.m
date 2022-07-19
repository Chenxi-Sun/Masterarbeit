% function [Result,R_mean,FWHM] = AM_C_PELDOR_RNA(sigma_y,nd,EP,zeit)
% 
nd=9;
nd=nd-7;
% sigma_y=0;
sigma_z=0;

%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

% %parameter for expression for C1-points at DNA from pymol
r0=8.284;
h0=29.8476;
b=0.5709;
c1_x=3.468;
c1_y=1.908;
c2_x=-1.457;
c2_y=3.265;
d=2.712;
e1=-1.527;
e2=-3.944;

%helix.A
x1=r0*sin(b.*n_bp+c1_x);
y1=r0*sin(b.*n_bp+c1_y);
z1=d.*n_bp+e1;  

%helix.B
x2=r0*sin(b.*n_bp+c2_x);
y2=r0*sin(b.*n_bp+c2_y);
z2=d.*n_bp+e2;


%z_axis
z_axis=[0 0 1];

%drehen die koordinatensystem
Z_axis=[0 0 1];
X_axis=[1 0 0];
Y_axis=[0 1 0];

%%rotationsangle 1st SL
alpha_1=74.42/360*2*pi;
beta_1=5.2/360*2*pi; 
%%rotationsangle 2nd SL 8-15 bp; (5-7)samples dont exist
alpha_2=[74.4187 74.6188 76.6019 78.0654 77.2759 75.1897 74.2406 75.5165]./360.*2.*pi;
beta_2=[5.2154 5.0264 4.3064 3.6792 3.6413 4.0283 4.3597 4.3319]./360.*2.*pi;

%bending direction
for k=1:361
theta=(k-1);
[new_X, new_Y, new_Z]=AxisAngleRotate(X_axis,Y_axis,Z_axis,Z_axis,theta);
C1a_dre=[new_X;new_Y;new_Z]\[x1;y1;z1];
C1a_dre=C1a_dre';
x1_dre=C1a_dre(1:20,1);
y1_dre=C1a_dre(1:20,2);
z1_dre=C1a_dre(1:20,3);
C1b_dre=[new_X;new_Y;new_Z]\[x2;y2;z2];
C1b_dre=C1b_dre';
x2_dre=C1b_dre(1:20,1);
y2_dre=C1b_dre(1:20,2);
z2_dre=C1b_dre(1:20,3);

%bending angle
alpha=28/360*2*pi;
R1=(z1(20)-z1(1))/alpha;
R2=(z2(20)-z2(1))/alpha;

for i=1:20
phi_1(i)=(z1_dre(i)-z1_dre(1))/(z1_dre(20)-z1_dre(1))*alpha;
phi_2(i)=(z2_dre(i)-z2_dre(1))/(z2_dre(20)-z2_dre(1))*alpha;
new_x1(i)=x1_dre(i)*cos(phi_1(i));
new_z1(i)=(R1+x1_dre(i))*sin(phi_1(i))+z1(1);
new_x2(i)=x2_dre(i)*cos(phi_2(i));
new_z2(i)=(R2+x2_dre(i))*sin(phi_2(i))+z1(1);
end

%coordiante after bending 
new_y1=y1_dre;
new_y2=y2_dre;

C1a(:,1)=new_x1;C1a(:,2)=new_y1;C1a(:,3)=new_z1;
C1b(:,1)=new_x2;C1b(:,2)=new_y2;C1b(:,3)=new_z2;

%vector calculation 
x_axis1_Spin1=(C1b(3,:)-C1a(3,:))/norm(C1b(3,:)-C1a(3,:)); 
z_proj=(dot(z_axis,x_axis1_Spin1)/norm(x_axis1_Spin1)^2)*x_axis1_Spin1;
z_axis_Spin1=z_axis-z_proj;
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis_Spin1=cross(x_axis1_Spin1,z_axis_Spin1);
y_axis_Spin1=y_axis_Spin1/norm(y_axis_Spin1);

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)-y_axis_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)-z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A


%C Label position rechnen
Z_Spin1(k,:)=cross(x_axis3_Spin1,x_axis1_Spin1);
Z_Spin1(k,:)=Z_Spin1(k,:)/norm(Z_Spin1(k,:));
y_Spin1(k,:)=cross(x_axis1_Spin1,Z_Spin1(k,:));
y_Spin1(k,:)=y_Spin1(k,:)/norm(y_Spin1(k,:));
N1_1(k,:)=C1a(3,:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin1+1.5*sin(54.2/360*2*pi)*y_Spin1(k,:);
C2_1(k,:)=C1a(3,:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin1+2.5*sin(24.3/360*2*pi)*y_Spin1(k,:);
C4_1(k,:)=C1a(3,:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin1+4.2*sin(52.6/360*2*pi)*y_Spin1(k,:);
C5_1(k,:)=C1a(3,:)+3.8*cos(72.2/360*2*pi)*x_axis1_Spin1+3.8*sin(72.2/360*2*pi)*y_Spin1(k,:);
% C45_1(k,:)=(C4_1(k,:)+C5_1(k,:))/2;
N1C2_1(k,:)=(N1_1(k,:)+C2_1(k,:))/2;
X_Spin1(k,:)=C5_1(k,:)-N1_1(k,:);
X_Spin1(k,:)=X_Spin1(k,:)/norm(X_Spin1(k,:));
Y_Spin1(k,:)=cross(Z_Spin1(k,:),X_Spin1(k,:));
Y_Spin1(k,:)=Y_Spin1(k,:)/norm(Y_Spin1(k,:));

% % %%rotation around N-O axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_1(k,:)=normrnd(0,sigma_y); 
[X1(k,:), Y1(k,:), Z1(k,:)] = AxisAngleRotate(X_Spin1(k,:),Y_Spin1(k,:),Z_Spin1(k,:),Y_Spin1(k,:),rotate_theta1_1(k,:));
%second rotation around z-axis with 5 grad
rotate_theta2_1(k,:)=normrnd(0,sigma_z); 
[X2(k,:), Y2(k,:), Z2(k,:)] = AxisAngleRotate(X1(k,:),Y1(k,:),Z1(k,:),Z1(k,:),rotate_theta2_1(k,:));
%find electron center
M(k,:)=N1C2_1(k,:)+X2(k,:)*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

x_axis1_Spin2=(C1a(10+(nd-1),:)-C1b(10+(nd-1),:))/norm(C1a(10+(nd-1),:)-C1b(10+(nd-1),:));
z_proj_Spin2=(dot(z_axis,x_axis1_Spin2)/norm(x_axis1_Spin2)^2)*x_axis1_Spin2;
z_axis_Spin2=z_axis-z_proj_Spin2;
z_axis_Spin2=z_axis_Spin2/norm(z_axis_Spin2);
y_axis_Spin2=cross(x_axis1_Spin2,z_axis_Spin2);
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);

%%first rotation
% x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_1)+y_axis_Spin2.*sin(alpha_1);
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))+y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
% x_axis3_Spin2=x_axis2_Spin2.*cos(beta_1)+z_axis_Spin2.*sin(beta_1);
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))+z_axis_Spin2.*sin(beta_2(nd));
%%find Spin2
Spin_2=C1b(10+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A



%Cm Label position rechnen
Z_Spin2(k,:)=cross(x_axis3_Spin2,x_axis1_Spin2);
Z_Spin2(k,:)=Z_Spin2(k,:)/norm(Z_Spin2(k,:));
y_Spin2(k,:)=cross(x_axis1_Spin2,Z_Spin2(k,:));
y_Spin2(k,:)=y_Spin2(k,:)/norm(y_Spin2(k,:));
N1_2(k,:)=C1b(10+(nd-1),:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin2+1.5*sin(54.2/360*2*pi)*y_Spin2(k,:);
C2_2(k,:)=C1b(10+(nd-1),:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin2+2.5*sin(24.3/360*2*pi)*y_Spin2(k,:);
C4_2(k,:)=C1b(10+(nd-1),:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin2+4.2*sin(52.6/360*2*pi)*y_Spin2(k,:);
C5_2(k,:)=C1b(10+(nd-1),:)+3.8*cos(72.9/360*2*pi)*x_axis1_Spin2+3.8*sin(72.9/360*2*pi)*y_Spin2(k,:);
C45_2(k,:)=(C4_2(k,:)+C5_2(k,:))/2;
N1C2_2(k,:)=(N1_2(k,:)+C2_2(k,:))/2;
X_Spin2(k,:)=C5_2(k,:)-N1_2(k,:);
X_Spin2(k,:)=X_Spin2(k,:)/norm(X_Spin2(k,:));
Y_Spin2(k,:)=cross(Z_Spin2(k,:),X_Spin2(k,:));
Y_Spin2(k,:)=Y_Spin2(k,:)/norm(Y_Spin2(k,:));

% % % %%rotation around N-O axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_2(k,:)=normrnd(0,sigma_y); %range -5 to 5grad
[X1_2(k,:), Y1_2(k,:), Z1_2(k,:)] = AxisAngleRotate(X_Spin2(k,:),Y_Spin2(k,:),Z_Spin2(k,:),Y_Spin2(k,:),rotate_theta1_2(k,:));
%second rotaion around z-axis with 5 grad
rotate_theta2_2(k,:)=normrnd(0,sigma_z); %range -5 to 5grad
[X2_2(k,:), Y2_2(k,:), Z2_2(k,:)] = AxisAngleRotate(X1_2(k,:),Y1_2(k,:),Z1_2(k,:),Z1_2(k,:),rotate_theta2_2(k,:));
%find electron center
M2(k,:)=N1C2_2(k,:)+X2_2(k,:)*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);

end 
R=M2-M;
[o1,o2,r]=RuCoordv1(R,X2,Y2,Z2,X2_2,Y2_2,Z2_2);

Conformers.M=M;
Conformers.M2=M2;
Conformers.EulerAngles.R1=[o1(:,:)];
Conformers.EulerAngles.R2=[o2(:,:)];
Conformers.Distance = r/10;
pd=fitdist(Conformers.Distance,'Normal'); %fit distances to a normal distribution
R_mean=pd.mu;
sigma=pd.sigma;  %sigma value of normal distribution
FWHM=2.3548*sigma;  %FWHM 
% zeiten = zeit*1000;
% Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the programr
% Result = MainPELDORtime_modAC(EP,Conformers,zeiten,6.5); %for G-band
% end 
