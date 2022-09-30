function [Result,ymax] = AM_C_PELDOR_CdotDNA(sigma_y,nd,EP,zeit)
% function [Result,R_mean,FWHM] = AM_C_PELDOR_DNA(sigma_y,nd,EP,zeit)
% function [R_mean,FWHM] = AM_C_PELDOR_DNA(sigma_y,nd,EP,zeit)
% function [Result] = AM_C_PELDOR_DNA(sigma_y,nd,EP,zeit)
% 
% for nd=5:14
nd=9;
nd=nd-4;
sigma_y=6;
sigma_z=6;

%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

% %parameter for expression for C1-points at DNA from pymol
r0=5.68;
h0=33.7242; 
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


%z_axis
z_axis=[0 0 1];

%drehen die koordinatensystem
Z_axis=[0 0 1];
X_axis=[1 0 0];
Y_axis=[0 1 0];

%%rotationsangle 1st SL
alpha_1=76.20/360*2*pi;
beta_1=6.64/360*2*pi;
%%rotationsangle 2nd SL 5-14 bp;
alpha_2=[76.1446 76.1787 76.1766 76.2030 76.1871 76.1902 76.1983 76.1741 76.1634 76.1681]./360.*2.*pi;
beta_2=[6.6477 6.6232 6.6401 6.6957 6.7435 6.7840 6.7882 6.7712 6.7300 6.6836]./360.*2.*pi;

%bending direction
for w=1:36
theta=(w-1)*10;
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
for k=1:25
alpha=abs(normrnd(0,28)/360*2*pi);
% alpha=28/360*2*pi;
if alpha == 0
    R1=0;R2=0;
else 
    R1=(z1(20)-z1(1))/alpha;  %%Radius value of torus
    R2=(z2(20)-z2(1))/alpha;
end 

for i=1:20
phi_1(i)=(z1_dre(i)-z1_dre(1))/(z1_dre(20)-z1_dre(1))*alpha;
phi_2(i)=(z2_dre(i)-z2_dre(1))/(z2_dre(20)-z2_dre(1))*alpha;
new_x1(i)=x1_dre(i)*cos(phi_1(i));  %(sin(alpha)/abs(sin(alpha))) for sign
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
z_axis_Spin1=z_axis-z_proj;   %Z_loc=Z_glob-Z_proj
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis_Spin1=cross(x_axis1_Spin1,z_axis_Spin1);
y_axis_Spin1=y_axis_Spin1/norm(y_axis_Spin1);

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)+y_axis_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%C Label position rechnen
Z_Spin1(k,:,w)=cross(x_axis3_Spin1,x_axis1_Spin1);
Z_Spin1(k,:,w)=Z_Spin1(k,:,w)/norm(Z_Spin1(k,:,w));
y_Spin1(k,:,w)=cross(x_axis1_Spin1,Z_Spin1(k,:,w));
y_Spin1(k,:,w)=y_Spin1(k,:,w)/norm(y_Spin1(k,:,w));
N1_1(k,:,w)=C1a(3,:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin1+1.5*sin(54.2/360*2*pi)*y_Spin1(k,:,w);
C2_1(k,:,w)=C1a(3,:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin1+2.5*sin(24.3/360*2*pi)*y_Spin1(k,:,w);
C4_1(k,:,w)=C1a(3,:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin1+4.2*sin(52.6/360*2*pi)*y_Spin1(k,:,w);
C5_1(k,:,w)=C1a(3,:)+3.8*cos(72.2/360*2*pi)*x_axis1_Spin1+3.8*sin(72.2/360*2*pi)*y_Spin1(k,:,w);
% C45_1(k,:,w)=(C4_1(k,:,w)+C5_1(k,:,w))/2;
N1C2_1(k,:,w)=(N1_1(k,:,w)+C2_1(k,:,w))/2;
X_Spin1(k,:,w)=C5_1(k,:,w)-N1_1(k,:,w);
X_Spin1(k,:,w)=X_Spin1(k,:,w)/norm(X_Spin1(k,:,w));
Y_Spin1(k,:,w)=cross(Z_Spin1(k,:,w),X_Spin1(k,:,w));
Y_Spin1(k,:,w)=Y_Spin1(k,:,w)/norm(Y_Spin1(k,:,w));

% % %%rotation around molecular axis 
%first rotaion around y-axis
rotate_theta1_1(k,:,w)=normrnd(0,sigma_y); 
[X1(k,:,w), Y1(k,:,w), Z1(k,:,w)] = AxisAngleRotate(X_Spin1(k,:,w),Y_Spin1(k,:,w),Z_Spin1(k,:,w),Y_Spin1(k,:,w),rotate_theta1_1(k,:,w));
%second rotation around z-axis
rotate_theta2_1(k,:,w)=normrnd(0,sigma_z); 
[X2(k,:,w), Y2(k,:,w) Z2(k,:,w)] = AxisAngleRotate(X1(k,:,w),Y1(k,:,w),Z1(k,:,w),Z1(k,:,w),rotate_theta2_1(k,:,w));
%second coordinatesystem for Cdot Spin label
[xaxis1_dre(k,:,w),yaxis1_dre(k,:,w),zaxis1_dre(k,:,w)] = AxisAngleRotate(X2(k,:,w),Y2(k,:,w),Z2(k,:,w),Z2(k,:,w),27);
%find electron center
M(k,:,w)=N1C2_1(k,:,w)+X2(k,:,w)*(2.3+0.688*1.5)+xaxis1_dre(k,:,w)*(3.959*1.5+0.7);

%%Spin2 position 
x_axis1_Spin2=(C1a(7+(nd-1),:)-C1b(7+(nd-1),:))/norm(C1a(7+(nd-1),:)-C1b(7+(nd-1),:));
z_proj_Spin2=(dot(z_axis,x_axis1_Spin2)/norm(x_axis1_Spin2)^2)*x_axis1_Spin2;
z_axis_Spin2=z_axis-z_proj_Spin2;
z_axis_Spin2=z_axis_Spin2/norm(z_axis_Spin2);
y_axis_Spin2=cross(x_axis1_Spin2,z_axis_Spin2);
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);

%%first rotation
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))-y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))-z_axis_Spin2.*sin(beta_2(nd));
%%find Spin2
Spin_2=C1b(7+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

%C Label position rechnen
Z_Spin2(k,:,w)=cross(x_axis3_Spin2,x_axis1_Spin2);
Z_Spin2(k,:,w)=Z_Spin2(k,:,w)/norm(Z_Spin2(k,:,w));
y_Spin2(k,:,w)=cross(x_axis1_Spin2,Z_Spin2(k,:,w));
y_Spin2(k,:,w)=y_Spin2(k,:,w)/norm(y_Spin2(k,:,w));
N1_2(k,:,w)=C1b(7+(nd-1),:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin2+1.5*sin(54.2/360*2*pi)*y_Spin2(k,:,w);
C2_2(k,:,w)=C1b(7+(nd-1),:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin2+2.5*sin(24.3/360*2*pi)*y_Spin2(k,:,w);
C4_2(k,:,w)=C1b(7+(nd-1),:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin2+4.2*sin(52.6/360*2*pi)*y_Spin2(k,:,w);
C5_2(k,:,w)=C1b(7+(nd-1),:)+3.8*cos(72.9/360*2*pi)*x_axis1_Spin2+3.8*sin(72.9/360*2*pi)*y_Spin2(k,:,w);
C45_2(k,:,w)=(C4_2(k,:,w)+C5_2(k,:,w))/2;
N1C2_2(k,:,w)=(N1_2(k,:,w)+C2_2(k,:,w))/2;
X_Spin2(k,:,w)=C5_2(k,:,w)-N1_2(k,:,w);
X_Spin2(k,:,w)=X_Spin2(k,:,w)/norm(X_Spin2(k,:,w));
Y_Spin2(k,:,w)=cross(Z_Spin2(k,:,w),X_Spin2(k,:,w));
Y_Spin2(k,:,w)=Y_Spin2(k,:,w)/norm(Y_Spin2(k,:,w));

% % % %%rotation around molecular axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_2(k,:,w)=normrnd(0,sigma_y); 
[X1_2(k,:,w), Y1_2(k,:,w), Z1_2(k,:,w)] = AxisAngleRotate(X_Spin2(k,:,w),Y_Spin2(k,:,w),Z_Spin2(k,:,w),Y_Spin2(k,:,w),rotate_theta1_2(k,:,w));
%second rotaion around z-axis with 5 grad
rotate_theta2_2(k,:,w)=normrnd(0,sigma_z); 
[X2_2(k,:,w), Y2_2(k,:,w), Z2_2(k,:,w)] = AxisAngleRotate(X1_2(k,:,w),Y1_2(k,:,w),Z1_2(k,:,w),Z1_2(k,:,w),rotate_theta2_2(k,:,w));
%second coordinatesystem for Cdot Spin label
[xaxis2_dre(k,:,w),yaxis2_dre(k,:,w),zaxis2_dre(k,:,w)] = AxisAngleRotate(X2_2(k,:,w),Y2_2(k,:,w),Z2_2(k,:,w),Z2_2(k,:,w),27);
%find electron center
M2(k,:,w)=N1C2_2(k,:,w)+X2_2(k,:,w)*(2.3+0.688*1.5)+xaxis2_dre(k,:,w)*(3.959*1.5+0.7);

end 

end 
M2_2d=M2(:,:,1);
M_2d=M(:,:,1);
X2_2d=X2(:,:,1);
Y2_2d=Y2(:,:,1);
Z2_2d=Z2(:,:,1);
X2_2_2d=X2_2(:,:,1);
Y2_2_2d=Y2_2(:,:,1);
Z2_2_2d=Z2_2(:,:,1);

for s=1:w-1
M2_2d=vertcat(M2_2d,M2(:,:,s+1));M_2d=vertcat(M_2d,M(:,:,s+1));
X2_2d=vertcat(X2_2d,X2(:,:,s+1));Y2_2d=vertcat(Y2_2d,Y2(:,:,s+1));Z2_2d=vertcat(Z2_2d,Z2(:,:,s+1));
X2_2_2d=vertcat(X2_2_2d,X2_2(:,:,s+1));Y2_2_2d=vertcat(Y2_2_2d,Y2_2(:,:,s+1));Z2_2_2d=vertcat(Z2_2_2d,Z2_2(:,:,s+1));
end 


R=M2_2d-M_2d;
[o1,o2,r]=RuCoordv1(R,X2_2d,Y2_2d,Z2_2d,X2_2_2d,Y2_2_2d,Z2_2_2d);

Conformers.M=M_2d;
Conformers.M2=M2_2d;
Conformers.EulerAngles.R1=[o1(:,:)];
Conformers.EulerAngles.R2=[o2(:,:)];
Conformers.Distance = r/10;
pd=fitdist(Conformers.Distance,'Normal'); %fit distances to a normal distribution
R_mean=pd.mu;
sigma=pd.sigma;  %sigma value of normal distribution
FWHM=2.3548*sigma;  %FWHM 
zeiten = zeit*1000;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the programr
% Result = MainPELDORtime_modAC(EP,Conformers,zeiten,6.5); %for G-band
h=histfit(Conformers.Distance);   %plot distribution
h(1).FaceColor = [0.43921568627451 0.631372549019608 0.76078431372549];
h(2).Color = [0 0.447058823529412 0.741176470588235];
y=get(h,'YData');  
y1=cell2mat(y(1));
ymax=max(y1);     %%max value of histogramm
% Distance(nd,:)=r/10;
end 