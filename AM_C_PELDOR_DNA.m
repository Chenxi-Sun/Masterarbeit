function [Result] = AM_C_PELDOR_DNA(nd,EP,zeit)

nd=9;
nd=nd-4;
sigma_y=0;
sigma_z=0;

%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

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

%parameter for expression for C1-points at DNA for AM
% r0=5.85;
% h0=33;
% b=0.6253;
% c1_x=-1.6875;
% c1_y=-0.1261;
% c2_x=0.7488;
% c2_y=2.3275;
% d=3.2814;
% e1=-69.1064;
% e2=-69.1862;

%helix.A
x1=r0*sin(b.*n_bp+c1_x);
y1=r0*sin(b.*n_bp+c1_y);
z1=d.*n_bp+e1;  
C1a=[x1;y1;z1]';
%helix.B
x2=r0*sin(b.*n_bp+c2_x);
y2=r0*sin(b.*n_bp+c2_y);
z2=d.*n_bp+e2;
C1b=[x2;y2;z2]';

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
for w=1:20
theta=(w-1)*20/360*360;
[new_X new_Y new_Z]=AxisAngleRotate(X_axis,Y_axis,Z_axis,Z_axis,theta);
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

%%intial position after determing bending direction
%vector calculation 
x_axis1_Spin1=(C1b_dre(3,:)-C1a_dre(3,:))/norm(C1b_dre(3,:)-C1a_dre(3,:)); 
z_proj=(dot(z_axis,x_axis1_Spin1)/norm(x_axis1_Spin1)^2)*x_axis1_Spin1;
z_axis_Spin1=z_axis-z_proj;
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis_Spin1=cross(x_axis1_Spin1,z_axis_Spin1);
y_axis_Spin1=y_axis_Spin1/norm(y_axis_Spin1);

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)+y_axis_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a_dre(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%C Label position rechnen
Z_Spin1=cross(x_axis3_Spin1,x_axis1_Spin1);
Z_Spin1=Z_Spin1/norm(Z_Spin1);
y_Spin1=cross(x_axis1_Spin1,Z_Spin1);
y_Spin1=y_Spin1/norm(y_Spin1);
N1_1=C1a(3,:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin1+1.5*sin(54.2/360*2*pi)*y_Spin1;
C2_1=C1a(3,:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin1+2.5*sin(24.3/360*2*pi)*y_Spin1;
C4_1=C1a(3,:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin1+4.2*sin(52.6/360*2*pi)*y_Spin1;
C5_1=C1a(3,:)+3.8*cos(72.2/360*2*pi)*x_axis1_Spin1+3.8*sin(72.2/360*2*pi)*y_Spin1;
N1C2_1=(N1_1+C2_1)/2;
X_Spin1=C5_1-N1_1;
X_Spin1=X_Spin1/norm(X_Spin1);
Y_Spin1=cross(Z_Spin1,X_Spin1);
Y_Spin1=Y_Spin1/norm(Y_Spin1);

%%Spin2 position 
x_axis1_Spin2=(C1a_dre(7+(nd-1),:)-C1b_dre(7+(nd-1),:))/norm(C1a_dre(7+(nd-1),:)-C1b_dre(7+(nd-1),:));
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
Spin_2=C1b_dre(7+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

%C Label position rechnen
Z_Spin2=cross(x_axis3_Spin2,x_axis1_Spin2);
Z_Spin2=Z_Spin2/norm(Z_Spin2);
y_Spin2=cross(x_axis1_Spin2,Z_Spin2);
y_Spin2=y_Spin2/norm(y_Spin2);
N1_2=C1b(7+(nd-1),:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin2+1.5*sin(54.2/360*2*pi)*y_Spin2;
C2_2=C1b(7+(nd-1),:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin2+2.5*sin(24.3/360*2*pi)*y_Spin2;
C4_2=C1b(7+(nd-1),:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin2+4.2*sin(52.6/360*2*pi)*y_Spin2;
C5_2=C1b(7+(nd-1),:)+3.8*cos(72.9/360*2*pi)*x_axis1_Spin2+3.8*sin(72.9/360*2*pi)*y_Spin2;
C45_2=(C4_2+C5_2)/2;
N1C2_2=(N1_2+C2_2)/2;
X_Spin2=C5_2-N1_2;
X_Spin2=X_Spin2/norm(X_Spin2);
Y_Spin2=cross(Z_Spin2,X_Spin2);
Y_Spin2=Y_Spin2/norm(Y_Spin2);

%middle axis cross the zylinder
for k=1:25
delta_alpha=normrnd(0,28)/360*2*pi;

m_axisA_x=zeros(1,20);
m_axisA_y=zeros(1,20);
m_axisA_z=z1;

m_axisB_x=zeros(1,20);
m_axisB_y=zeros(1,20);
m_axisB_z=z2;

%%after bending   
if delta_alpha>0   %>0 means bending in left side
    r_bend=(z1(20)-z1(1))/delta_alpha;  %%2pi*r/2pi=middle axis length=z1(20)-z1(1)
    middle_punkt_A=[0 0 (z1(20)+z1(1))/2];  %middle point of the middle axis 
    middle_punkt_B=[0 0 (z2(20)+z2(1))/2];
    null_punkt_A=[-r_bend 0 (z1(20)+z1(1))/2]; %middle point of the circle
    null_punkt_B=[-r_bend 0 (z2(20)+z2(1))/2];
    
    for s=1:20
    Beta(s)=(z1(s)-z1(1))/r_bend;   %
    m_axisA_x(s)=-(r_bend-r_bend*cos(Beta(s)-delta_alpha/2));
    m_axisA_y(s)=0;
    m_axisA_z(s)=r_bend*sin(Beta(s)-delta_alpha/2)+((z1(20)-z1(1))/2+z1(1));
    m_axisA(s,:)=[m_axisA_x(s) m_axisA_y(s) m_axisA_z(s)];
    
    m_axisB_x(s)=-(r_bend-r_bend*cos(Beta(s)-delta_alpha/2));
    m_axisB_y(s)=0;
    m_axisB_z(s)=r_bend*sin(Beta(s)-delta_alpha/2)+((z2(20)-z2(1))/2+z2(1));
    m_axisB(s,:)=[m_axisB_x(s) m_axisB_y(s) m_axisB_z(s)];
    
    new_x_vectorA(s,:)=m_axisA(s,:)-null_punkt_A;
    new_x_vectorA(s,:)=new_x_vectorA(s,:)/norm(new_x_vectorA(s,:));
    new_x_vectorB(s,:)=m_axisB(s,:)-null_punkt_B;
    new_x_vectorB(s,:)=new_x_vectorB(s,:)/norm(new_x_vectorB(s,:));
    
    m_proj_A(s,:)=(dot(z_axis,new_x_vectorA(s,:))/norm(new_x_vectorA(s,:))^2)*new_x_vectorA(s,:);
    new_zaxisA(s,:)=z_axis-m_proj_A(s,:);

    m_proj_B(s,:)=(dot(z_axis,new_x_vectorB(s,:))/norm(new_x_vectorB(s,:))^2)*new_x_vectorB(s,:);
    new_zaxisB(s,:)=z_axis-m_proj_B(s,:);
    
    new_positionA(s,:)=m_axisA(s,:)+x1_dre(s)*new_x_vectorA(s,:);
    new_positionB(s,:)=m_axisB(s,:)+x2_dre(s)*new_x_vectorB(s,:);
    
    
    new_x1(s)=new_positionA(s,1);
    new_y1(s)=y1_dre(s);
    new_z1(s)=new_positionA(s,3);
    
    new_x2(s)=new_positionB(s,1);
    new_y2(s)=y2_dre(s);
    new_z2(s)=new_positionB(s,3);
    
    end 
elseif delta_alpha<0  %<0 means bending in right side
    r_bend=(z1(20)-z1(1))/delta_alpha;
    middle_punkt_A=[0 0 (z1(20)+z1(1))/2];
    middle_punkt_B=[0 0 (z2(20)+z2(1))/2];
    null_punkt_A=[r_bend 0 (z1(20)+z1(1))/2];

    null_punkt_B=[r_bend 0 (z2(20)+z2(1))/2];
   
    for s=1:20
    Beta(s)=(z1(s)-z1(1))/r_bend;
    
    m_axisA_x(s)=(r_bend-r_bend*cos(Beta(s)-delta_alpha/2));
    m_axisA_y(s)=0;
    m_axisA_z(s)=r_bend*sin(Beta(s)-delta_alpha/2)+((z1(20)-z1(1))/2+z1(1));
    m_axisA(s,:)=[m_axisA_x(s) m_axisA_y(s) m_axisA_z(s)];
    
    m_axisB_x(s)=(r_bend-r_bend*cos(Beta(s)-delta_alpha/2));
    m_axisB_y(s)=0;
    m_axisB_z(s)=r_bend*sin(Beta(s)-delta_alpha/2)+((z2(20)-z2(1))/2+z2(1));
    m_axisB(s,:)=[m_axisB_x(s) m_axisB_y(s) m_axisB_z(s)];
    
    
    
    new_x_vectorA(s,:)=m_axisA(s,:)-null_punkt_A;
    new_x_vectorA(s,:)=new_x_vectorA(s,:)/norm(new_x_vectorA(s,:));
    new_x_vectorB(s,:)=m_axisB(s,:)-null_punkt_B;
    new_x_vectorB(s,:)=new_x_vectorB(s,:)/norm(new_x_vectorB(s,:));
    
    m_proj_A(s,:)=(dot(z_axis,new_x_vectorA(s,:))/norm(new_x_vectorA(s,:))^2)*new_x_vectorA(s,:);
    new_zaxisA(s,:)=z_axis-m_proj_A(s,:);

    m_proj_B(s,:)=(dot(z_axis,new_x_vectorB(s,:))/norm(new_x_vectorB(s,:))^2)*new_x_vectorB(s,:);
    new_zaxisB(s,:)=z_axis-m_proj_B(s,:);
    
    new_positionA(s,:)=m_axisA(s,:)+x1_dre(s)*new_x_vectorA(s,:);
    new_positionB(s,:)=m_axisB(s,:)+x2_dre(s)*new_x_vectorB(s,:);
    
    
    new_x1(s)=new_positionA(s,1);
    new_y1(s)=y1_dre(s);
    new_z1(s)=new_positionA(s,3);
    
    new_x2(s)=new_positionB(s,1);
    new_y2(s)=y2_dre(s);
    new_z2(s)=new_positionB(s,3);
    end
else
new_x1=x1_dre';
new_y1=y1_dre';
new_z1=z1_dre';
new_x2=x2_dre';
new_y2=y2_dre';
new_z2=z2_dre';
for i=1:20
new_zaxisA(i,:)=z_axis;
new_zaxisB(i,:)=z_axis;
end
end


%%%%new coordinate calculation
%new C1a/C1b position
new_C1a=[new_x1;new_y1;new_z1]';
new_C1b=[new_x2;new_y2;new_z2]';

%vector calculation 

x_axis1_Spin1=(new_C1b(3,:)-new_C1a(3,:))/norm(new_C1b(3,:)-new_C1a(3,:)); 
y_axis1_Spin1=cross(x_axis1_Spin1,new_zaxisA(3,:));
y_axis1_Spin1=y_axis1_Spin1/norm(y_axis1_Spin1);

z_proj=(dot(new_zaxisA(3,:),x_axis1_Spin1)/norm(x_axis1_Spin1)^2)*x_axis1_Spin1;
z_axis_Spin1=z_axis-z_proj;
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis_Spin1=cross(x_axis1_Spin1,z_axis);
y_axis_Spin1=y_axis_Spin1/norm(y_axis_Spin1);

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)+y_axis1_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+new_zaxisA(3,:).*sin(beta_1);
%%find Spin1
Spin_1=new_C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%C Label position rechnen
Z_Spin1=cross(x_axis3_Spin1,x_axis1_Spin1);
Z_Spin1=Z_Spin1/norm(Z_Spin1);
y_Spin1=cross(x_axis1_Spin1,Z_Spin1);
y_Spin1=y_Spin1/norm(y_Spin1);
N1_1=C1a(3,:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin1+1.5*sin(54.2/360*2*pi)*y_Spin1;
C2_1=C1a(3,:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin1+2.5*sin(24.3/360*2*pi)*y_Spin1;
C4_1=C1a(3,:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin1+4.2*sin(52.6/360*2*pi)*y_Spin1;
C5_1=C1a(3,:)+3.8*cos(72.2/360*2*pi)*x_axis1_Spin1+3.8*sin(72.2/360*2*pi)*y_Spin1;
% C45_1=(C4_1+C5_1)/2;
N1C2_1=(N1_1+C2_1)/2;
X_Spin1=C5_1-N1_1;
X_Spin1=X_Spin1/norm(X_Spin1);
Y_Spin1=cross(Z_Spin1,X_Spin1);
Y_Spin1=Y_Spin1/norm(Y_Spin1);

% % %%rotation around molecular axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_1=normrnd(0,sigma_y); 
[X1, Y1 Z1] = AxisAngleRotate(X_Spin1,Y_Spin1,Z_Spin1,Y_Spin1,rotate_theta1_1);
%second rotation around z-axis with 5 grad
rotate_theta2_1=normrnd(0,sigma_z); 
[X2, Y2 Z2] = AxisAngleRotate(X1,Y1,Z1,Z1,rotate_theta2_1);
%find electron center
M=N1C2_1+X2*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction
x_axis1_Spin2=(new_C1a(7+(nd-1),:)-new_C1b(7+(nd-1),:))/norm(new_C1a(7+(nd-1),:)-new_C1b(7+(nd-1),:));
y_axis_Spin2=cross(x_axis1_Spin2,new_zaxisB(7+(nd-1),:));
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);


%%first rotation
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))-y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))-new_zaxisB(7+(nd-1),:).*sin(beta_2(nd));
%%find Spin2
Spin_2=new_C1b(7+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

%C Label position rechnen
Z_Spin2=cross(x_axis3_Spin2,x_axis1_Spin2);
Z_Spin2=Z_Spin2/norm(Z_Spin2);
y_Spin2=cross(x_axis1_Spin2,Z_Spin2);
y_Spin2=y_Spin2/norm(y_Spin2);
N1_2=C1b(7+(nd-1),:)+1.5*cos(54.2/360*2*pi)*x_axis1_Spin2+1.5*sin(54.2/360*2*pi)*y_Spin2;
C2_2=C1b(7+(nd-1),:)+2.5*cos(24.3/360*2*pi)*x_axis1_Spin2+2.5*sin(24.3/360*2*pi)*y_Spin2;
C4_2=C1b(7+(nd-1),:)+4.2*cos(52.6/360*2*pi)*x_axis1_Spin2+4.2*sin(52.6/360*2*pi)*y_Spin2;
C5_2=C1b(7+(nd-1),:)+3.8*cos(72.9/360*2*pi)*x_axis1_Spin2+3.8*sin(72.9/360*2*pi)*y_Spin2;
C45_2=(C4_2+C5_2)/2;
N1C2_2=(N1_2+C2_2)/2;
X_Spin2=C5_2-N1_2;
X_Spin2=X_Spin2/norm(X_Spin2);
Y_Spin2=cross(Z_Spin2,X_Spin2);
Y_Spin2=Y_Spin2/norm(Y_Spin2);

% % % %%rotation around molecular axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_2=normrnd(0,sigma_y); 
[X1_2, Y1_2 Z1_2] = AxisAngleRotate(X_Spin2,Y_Spin2,Z_Spin2,Y_Spin2,rotate_theta1_2);
%second rotaion around z-axis with 5 grad
rotate_theta2_2=normrnd(0,sigma_z); 
[X2_2, Y2_2 Z2_2] = AxisAngleRotate(X1_2,Y1_2,Z1_2,Z1_2,rotate_theta2_2);
%find electron center
M2=N1C2_2+X2_2*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);


end
R((w-1)*100+1:w*100,:)=M2-M;

[o1((w-1)*100+1:w*100,:),o2((w-1)*100+1:w*100,:),r((w-1)*100+1:w*100,:)]=RuCoordv1(R((w-1)*100+1:w*100,:),X2,Y2,Z2,X2_2,Y2_2,Z2_2);

end 
Conformers.M=M;
Conformers.M2=M2;
% Conformers.R=R;
Conformers.EulerAngles.R1=[o1(:,:)];
% Conformers.EulerAngles.R1(:,2)=Conformers.EulerAngles.R1(:,2)+ones(length(Conformers.EulerAngles.R1(:,2)),1)*deg2rad(10);
Conformers.EulerAngles.R2=[o2(:,:)];
% Conformers.EulerAngles.R2(:,2)=Conformers.EulerAngles.R2(:,2)+ones(length(Conformers.EulerAngles.R2(:,2)),1)*deg2rad(10);
Conformers.Distance = r/10;
zeiten = zeit*1000;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the program
pd=fitdist(Conformers.Distance,'Normal');
R_mean=pd.mu;
sigma=pd.sigma;
FWHM=2.3548*sigma;
h=histfit(Conformers.Distance);
h(1).FaceColor = [.3 .75 .93];
h(2).Color = [.0 .0 1];
end 


