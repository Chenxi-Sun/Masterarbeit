clear;
load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')

n_bp=1:20;

% %parameter for expression for C1-points at ApriRNA
% r0=8.837;
% h0=35.5686;
% b=0.522;
% c1_x=3.565;
% c1_y=1.983;
% c2_x=-1.468;
% c2_y=3.273;
% d=2.955;
% e1=-2.177;
% e2=-3.728;

%parameter for expression for C1-points at A-RNA
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
n_turns=b.*19/2/pi;  %the range of b.*n_bp+c1 = how much 2pi = how much turns

%z_axis
z_axis=[0 0 1];

%
deltar=0;
r=r0+deltar;
deltaL=-deltar*20/3.8;  %AM paper
L=L0+deltaL; 

new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r)^2+h0^2); %countor length=n_turns*sqrt((2pir)^2+h^2)

new_b=new_turns*2*pi/19; %the range of b.*n_bp+c1 = how much 2pi = how much turns<<<<<for r=5.85A
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

%%rotationsangle 1st SL
alpha_1=74.42/360*2*pi;
beta_1=5.2/360*2*pi; 
%%rotationsangle 2nd SL 8-15 bp; (5-7)samples dont exist
alpha_2=[74.4187 74.6188 76.6019 78.0654 77.2759 75.1897 74.2406 75.5165]./360.*2.*pi;
beta_2=[5.2154 5.0264 4.3064 3.6792 3.6413 4.0283 4.3597 4.3319]./360.*2.*pi;

for i=1:30
ratio=1+(i-1)*0.01;
new_z1_1=new_z1*ratio;
new_z2_1=new_z2*ratio;
%new C1a/C1b position
C1a=[new_x1;new_y1;new_z1_1]';
C1b=[new_x2;new_y2;new_z2_1]';

for nd=8:15
nd=nd-7;
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
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

x_axis1_Spin2=(C1a(10+(nd-1),:)-C1b(10+(nd-1),:))/norm(C1a(10+(nd-1),:)-C1b(10+(nd-1),:));
z_proj_Spin2=(dot(z_axis,x_axis1_Spin2)/norm(x_axis1_Spin2)^2)*x_axis1_Spin2;
z_axis_Spin2=z_axis-z_proj_Spin2;
z_axis_Spin2=z_axis_Spin2/norm(z_axis_Spin2);
y_axis_Spin2=cross(x_axis1_Spin2,z_axis_Spin2);
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);

%%first rotation
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))+y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))-z_axis_Spin2.*sin(beta_2(nd));
%%find Spin2
Spin_2=C1b(10+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

Distance(nd,:)=norm(Spin_1-Spin_2)/10;

end 

rmsd(i)=sum((expDistances.RNA.meandist-Distance').^2,2)/8;

end 
[m,n]=min(rmsd);
best_ratio=1+(n-1)*0.01;
% best_ratio=1
new_z1_2=new_z1*best_ratio;
new_z2_2=new_z2*best_ratio;
%new C1a/C1b position
C1a=[new_x1;new_y1;new_z1_2]';
C1b=[new_x2;new_y2;new_z2_2]';

for nd=8:15
nd=nd-7;
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
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)+z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A

%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

x_axis1_Spin2=(C1a(10+(nd-1),:)-C1b(10+(nd-1),:))/norm(C1a(10+(nd-1),:)-C1b(10+(nd-1),:));
z_proj_Spin2=(dot(z_axis,x_axis1_Spin2)/norm(x_axis1_Spin2)^2)*x_axis1_Spin2;
z_axis_Spin2=z_axis-z_proj_Spin2;
z_axis_Spin2=z_axis_Spin2/norm(z_axis_Spin2);
y_axis_Spin2=cross(x_axis1_Spin2,z_axis_Spin2);
y_axis_Spin2=y_axis_Spin2/norm(y_axis_Spin2);

%%first rotation
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_2(nd))+y_axis_Spin2.*sin(alpha_2(nd));
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_2(nd))-z_axis_Spin2.*sin(beta_2(nd));
%%find Spin2
Spin_2=C1b(10+(nd-1),:)+x_axis3_Spin2.*11.4;  %length of C1 and spinlabel is 11.5A

Distance(nd,:)=norm(Spin_1-Spin_2)/10;

end 

plot(8:15,Distance,'Linewidth',3,'Marker','o');
hold on 
plot(8:15,expDistances.RNA.meandist,'Linewidth',3,'Marker','o')

