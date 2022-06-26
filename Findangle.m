clear;
str='B';
%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

%%initial position for C from pymol

% %parameter for expression for C1-points at DNA from pymol
% r0=5.68;
% h0=33.75; %h0 for C
% b=0.6288;
% c1_x=-1.691;
% c1_y=-0.1296;
% c2_x=0.7453;
% c2_y=2.324;
% d=3.375;
% e1=-69.2;
% e2=-69.17;
load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')
%parameter for expression for C1-points at DNA for AM
r0=5.85;
h0=33;
b=0.6253;
c1_x=-1.6875;
c1_y=-0.1261;
c2_x=0.7488;
c2_y=2.3275;
d=3.2814;
e1=-69.1064;
e2=-69.1862;

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

%%range for both angle 
alpha_range=76.2;
beta_range=[0.1:0.01:10];

p=991;
for k=1:p
switch (str)
    case 'B'
% parameter
% deltar=normrnd(0,0.65);   %standard deviation of radius=0.65 
deltar=0;
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


% %%%%% atan2(norm(cross(vector1,vector2)), dot(vector1,vector2))/2/pi*360 %angle cal
% %%%%% vector2_proj=(dot(vector2,z_axis)/norm(z_axis))*z_axis; %projection
% 
% %



%%
   case 'A'
% Model A
% % %parameter
% deltah=normrnd(0,3.85);   %standard deviation of pitch height= 3.85 
deltah=0;
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


alpha=alpha_range/360*2*pi;               %first rotation angle alpha
beta(k)=beta_range(k)/360*2*pi;                %second rotation angle beta

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha)+y_axis_Spin1.*sin(alpha);
%%second rotation
% x_axis3_Spin1=x_axis2_Spin1.*cos(beta)+z_axis.*sin(beta);
x_axis3_Spin1=x_axis2_Spin1.*cos(beta(k))+z_axis_Spin1.*sin(beta(k));
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A
% 
% theta1=rand(1,1)*2*pi;
% phi1=rand(1,1)*2*pi;
% Spin_1=Spin_1+[sin(theta1)*cos(phi1) sin(theta1)*sin(phi1) cos(theta1)];
%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

for i=1:10
x_axis1_Spin2(i,:)=(C1a(7+(i-1),:)-C1b(7+(i-1),:))/norm(C1a(7+(i-1),:)-C1b(7+(i-1),:));
z_proj_Spin2(i,:)=(dot(z_axis,x_axis1_Spin2(i,:))/norm(x_axis1_Spin2(i,:))^2)*x_axis1_Spin2(i,:);
z_axis_Spin2(i,:)=z_axis-z_proj_Spin2(i,:);
z_axis_Spin2(i,:)=z_axis_Spin2(i,:)/norm(z_axis_Spin2(i,:));
y_axis_Spin2(i,:)=cross(x_axis1_Spin2(i,:),z_axis);
y_axis_Spin2(i,:)=y_axis_Spin2(i,:)/norm(y_axis_Spin2(i,:));

%%first rotation
x_axis2_Spin2(i,:)=x_axis1_Spin2(i,:).*cos(alpha)-y_axis_Spin2(i,:).*sin(alpha);
%%second rotation
% x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta)-z_axis.*sin(beta);
x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta(k))-z_axis_Spin2(i,:).*sin(beta(k));
%%find Spin2
Spin_2(i,:)=C1b(7+(i-1),:)+x_axis3_Spin2(i,:).*11.4;  %length of C1 and spinlabel is 11.5A
% theta2(i)=rand(1,1)*2*pi;
% phi2(i)=rand(1,1)*2*pi;
% Spin_2(i,:)=Spin_2(i,:)+[sin(theta2(i))*cos(phi2(i)) sin(theta2(i))*sin(phi2(i)) cos(theta2(i))];
end


R(:,k)=sqrt(sum((Spin_2-Spin_1).^2,2))./10;

% R0(k)=sqrt(sum((R(:,k)'-expDistances.DNA.meandist(1:10).*10).^2,2));
R0(k)=R(5,k)'-expDistances.DNA.meandist(5);

end 


% rmsd=sqrt(sum((R-expDistances.DNA.meandist.*10).^2,2)./p);
[m,n]=min(abs(R0));
beta_1=beta(n)/2/pi*360
R(5,n)
R0(n)