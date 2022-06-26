%% 

clear;
str='B';
%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair

%%initial position for C from pymol

%parameter for expression for C1-points at RNA
r0=8.284;
h0=29.8320;
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

%end-end length
L0=z2(20)-z1(1);

%contour length C
C=sqrt((2*pi*r0)^2*L0^2/h0^2+L0^2); %countor length (AM paper)

%n_turns
n_turns=b.*19/2/pi;  %the range of b.*n_bp+c1 = how much 2pi = how much turns

%z_axis
% Z=[x1(11)-x1(1) y1(11)-y1(1) z1(11)-z1(1)];
% z_axis=Z/norm(Z);
z_axis=[0 0 1];


%%rotationsangle 1st SL
alpha_1=74.42/360*2*pi;
beta_1=5.2/360*2*pi; 
%%rotationsangle 2nd SL 8-15 bp; (5-7)samples dont exist
alpha_2=[74.4187 74.6188 76.6019 78.0654 77.2759 75.1897 74.2406 75.5165]./360.*2.*pi;
beta_2=[5.2154 5.0264 4.3064 3.6792 3.6413 4.0283 4.3597 4.3319]./360.*2.*pi;


p=500;
for k=1:p
switch (str)
    case 'B'
% parameter
deltar=normrnd(0,0.65);   %standard deviation of radius=0.65 
% deltar=0;
r=r0+deltar;
deltaL=-deltar*20/3.2;  %AM paper
L=L0+deltaL; 

new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r)^2+h0^2); %countor length=n_turns*sqrt((2pir)^2+h^2)

new_b=new_turns*2*pi/19; %the range of b.*n_bp+c1 = how much 2pi = how much turns<<<<<for r=5.85A
new_c1_x=b+c1_x-new_b;  %keep the first position for C1a(1) (x1,y1,z1) same
new_c1_y=b+c1_y-new_b;
new_d=new_turns*h0/19+0.0014;  
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
% % % %parameter
deltah=normrnd(0,3.85);   %standard deviation of pitch height= 3.85 
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
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha_1)-y_axis_Spin1.*sin(alpha_1);
%%second rotation
x_axis3_Spin1=x_axis2_Spin1.*cos(beta_1)-z_axis_Spin1.*sin(beta_1);
%%find Spin1
Spin_1=C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A
% 
% theta1=rand(1,1)*pi;
% phi1=rand(1,1)*2*pi;
% M=Spin_1+[sin(theta1)*cos(phi1) sin(theta1)*sin(phi1) cos(theta1)];
% x_axis4_Spin1=M-C1a(3,:);
% x_axis4_Spin1=x_axis4_Spin1/norm(x_axis4_Spin1);
% M=C1a(3,:)+x_axis4_Spin1*11.4;


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction



for i=1:8
x_axis1_Spin2(i,:)=(C1a(10+(i-1),:)-C1b(10+(i-1),:))/norm(C1a(10+(i-1),:)-C1b(10+(i-1),:));
z_proj_Spin2(i,:)=(dot(z_axis,x_axis1_Spin2(i,:))/norm(x_axis1_Spin2(i,:))^2)*x_axis1_Spin2(i,:);
z_axis_Spin2(i,:)=z_axis-z_proj_Spin2(i,:);
z_axis_Spin2(i,:)=z_axis_Spin2(i,:)/norm(z_axis_Spin2(i,:));
y_axis_Spin2(i,:)=cross(x_axis1_Spin2(i,:),z_axis);
y_axis_Spin2(i,:)=y_axis_Spin2(i,:)/norm(y_axis_Spin2(i,:));

%%first rotation
x_axis2_Spin2(i,:)=x_axis1_Spin2(i,:).*cos(alpha_2(i))+y_axis_Spin2(i,:).*sin(alpha_2(i));
%%second rotation
x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta_1)+z_axis_Spin2(i,:).*sin(beta_1);
%%find Spin2
Spin_2(i,:)=C1b(10+(i-1),:)+x_axis3_Spin2(i,:).*11.4;  %length of C1 and spinlabel is 11.5A
% theta2(i)=rand(1,1)*pi;
% phi2(i)=rand(1,1)*2*pi;
% M2(i,:)=Spin_2(i,:)+[sin(theta2(i))*cos(phi2(i)) sin(theta2(i))*sin(phi2(i)) cos(theta2(i))];
% x_axis4_Spin2(i,:)=M2(i,:)-C1b(10+(i-1),:);
% x_axis4_Spin2(i,:)=x_axis4_Spin2(i,:)/norm(x_axis4_Spin2(i,:));
% M2(i,:)=C1b(10+(i-1),:)+x_axis4_Spin2(i,:).*11.4;
end

R(:,k)=sqrt(sum((Spin_2-Spin_1).^2,2));
% R(:,k)=sqrt(sum((M2-M).^2,2));
end 
mean_R=sum(R,2)./500;
% load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')
load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')
figure(1)
plot(8:15,expDistances.RNA.meandist,'Linewidth',2,'Marker','o')
hold on 
plot(8:15,mean_R/10,'Linewidth',2,'Marker','o')
ylim([1.5 4.5])
xlabel('position of 2^{nd} spin label')
ylabel('distance [nm]')
legend('PELDOR','simulated')