n_bp=1:1:20; %20base pair

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



load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\distance_compared_delta=0.mat')
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

%drehen die koordinatensystem
Z_axis=[0 0 1];
X_axis=[1 0 0];
Y_axis=[0 1 0];

%bending winkel
delta_alpha=28/360*2*pi;

p=37;
for o=1:p
    
theta=(o-1)/(p-1)*2*pi;
%%C1 points nach der Drehnung
[new_X new_Y new_Z]=drehnung(X_axis,Y_axis,Z_axis,theta);
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

%%
m_axisA_x=zeros(1,20);
m_axisA_y=zeros(1,20);
m_axisA_z=z1;

m_axisB_x=zeros(1,20);
m_axisB_y=zeros(1,20);
m_axisB_z=z2;

r_bend=(z1(20)-z1(1))/delta_alpha; 
middle_punkt_A=[0 0 (z1(20)+z1(1))/2];
middle_punkt_B=[0 0 (z2(20)+z2(1))/2];

null_punkt_A=[-r_bend 0 (z1(20)+z1(1))/2];
null_punkt_B=[-r_bend 0 (z2(20)+z2(1))/2];

for s=1:20
Beta(s)=(z1(s)-z1(1))/r_bend;

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

%%%%new coordinate calculation
%new C1a/C1b position
new_C1a=[new_x1;new_y1;new_z1]';
new_C1b=[new_x2;new_y2;new_z2]';

%vector calculation 

x_axis1_Spin1=(new_C1b(3,:)-new_C1a(3,:))/norm(new_C1b(3,:)-new_C1a(3,:)); 
y_axis1_Spin1=cross(x_axis1_Spin1,new_zaxisA(3,:));
y_axis1_Spin1=y_axis1_Spin1/norm(y_axis1_Spin1);


alpha=75.5/360*2*pi;               %first rotation angle alpha
beta=6.8/360*2*pi;               %second rotation angle beta

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha)+y_axis1_Spin1.*sin(alpha);
%%second rotation
% x_axis3_Spin1=x_axis2_Spin1.*cos(beta)+z_axis.*sin(beta);
x_axis3_Spin1=x_axis2_Spin1.*cos(beta)+new_zaxisA(3,:).*sin(beta);
%%find Spin1
Spin_1=new_C1a(3,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A
theta1=rand(1,1)*2*pi;
phi1=rand(1,1)*2*pi;
Spin_1=Spin_1+[sin(theta1)*cos(phi1) sin(theta1)*sin(phi1) cos(theta1)];
% axis4=[-1 + 2.*rand(1,1) -1 + 2.*rand(1,1) -1 + 2.*rand(1,1)]


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

for i=1:10
x_axis1_Spin2(i,:)=(new_C1a(7+(i-1),:)-new_C1b(7+(i-1),:))/norm(new_C1a(7+(i-1),:)-new_C1b(7+(i-1),:));
y_axis_Spin2(i,:)=cross(x_axis1_Spin2(i,:),new_zaxisB(7+(i-1),:));
y_axis_Spin2(i,:)=y_axis_Spin2(i,:)/norm(y_axis_Spin2(i,:));


%%first rotation
x_axis2_Spin2(i,:)=x_axis1_Spin2(i,:).*cos(alpha)-y_axis_Spin2(i,:).*sin(alpha);
%%second rotation
% x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta)-z_axis.*sin(beta);
x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta)-new_zaxisB(7+(i-1),:).*sin(beta);
%%find Spin2
Spin_2(i,:)=new_C1b(7+(i-1),:)+x_axis3_Spin2(i,:).*11.4;  %length of C1 and spinlabel is 11.5A
theta2(i)=rand(1,1)*2*pi;
phi2(i)=rand(1,1)*2*pi;
Spin_2(i,:)=Spin_2(i,:)+[sin(theta2(i))*cos(phi2(i)) sin(theta2(i))*sin(phi2(i)) cos(theta2(i))];
end


R(:,o)=sqrt(sum((Spin_2-Spin_1).^2,2));
end

rmsd=sqrt((R-Distance_compare.*10).^2);
Distance=R_avg./10;

for w=1:10
bpd=w;
pd=fitdist(Distance(bpd,:)','Normal');
mu(w,o)=pd.mu;
sigma(w,o)=pd.sigma;
end 



rmsd=sum(rmsd,2)/p;
mu=sum(mu,2)/p;
sigma=sum(sigma,2)/p;
% Distance=R./10;
figure(1)
plot(5:14,rmsd,'Marker','o','LineWidth',2)   
xlabel('position of 2^{nd} spin label')
ylabel('r.m.s.d [A]') 
ylim([0 3.5])

figure(2)
x=0:0.01:5;
for q=1:10
y(q,:)=normpdf(x,mu(q),sigma(q));
plot(x,y(q,:),'LineWidth',2)
xlim([0 5])
hold on 
end 
