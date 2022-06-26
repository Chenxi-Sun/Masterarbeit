% function [pd]=AMmodelC1_DNA

% prompt2='Which base pair difference? (4-13):';
% bpd=input(prompt2)-3;

%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair
r0=8.50; %r0 for C1 
h0=28.028; %h0 for C
%%initial position for C from pymol

%parameter for expression for C1-points at DNA
b=0.5707;
c1_x=3.352;
c1_y=1.781;
c2_x=2.094;
c2_y=0.5231;
d=2.548;
e1=-52.84;
e2=-49.08;

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
Z=[x1(12)-x1(1) y1(12)-y1(1) z1(12)-z1(1)];
z_axis=Z/norm(Z);

C1a(:,1)=x1;
C1a(:,2)=y1;
C1a(:,3)=z1;
C1b(:,1)=x2;
C1b(:,2)=y2;
C1b(:,3)=z2;

%drehen die koordinatensystem
Z_axis=[0 0 1];
X_axis=[1 0 0];
Y_axis=[0 1 0];

for o=1:37
    
theta=(o-1)*10/360*2*pi;
% theta=0/360*2*pi;

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


for i=1:20
C_vector(i,:)=C1a_dre(i,:)-C1b_dre(i,:);
C_vector(i,:)=C_vector(i,:)/norm(C_vector(i,:));
zaxis_proj(i,:)=(dot(z_axis,C_vector(i,:))/norm(C_vector(i,:))^2)*C_vector(i,:);
new_zaxis(i,:)=z_axis-zaxis_proj(i,:);
new_zaxis(i,:)=new_zaxis(i,:)/norm(new_zaxis(i,:));

end 

%middle axis cross the zylinder
for k=1:500
delta_alpha=randi([-28 28],1,1)/360*2*pi;
% delta_alpha=normrnd(0,28);
% delta_alpha=28/360*2*pi;

m_axisA_x=zeros(1,20);
m_axisA_y=zeros(1,20);
m_axisA_z=z1;

m_axisB_x=zeros(1,20);
m_axisB_y=zeros(1,20);
m_axisB_z=z2;

%%after bending 
if delta_alpha>0
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
    new_zaxisA(s,:)=new_zaxis(s,:)-m_proj_A(s,:);

    m_proj_B(s,:)=(dot(z_axis,new_x_vectorB(s,:))/norm(new_x_vectorB(s,:))^2)*new_x_vectorB(s,:);
    new_zaxisB(s,:)=new_zaxis(s,:)-m_proj_B(s,:);
    
    new_positionA(s,:)=m_axisA(s,:)+x1_dre(s)*new_x_vectorA(s,:);
    new_positionB(s,:)=m_axisB(s,:)+x2_dre(s)*new_x_vectorB(s,:);
    
    
    new_x1(s)=new_positionA(s,1);
    new_y1(s)=y1_dre(s);
    new_z1(s)=new_positionA(s,3);
    
    new_x2(s)=new_positionB(s,1);
    new_y2(s)=y2_dre(s);
    new_z2(s)=new_positionB(s,3);
    
    end 
elseif delta_alpha<0
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
    new_zaxisA(s,:)=new_zaxis(s,:)-m_proj_A(s,:);
    
    m_proj_B(s,:)=(dot(z_axis,new_x_vectorB(s,:))/norm(new_x_vectorB(s,:))^2)*new_x_vectorB(s,:);
    new_zaxisB(s,:)=new_zaxis(s,:)-m_proj_B(s,:);
    
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
new_x1=x1;
new_y1=y1;
new_z1=z1;
new_x2=x2;
new_y2=y2;
new_z2=z2;
end


%%%%new coordinate calculation
%new C1a/C1b position
new_C1a=[new_x1;new_y1;new_z1]';
new_C1b=[new_x2;new_y2;new_z2]';

%vector calculation 

x_axis1_Spin1=(new_C1b(18,:)-new_C1a(18,:))/norm(new_C1b(18,:)-new_C1a(18,:)); 
z_axis_Spin1=new_zaxisA(18,:);
z_axis_Spin1=z_axis_Spin1/norm(z_axis_Spin1);
y_axis1_Spin1=cross(z_axis_Spin1,x_axis1_Spin1);
y_axis1_Spin1=y_axis1_Spin1/norm(y_axis1_Spin1);


alpha=74.59/360*2*pi;               %first rotation angle alpha
beta=4.46/360*2*pi;               %second rotation angle beta

%%first rotation
x_axis2_Spin1=x_axis1_Spin1.*cos(alpha)-y_axis1_Spin1.*sin(alpha);
%%second rotation
% x_axis3_Spin1=x_axis2_Spin1.*cos(beta)+z_axis.*sin(beta);
x_axis3_Spin1=x_axis2_Spin1.*cos(beta)-z_axis_Spin1.*sin(beta);
%%find Spin1
Spin_1=new_C1a(18,:)+x_axis3_Spin1.*11.4;  %length of C1 and spinlabel is 11.4A
% axis4=[-1 + 2.*rand(1,1) -1 + 2.*rand(1,1) -1 + 2.*rand(1,1)]


%%Spin2 position 
%%rotation equation:x2=x*cos(theta)+y*sin(theta),y is direction

for i=1:10
x_axis1_Spin2(i,:)=(new_C1a(13-(i-1),:)-new_C1b(13-(i-1),:))/norm(new_C1a(13-(i-1),:)-new_C1b(13-(i-1),:));
z_axis_Spin2(i,:)=new_zaxisB(13-(i-1),:);
z_axis_Spin2(i,:)=z_axis_Spin2(i,:)/norm(z_axis_Spin2(i,:));
y_axis_Spin2(i,:)=cross(z_axis_Spin2(i,:),x_axis1_Spin2(i,:));
y_axis_Spin2(i,:)=y_axis_Spin2(i,:)/norm(y_axis_Spin2(i,:));


%%first rotation
x_axis2_Spin2(i,:)=x_axis1_Spin2(i,:).*cos(alpha)+y_axis_Spin2(i,:).*sin(alpha);
%%second rotation
x_axis3_Spin2(i,:)=x_axis2_Spin2(i,:).*cos(beta)+z_axis_Spin2(i,:).*sin(beta);
%%find Spin2
Spin_2(i,:)=new_C1b(13-(i-1),:)+x_axis3_Spin2(i,:).*11.4;  %length of C1 and spinlabel is 11.5A
end


R(:,k)=sqrt(sum((Spin_2-Spin_1).^2,2));
end
load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\distance_pymol.mat')
rmsd(:,o)=sqrt(sum((R-distance_pymol_RNA_5to14).^2,2)./500);

Distance=R./10;
for w=1:10
bpd=w;
pd=fitdist(Distance(bpd,:)','Normal');
mu(w,o)=pd.mu;
sigma(w,o)=pd.sigma;
end 
end 

rmsd=sum(rmsd,2)/37;

mu=sum(mu,2)/37;
sigma=sum(sigma,2)/37;

figure(1)
plot(5:14,rmsd,'Marker','o','LineWidth',2)   
xlabel('position of 2^{nd} spin label')
ylabel('r.m.s.d [A]') 
ylim([0 3.5])
% 


figure(2)
x=0:0.01:5;
for p=1:10
y(p,:)=normpdf(x,mu(p),sigma(p));
plot(x,y(p,:),'LineWidth',2)
xlim([0 5])
hold on 
end 
xlabel('distance [nm]')