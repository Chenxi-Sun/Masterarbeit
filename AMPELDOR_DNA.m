function [Result,R_mean,FWHM,ymax] = AMPELDOR_DNA(sigma_y,nd,EP,zeit,str)
% function [Result,R_mean,FWHM] = AMPELDOR_DNA(sigma_y,nd,EP,zeit)
% function [Result] = AMPELDOR_DNA(nd,sigma_r,EP,zeit)
% nd=5;
% for nd=5:14
% str='B';
sigma_z=5;
% sigma_y=6;
nd=nd-4;

n_bp=1:1:20; %20base pair
%%Parameter extracted from pymol and fitted by cftool
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
%%Parameter extracted from Nicole's skript
% r0=5.85;
% h0=33;
% b=0.6169;
% c1_x=-0.554;
% c1_y=1.0168;
% c2_x=-4.477;
% c2_y=-2.9062;
% d=3.2401;
% e1=1.18;
% e2=2.19;

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

%%rotationsangle 1st SL
alpha_1=76.20/360*2*pi;
beta_1=6.65/360*2*pi;
% beta_1=10/360*2*pi;
%%rotationsangle 2nd SL 5-14 bp;
alpha_2=[76.1446 76.1787 76.1766 76.2030 76.1871 76.1902 76.1983 76.1741 76.1634 76.1681]./360.*2.*pi;
beta_2=[6.6477 6.6232 6.6401 6.6957 6.7435 6.7840 6.7882 6.7712 6.7300 6.6836]./360.*2.*pi;


%%iteration begin
for k=1:500
switch (str)
    case 'B'
% parameter
deltar=normrnd(0,0.65);   %standard deviation of radius=0.65 
% deltar=0;
r=r0+deltar;
deltaL=-deltar*20/3.2;  %AM paper
L=L0+deltaL; 

new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r)^2+h0^2); %countor length=n_turns*sqrt((2pir)^2+h^2)
new_b=new_turns*2*pi/19; %the range of b.*n_bp+c1 = how much 2pi = how much turns
new_c1_x=b+c1_x-new_b;  %keep the first position for C1a(1) (x1,y1,z1) same
new_c1_y=b+c1_y-new_b;
new_d=new_turns*h0/19;  
new_e1=d+e1-new_d;
new_c2_x=b+c2_x-new_b;
new_c2_y=b+c2_y-new_b;
new_e2=20*d-20*new_d+e2+deltaL; %the last position for C1b(20) is (x2,y2,z2+deltaL)

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
% deltah=0;
h=h0+deltah;
deltaL=deltah*20/18.2;  %AM paper
L=L0+deltaL; 
r=r0;
new_turns=n_turns.*sqrt((2*pi*r0)^2+h0^2)/sqrt((2*pi*r0)^2+h^2); %countor length=n_turns*sqrt((2pir)^2+h^2)
new_b=new_turns*2*pi/19; %the range of b.*n_bp+c1 = how much 2pi = how much turns; -0.0138: when deltar=0,keep new_b = b
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

% % %%rotation around molecular axis 
%first rotaion around y-axis
rotate_theta1_1(k,:)=normrnd(0,sigma_y); 
[X1(k,:), Y1(k,:) Z1(k,:)] = AxisAngleRotate(X_Spin1(k,:),Y_Spin1(k,:),Z_Spin1(k,:),Y_Spin1(k,:),rotate_theta1_1(k,:));
%second rotation around z-axis
rotate_theta2_1(k,:)=normrnd(0,sigma_z); 
[X2(k,:), Y2(k,:) Z2(k,:)] = AxisAngleRotate(X1(k,:),Y1(k,:),Z1(k,:),Z1(k,:),rotate_theta2_1(k,:));
%find electron center
M(k,:)=N1C2_1(k,:)+X2(k,:)*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);

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

% % % %%rotation around molecular axis 
%first rotaion around y-axis with 6 grad
rotate_theta1_2(k,:)=normrnd(0,sigma_y); 
[X1_2(k,:), Y1_2(k,:) Z1_2(k,:)] = AxisAngleRotate(X_Spin2(k,:),Y_Spin2(k,:),Z_Spin2(k,:),Y_Spin2(k,:),rotate_theta1_2(k,:));
%second rotaion around z-axis with 5 grad
rotate_theta2_2(k,:)=normrnd(0,sigma_z); 
[X2_2(k,:), Y2_2(k,:) Z2_2(k,:)] = AxisAngleRotate(X1_2(k,:),Y1_2(k,:),Z1_2(k,:),Z1_2(k,:),rotate_theta2_2(k,:));
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
zeiten = zeit*1000;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the programr
% Result = MainPELDORtime_modAC(EP,Conformers,zeiten,6.5); %for G-band
h=histfit(Conformers.Distance);   %plot distribution
h(1).FaceColor = [.3 .75 .93];
h(2).Color = [.0 .0 1];
y=get(h,'YData');  
y1=cell2mat(y(1));
ymax=max(y1);     %%max value of histogramm
% Distance(nd,:)=r/10;
end 
% 
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')
% plot(5:14,Distance(:,1),'Linewidth',3,'Marker','o')
% hold on 
% plot(5:14,expDistances.DNA.meandist(1:10),'Linewidth',3,'Marker','o')

% for s=1:10
% bpd=s;
% pd=fitdist(Distance(bpd,:)','Normal');
% mu(s)=pd.mu;
% sigma(s)=pd.sigma;
% end 
% x=1:0.01:5.5;
% for p=1:s
% y(p,:)=normpdf(x,mu(p),sigma(p));
% plot(x,y(p,:),'LineWidth',3)
% % xlim([0 5])
% hold on 
% end 
% xlim([1 5.5])
% set(gca,'FontSize',14,'FontWeight','bold','XTick',...
%     [1 2 3 4 5 6]);
% set(gca,'linewidth',1.5) 
% box off
% xlabel('Distance [nm]')
% ylabel('Normalised probability')