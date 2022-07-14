%% 
% function [Result,R_mean] = AM_PELDOR_ApriRNA(sigma_r,nd,EP,zeit)
% function [Result,R_mean] = AM_PELDOR_ApriRNA(sigma_h,nd,EP,zeit)
function [Result,R_mean,FWHM] = AM_PELDOR_ApriRNA(sigma_y,nd,EP,zeit,str)
% function [Result,R_mean,FWHM,ymax] = AM_PELDOR_ApriRNA(sigma_y,nd,EP,zeit)
% % function [Result,R_mean] = AM_PELDOR_ApriRNA(nd,EP,zeit)
% str='B'
% for nd=8:15
% sigma_y=0;
sigma_z=5;
%%Parameter extracted from pymol and fitted by cftool
n_bp=1:1:20; %20base pair
% nd=8;
% sigma_r=1;
nd=nd-7;
%%initial position for C from pymol

%parameter for expression for C1-points at RNA
r0=8.837;
h0=35.5686;
b=0.522;
c1_x=3.565;
c1_y=1.983;
c2_x=-1.468;
c2_y=3.273;
d=2.955;
e1=-2.177;
e2=-3.728;

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


%%rotationsangle 1st SL
alpha_1=74.42/360*2*pi;
% beta_1=5.2/360*2*pi; 
beta_1=3.6/360*2*pi; 

%%rotationsangle 2nd SL 8-15 bp; (5-7)samples dont exist
% alpha_2=[74.4187 74.6188 76.6019 78.0654 77.2759 75.1897 74.2406 75.5165]./360.*2.*pi;
% beta_2=[5.2154 5.0264 4.3064 3.6792 3.6413 4.0283 4.3597 4.3319]./360.*2.*pi;



p=500;
for k=1:p
switch (str)
    case 'B'
% parameter
% deltar=normrnd(0,sigma_r);   %standard deviation of radius=0.65 
deltar=normrnd(0,0.9);   %standard deviation of radius=0.65 
% deltar=0;
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


% %%%%% atan2(norm(cross(vector1,vector2)), dot(vector1,vector2))/2/pi*360 %angle cal
% %%%%% vector2_proj=(dot(vector2,z_axis)/norm(z_axis))*z_axis; %projection
% 
% %



%%
   case 'A'
% Model A
% % % %parameter
deltah=normrnd(0,5);   %standard deviation of pitch height= 3.85 
% deltah=normrnd(0,sigma_h);
% deltah=0;
h=h0+deltah;
deltaL=deltah*20/20.9712;  %AM paper
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
[X1(k,:), Y1(k,:) Z1(k,:)] = AxisAngleRotate(X_Spin1(k,:),Y_Spin1(k,:),Z_Spin1(k,:),Y_Spin1(k,:),rotate_theta1_1(k,:));
%second rotation around z-axis with 5 grad
rotate_theta2_1(k,:)=normrnd(0,sigma_z); 
[X2(k,:), Y2(k,:) Z2(k,:)] = AxisAngleRotate(X1(k,:),Y1(k,:),Z1(k,:),Z1(k,:),rotate_theta2_1(k,:));
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
x_axis2_Spin2=x_axis1_Spin2.*cos(alpha_1)+y_axis_Spin2.*sin(alpha_1);
%%second rotation
x_axis3_Spin2=x_axis2_Spin2.*cos(beta_1)-z_axis_Spin2.*sin(beta_1);
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
[X1_2(k,:), Y1_2(k,:) Z1_2(k,:)] = AxisAngleRotate(X_Spin2(k,:),Y_Spin2(k,:),Z_Spin2(k,:),Y_Spin2(k,:),rotate_theta1_2(k,:));
%second rotaion around z-axis with 5 grad
rotate_theta2_2(k,:)=normrnd(0,sigma_z); %range -5 to 5grad
[X2_2(k,:), Y2_2(k,:) Z2_2(k,:)] = AxisAngleRotate(X1_2(k,:),Y1_2(k,:),Z1_2(k,:),Z1_2(k,:),rotate_theta2_2(k,:));
%find electron center
M2(k,:)=N1C2_2(k,:)+X2_2(k,:)*((sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2+2.3);


end 
R=M2-M;
% [o1,o2,r]=RuCoordv1(R,X_Spin1,Y_Spin1,Z_Spin1,X_Spin2,Y_Spin2,Z_Spin2);
[o1,o2,r]=RuCoordv1(R,X2,Y2,Z2,X2_2,Y2_2,Z2_2);

Conformers.M=M;
Conformers.M2=M2;
% Conformers.R=R;
Conformers.EulerAngles.R1=[o1(:,:)];
% Conformers.EulerAngles.R1(:,2)=Conformers.EulerAngles.R1(:,2)+ones(length(Conformers.EulerAngles.R1(:,2)),1)*deg2rad(10);
Conformers.EulerAngles.R2=[o2(:,:)];
% Conformers.EulerAngles.R2(:,2)=Conformers.EulerAngles.R2(:,2)+ones(length(Conformers.EulerAngles.R2(:,2)),1)*deg2rad(10);
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\r_y_6_grad_rot\delta_Distance.mat')
% Conformers.Distance = r/10+delta_Distance(nd);
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\delta_Distance_sphere.mat')
% Conformers.Distance = r/10+delta_Distance_sphere(nd);
Conformers.Distance = r/10;
% R_mean=sum(Conformers.Distance)/length(Conformers.Distance);
pd=fitdist(Conformers.Distance,'Normal');
R_mean=pd.mu;
sigma=pd.sigma;
FWHM=2.3548*sigma;
zeiten = zeit*1000;
Result = MainPELDORtime(EP,Conformers,zeiten); %...time lets you set the time axis from outside the program
% Result = MainPELDORtime_modAC(EP,Conformers,zeiten,6.5); %for G-band
% mean_trend=mean(trend);
% h=histfit(Conformers.Distance);
% h(1).FaceColor = [1 0.411764705882353 0.16078431372549];
% h(2).Color = [1 0 0];
% y=get(h,'YData');
% y1=cell2mat(y(1));
% ymax=max(y1);
% Distance(nd,:)=r/10;
end

% for s=1:8
% bpd=s;
% pd=fitdist(Distance(bpd,:)','Normal');
% mu(s)=pd.mu;
% sigma(s)=pd.sigma;
% end 
% x=1:0.01:5.5;
% for p=1:s
% y(p,:)=normpdf(x,mu(p),sigma(p));
% plot(x,y(p,:),'LineWidth',2)
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