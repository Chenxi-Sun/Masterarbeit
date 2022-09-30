nd=6;
nd=nd-4;
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
for w=1:361
theta=(w-1);
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
for k=1:50
alpha=normrnd(0,28)/360*2*pi;
if alpha<0
    alpha=-alpha;
elseif alpha==0
    alpha=normrnd(0,28)/360*2*pi;
else
    alpha=alpha;
% alpha=28/360*2*pi;
end 

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

plot3(new_x1,new_y1,new_z1,'Marker','x','LineWidth',2)
plot3(new_x2,new_y2,new_z2,'Marker','x','LineWidth',2)
hold on 
pause(.1)
end 
clf
end 
