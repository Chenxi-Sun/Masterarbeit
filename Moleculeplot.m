xyz=importdata('newClabelDFT - Kopie.xyz')

coords = xyz.data;
labels=xyz.textdata;

figure
set(gcf,'Color','w')
molecule3D(coords,labels);
cameratoolbar
hold on 
N = coords(6,:);
C5 = coords(3,:);
C2 = coords(5,:);
C4 = coords(2,:);
NC5 = C5-N;
NC2 = C2-N;
xaxis1 = NC5/norm(NC5); 
NC2 = NC2/norm(NC2);
zaxis1 = cross(xaxis1,NC2);
zaxis1 = zaxis1/norm(zaxis1);
yaxis1 = cross(zaxis1,xaxis1);  
Y = N+norm(NC5)*yaxis1;
Z = N+norm(NC5)*zaxis1;
vectarrow(N,C5); hold on
vectarrow(N,Y); hold on
vectarrow(N,Z); hold on
theta = 27/360*2*pi; 
[xaxis2 yaxis2 zaxis2] = drehnung(xaxis1,yaxis1,zaxis1,theta);
M1 = C5+0.5*(C4-C5);
Nullpunkt = M1+0.688*norm(C4-C5)*xaxis1;
X_1=Nullpunkt+(3.959*norm(C4-C5)+0.7)*xaxis2;
Y_1=Nullpunkt+norm(NC5)*yaxis2;
Z_1=Nullpunkt+norm(NC5)*zaxis2;
vectarrow(Nullpunkt,X_1); hold on
vectarrow(Nullpunkt,Y_1); hold on
vectarrow(Nullpunkt,Z_1); hold on