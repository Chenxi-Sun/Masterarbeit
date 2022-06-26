function Result = autoneurechner(MD1,MD2,Param,zeit)

%finding the principle axis of the Cspin-label from MD data of cytosin
%atoms. Also the position of cooridinate origin of the PA-system (middle of
%N-O-bond)


% format of MD file: N, C2, C4, C5; x,y,z for each one; one row=one pico sec
A = importdata(MD1);

[m,n]=size(A);
N  = [A(:,1) A(:,2) A(:,3)];
C2 = [A(:,4) A(:,5) A(:,6)];
C4 = [A(:,7) A(:,8) A(:,9)];
C5 = [A(:,10) A(:,11) A(:,12)];

NC5   = C5-N;
NC2   = C2-N;

xaxis1 = zeros(m,3);
yaxis1 = zeros(m,3);
zaxis1 = zeros(m,3);

for ra=1:m
    xaxis1(ra,:) = NC5(ra,:)/norm(NC5(ra,:));
    NC2(ra,:) = NC2(ra,:)/norm(NC2(ra,:));
    zaxis1(ra,:) = cross(xaxis1(ra,:),NC2(ra,:));
    zaxis1(ra,:) = zaxis1(ra,:)/norm(zaxis1(ra,:));
    yaxis1(ra,:) = cross(zaxis1(ra,:),xaxis1(ra,:));
end

% M = C5 + 0.5*(C4-C5) + xaxis1.*(sqrt(7^2-0.5*vecnorm((C4-C5)')'.^2)+sqrt(8.1^2-0.5*vecnorm((C4-C5)')'.^2))/2;
M = C5 + 0.5*(C4-C5) + xaxis1*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;
% %7: C4 zu N von N-O und 8.2: C4 zu O von N-O
% M = C5 + 0.5*(C4-C5) + 2*(C4-C2) + xaxis1*sqrt(2.3^2-0.71^2) + 0.6*xaxis1 %2.3 ist quer durch 5ring aus pymol, 0.71 ist hälfte von C4-C5 (durchschnitt), 0.6 (halbe N-O aus pymol)
% M = C5 + 0.5*(C4-C5) + 6*(1.5*sin(pi/3)) + 0.5*1.4 %1.5 für C-C bindung. fünf=sechsring angenommen. N-O 1.4
MN = C5 + 0.5*(C4-C5) + xaxis1*(sqrt(7^2-0.75^2));
MO = C5 + 0.5*(C4-C5) + xaxis1*(sqrt(8.2^2-0.75^2));


A2 = importdata(MD2);

[m2,n2]=size(A2);
N_2  = [A2(:,1) A2(:,2) A2(:,3)];
C2_2 = [A2(:,4) A2(:,5) A2(:,6)];
C4_2 = [A2(:,7) A2(:,8) A2(:,9)];
C5_2 = [A2(:,10) A2(:,11) A2(:,12)];

NC5_2   = C5_2-N_2;
NC2_2   = C2_2-N_2;

xaxis2 = zeros(m2,3);
yaxis2 = zeros(m2,3);
zaxis2 = zeros(m2,3);

for ra2=1:m2
    xaxis2(ra2,:) = NC5_2(ra2,:)/norm(NC5_2(ra2,:));
    NC2_2(ra2,:) = NC2_2(ra2,:)/norm(NC2_2(ra2,:));
    zaxis2(ra2,:) = cross(xaxis2(ra2,:),NC2_2(ra2,:));
    zaxis2(ra2,:) = zaxis2(ra2,:)/norm(zaxis2(ra2,:));
    yaxis2(ra2,:) = cross(zaxis2(ra2,:),xaxis2(ra2,:));
end

% M2 = C5_2 + 0.5*(C4_2-C5_2) + xaxis2.*(sqrt(7^2-0.5*vecnorm((C4_2-C5_2)')'.^2)+sqrt(8.1^2-0.5*vecnorm((C4_2-C5_2)')'.^2))/2;
M2 = C5_2 + 0.5*(C4_2-C5_2) + xaxis2*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;
% M2 = C5_2 + 0.5*(C4_2-C5_2) + 6*(1.5*sin(pi/3)) + 0.5*1.4 %1.5 für C-C bindung. fünf=sechsring angenommen. N-O 1.4

MN_2 = C5_2 + 0.5*(C4_2-C5_2) + xaxis2*(sqrt(7^2-0.75^2));
MO_2 = C5_2 + 0.5*(C4_2-C5_2) + xaxis2*(sqrt(8.2^2-0.75^2));

R = M2-M;
% R = R*1.11;

[o1,o2,r]=RuCoordv1(R,xaxis1,yaxis1,zaxis1,xaxis2,yaxis2,zaxis2);


Conformers.EulerAngles.R1=[o1(:,:)];

Conformers.EulerAngles.R2=[o2(:,:)];

Conformers.Distance = r/10;

Result=MainPELDORtime(Param,Conformers,zeit); 

