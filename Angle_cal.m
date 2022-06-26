m=20;
m2=20;

NC5   = C5_helixA_evenRNA-N1_helixA_evenRNA;
NC2   = C2_helixA_evenRNA-N1_helixA_evenRNA;

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

M = C5_helixA_evenRNA + 0.5*(C4_helixA_evenRNA-C5_helixA_evenRNA) + xaxis1*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;

NC5_2   = C5_helixB_evenRNA-N1_helixB_evenRNA;
NC2_2   = C2_helixB_evenRNA-N1_helixB_evenRNA;

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

M2 = C5_helixB_evenRNA + 0.5*(C4_helixB_evenRNA-C5_helixB_evenRNA) + xaxis2*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;