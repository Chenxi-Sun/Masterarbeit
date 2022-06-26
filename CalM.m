% NC5   = C5_SpinA_evenDNA-N1_A_even;
% NC2   = C2_SpinA_evenDNA-N1_A_even;
% m=20;
% xaxis1 = zeros(m,3);
% yaxis1 = zeros(m,3);
% zaxis1 = zeros(m,3);
% 
% for ra=1:m
%     xaxis1(ra,:) = NC5(ra,:)/norm(NC5(ra,:));
%     NC2(ra,:) = NC2(ra,:)/norm(NC2(ra,:));
%     zaxis1(ra,:) = cross(xaxis1(ra,:),NC2(ra,:));
%     zaxis1(ra,:) = zaxis1(ra,:)/norm(zaxis1(ra,:));
%     yaxis1(ra,:) = cross(zaxis1(ra,:),xaxis1(ra,:));
% end
% 
% M = C5_SpinA_evenDNA + 0.5*(C4_SpinA_evenDNA-C5_SpinA_evenDNA) + xaxis1*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;
% 
% %%
% NC5_2   = C5_SpinB_evenDNA-N1_B_even;
% NC2_2   = C2_SpinB_evenDNA-N1_B_even;
% m=20;
% xaxis1_2 = zeros(m,3);
% yaxis1_2 = zeros(m,3);
% zaxis1_2 = zeros(m,3);
% 
% for ra=1:m
%     xaxis1_2(ra,:) = NC5_2(ra,:)/norm(NC5_2(ra,:));
%     NC2_2(ra,:) = NC2_2(ra,:)/norm(NC2_2(ra,:));
%     zaxis1_2(ra,:) = cross(xaxis1_2(ra,:),NC2_2(ra,:));
%     zaxis1_2(ra,:) = zaxis1_2(ra,:)/norm(zaxis1_2(ra,:));
%     yaxis1_2(ra,:) = cross(zaxis1_2(ra,:),xaxis1_2(ra,:));
% end
% 
% M2 = C5_SpinB_evenDNA + 0.5*(C4_SpinB_evenDNA-C5_SpinB_evenDNA) + xaxis1_2*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;

%%%%%%odd
NC5   = C5_helixA_odd_DNA-N1_helixA_odd_DNA;
NC2   = C2_helixA_odd_DNA-N1_helixA_odd_DNA;
m=20;
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

M = C5_helixA_odd_DNA + 0.5*(C4_helixA_odd_DNA-C5_helixA_odd_DNA) + xaxis1*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;

%%
NC5_2   = C5_helixB_odd_DNA-N1_helixB_odd_DNA;
NC2_2   = C2_helixB_odd_DNA-N1_helixB_odd_DNA;
m=20;
xaxis1_2 = zeros(m,3);
yaxis1_2 = zeros(m,3);
zaxis1_2 = zeros(m,3);

for ra=1:m
    xaxis1_2(ra,:) = NC5_2(ra,:)/norm(NC5_2(ra,:));
    NC2_2(ra,:) = NC2_2(ra,:)/norm(NC2_2(ra,:));
    zaxis1_2(ra,:) = cross(xaxis1_2(ra,:),NC2_2(ra,:));
    zaxis1_2(ra,:) = zaxis1_2(ra,:)/norm(zaxis1_2(ra,:));
    yaxis1_2(ra,:) = cross(zaxis1_2(ra,:),xaxis1_2(ra,:));
end

M2 = C5_helixB_odd_DNA + 0.5*(C4_helixB_odd_DNA-C5_helixB_odd_DNA) + xaxis1_2*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;
