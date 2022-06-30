clear;
% DNA1_9 = importdata('DNA1_9.mat');
% CdotdsDNA.Qband = importdata('CdotdsDNA_4PPeldor_distr.dat');
% plot(CdotdsDNA.Qband(:,1),CdotdsDNA.Qband(:,2),'LineWidth',3)
% hold on 
% plot(DNA1_9.Qband.distr(:,1),DNA1_9.Qband.distr(:,2),'r','LineWidth',3)
% axis([1.4 6.5 0 0.07])
% xlabel('r(nm)')
% ylabel('P(nm^-1)')
% title('Distance distribution of dsDNA1\_9 molecules')
% legend('Cdot labeled','Ç labeled')

% [ymax1,p1] = max(CdotdsDNA.Qband(:,2));
% xmax1=CdotdsDNA.Qband(p1,1);
% 
% [ymax2,p2] = max(DNA1_9.Qband.distr(:,2));
% xmax2=DNA1_9.Qband.distr(p2,1);


RNA1_12 = importdata('RNA1_12_old.mat');
CmdotdsRNA.Qband = importdata('1_Cmdot_dsRNA_PELDOR_R_distr.dat');
plot(CmdotdsRNA.Qband(:,1),CmdotdsRNA.Qband(:,2),'LineWidth',3)
hold on 
plot(RNA1_12.Qband.distr(:,1),RNA1_12.Qband.distr(:,2),'r','LineWidth',3)
axis([1.4 6.5 0 0.02])
xlabel('r(nm)')
ylabel('P(nm^-1)')
title('Distance distribution of dsRNA1\_12 molecules')
legend('Cmdot labeled','Çm labeled') 

[ymax1,p1] = max(CmdotdsRNA.Qband(:,2));
xmax1=CmdotdsRNA.Qband(p1,1);

[ymax2,p2] = max(RNA1_12.Qband.distr(:,2));
xmax2=RNA1_12.Qband.distr(p2,1);

[ymax12,p12] = 0.5*max(CmdotdsRNA.Qband(:,2));
xmax12=CmdotdsRNA.Qband(p12,1);



% CmRNA = importdata('CmRNA_distr.dat');
% CmdotRNA = importdata('CmdotRNA_distr.dat');
% plot(CmdotRNA(:,1),CmdotRNA(:,2),'LineWidth',3)
% hold on 
% plot(CmRNA(:,1),CmRNA(:,2),'r','LineWidth',3)
% xlim([1.4 6.5])
% xlabel('r(nm)')
% ylabel('P(nm^-1)')
% title('Distance distribution of dsRNA1\_12 molecules')
% legend('Cmdot labeled','Çm labeled') 
% 
% [ymax1,p1] = max(CmdotRNA(:,2));
% xmax1=CmdotRNA(p1,1);
% 
% [ymax2,p2] = max(CmRNA(:,2));
% xmax2=CmRNA(p2,1);


% C_DNA = importdata('C_DNA_distr.dat');
% Cdot_DNA = importdata('Cdot_DNA_distr.dat');
% plot(Cdot_DNA(:,1),Cdot_DNA(:,2),'LineWidth',3)
% hold on 
% plot(C_DNA(:,1),C_DNA(:,2),'r','LineWidth',3)
% xlim([1.2 6.5])
% xlabel('r(nm)')
% ylabel('P(nm^-1)')
% title('Distance distribution of dsDNA1\_9 molecules')
% legend('Cdot labeled','Ç labeled')
% 
% [ymax1,p1] = max(Cdot_DNA(:,2));
% xmax1=Cdot_DNA(p1,1);
% 
% [ymax2,p2] = max(C_DNA(:,2));
% xmax2=C_DNA(p2,1);


C_DNA_1_5 = importdata('C_DNA_1_5_distr.dat');
C_DNA_1_6 = importdata('C_DNA_1_6_distr.dat');
C_DNA_1_7 = importdata('C_DNA_1_7_distr.dat');
C_DNA_1_8 = importdata('C_DNA_1_8_distr.dat');
C_DNA_1_9 = importdata('C_DNA_1_9_distr.dat');
C_DNA_1_10 = importdata('C_DNA_1_10_distr.dat');
C_DNA_1_11= importdata('C_DNA_1_11_distr.dat');
C_DNA_1_12 = importdata('C_DNA_1_12_distr.dat');
C_DNA_1_13 = importdata('C_DNA_1_13_distr.dat');
C_DNA_1_14 = importdata('C_DNA_1_14_distr.dat');

fitdist(C_DNA_1_5(:,2),'Normal')

plot(C_DNA_1_5(:,1),C_DNA_1_5(:,2),'LineWidth',2)
hold on 
plot(C_DNA_1_6(:,1),C_DNA_1_6(:,2),'LineWidth',2)
plot(C_DNA_1_7(:,1),C_DNA_1_7(:,2),'LineWidth',2)
plot(C_DNA_1_8(:,1),C_DNA_1_8(:,2),'LineWidth',2)
plot(C_DNA_1_9(:,1),C_DNA_1_9(:,2),'LineWidth',2)
plot(C_DNA_1_10(:,1),C_DNA_1_10(:,2),'LineWidth',2)
plot(C_DNA_1_11(:,1),C_DNA_1_11(:,2),'LineWidth',2)
plot(C_DNA_1_12(:,1),C_DNA_1_12(:,2),'LineWidth',2)
plot(C_DNA_1_13(:,1),C_DNA_1_13(:,2),'LineWidth',2)
plot(C_DNA_1_14(:,1),C_DNA_1_14(:,2),'LineWidth',2)

