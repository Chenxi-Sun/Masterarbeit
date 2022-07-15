%%CDNA
% nd='Which 2nd position? (5-14):';
% str=input(nd,'s');
% 
% importdata(['C_DNA_1_',str,'_distr.dat']);
% hold on 
% plot(ans(:,1),ymax/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])
% xlim([1 5])

%%CmRNA
% importdata('SumX_CmRNA1_13_distr.dat');
% hold on 
% plot(ans(:,1),67/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])

%%CdotDNA
% importdata('sumX_Cdot_DNA1_9_distr.dat');
% hold on 
% plot(ans(:,1),51/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])

%%CmdotRNA
importdata('SumX_CmdotRNA1_12_distr.dat');
hold on 
plot(ans(:,1),58/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
xlim([1 5])
xlabel('distance [nm]')
ylabel('no. of occurences')
xticks([1 2 3 4 5 6])