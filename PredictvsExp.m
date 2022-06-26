
%%CDNA
% importdata('C_DNA_1_14_distr.dat');
% hold on 
% plot(ans(:,1),48/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])
% % xlim([2 6])

%%CmRNA
% importdata('SumX_CmRNA1_14_distr.dat');
% hold on 
% plot(ans(:,1),56/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])

%%CdotDNA
% importdata('sumX_Cdot_DNA1_9_distr.dat');
% hold on 
% plot(ans(:,1),58/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])

%%CmdotRNA
importdata('SumX_CmdotRNA1_12_distr.dat');
hold on 
plot(ans(:,1),66/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
xlim([1 5])
xlabel('distance [nm]')
ylabel('no. of occurences')
xticks([1 2 3 4 5 6])