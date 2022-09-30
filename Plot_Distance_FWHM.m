load('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')


%%CDNA Model A
% load('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\FWHM_dsDNA_ModelA.mat')
% load('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\R_mean_dsDNA_ModelA.mat')
% load('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expTiknohovWidth.mat')
% 
%%CDNA Model B
% load('D:\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\FWHM_ModelB.mat')
% load('D:\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\mean_distance.mat')
% load('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expTiknohovWidth.mat')
% % % % 
% % % % 
% F1=figure(1)
% plot(5:14,R_mean_all,'Linewidth',3,'Marker','o')
% hold on 
% plot(5:14,expDistances.DNA.meandist(1:10),'Linewidth',3,'Marker','o')
% xlabel('Position of 2^{nd} spin label')
% ylabel('Distance [nm]')
% legend('predicted','exp.')
% set(legend,'Location','northwest');
% ylim([1.5 4.5])
% xlim([5 14])
% set(gca,'FontSize',14,'FontWeight','bold');
% set(gca,'linewidth',1.5) 
% 
% 
% F2=figure(2)
% plot(5:14,FWHM_all,'Linewidth',3,'Marker','o')
% hold on 
% plot(5:14,FWHM_Tikhonov,'Linewidth',3,'Marker','o')
% xlabel('Position of 2^{nd} spin label')
% ylabel('FWHM [nm]')
% legend('predicted','exp.')
% set(legend,'Location','northwest');
% ylim([0 1])
% xlim([5 14])
% set(gca,'FontSize',14,'FontWeight','bold');
% set(gca,'linewidth',1.5) 
% title('FWHM')

% savefig(F1,['D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\Mean_Distance_dsDNA.fig'])
% savePDF('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\','Mean_Distance_dsDNA.pdf')
% 
% savefig(F1,['D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\Mean_Distance_dsDNA.fig'])
% savePDF('D:\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\','Mean_Distance_dsDNA.pdf')

% %%CRNA ModelA
% load('D:\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\R_mean_modelA_dsRNA.mat')
% % load('D:\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\FWHM_modelA_dsRNA.mat')
% load('D:\ChSun\Masterarbeit\AMmodel_RNA\RNA_exp_FWHM_Tikhonov.mat')
% 
% %%CRNA ModelB
% % load('D:\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\R_mean_modelB_dsRNA.mat')
% % load('D:\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\FWHM_modelB_dsRNA.mat')
load('D:\ChSun\Masterarbeit\AMmodel_RNA\RNA_exp_FWHM_Tikhonov.mat')
% % 
% F1=figure(1);
% plot(8:15,R_mean_all,'Linewidth',3,'Marker','o')
% % plot(8:15,Distance(:,1),'Linewidth',3,'Marker','o')
% hold on 
% plot(8:15,expDistances.RNA.meandist,'Linewidth',3,'Marker','o')
% xlabel('Position of 2^{nd} spin label')
% ylabel('Distance [nm]')
% legend('predicted','exp.')
% set(legend,'Location','northwest');
% ylim([1.5 4.5])
% xlim([8 15])
% set(gca,'FontSize',14,'FontWeight','bold');
% set(gca,'linewidth',1.5) 


F2=figure(2);
% plot(8:15,FWHM_all,'Linewidth',3,'Marker','o','Color',[0.466666666666667 0.674509803921569 0.188235294117647])
% plot(8:15,FWHM_all,'Linewidth',3,'Marker','o','Color',[0.466666666666667 0.674509803921569 0.188235294117647])
plot(8:15,FWHM_all,'Linewidth',3,'Marker','o','Color',[1 0.509803921568627 0.811764705882353])
hold on 
plot(8:15,FWHM_Tikhonov,'Linewidth',3,'Marker','o','Color',[0.149019607843137 0.149019607843137 0.149019607843137])
% plot(8:15,FWHM_Q,'Linewidth',3,'Marker','o')
xlabel('Position of 2^{nd} spin label')
ylabel('FWHM [nm]')
legend('predicted','exp.')
set(legend,'Location','northwest');
ylim([0 1])
xlim([8 15])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'linewidth',1.5) 
box off