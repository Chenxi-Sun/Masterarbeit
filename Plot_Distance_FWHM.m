load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expDistances.mat')


%%CDNA Model A
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\FWHM_dsDNA_ModelA.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\R_mean_dsDNA_ModelA.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expTiknohovWidth.mat')

%%CDNA Model B
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\helixdy_SL\FWHM_dsDNA_ModelB.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\helixdy_SL\R_mean_dsDNA_ModelB.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\expTiknohovWidth.mat')
% % 
% % 
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

% savefig(F1,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\Mean_Distance_dsDNA.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\','Mean_Distance_dsDNA.pdf')
% 
% savefig(F1,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\Mean_Distance_dsDNA.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\','Mean_Distance_dsDNA.pdf')

% %%CRNA ModelA
load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\R_mean_modelA_dsRNA.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\FWHM_modelA_dsRNA.mat')
load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\RNA_exp_FWHM_Tikhonov.mat')

%%CRNA ModelB
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\R_mean_modelB_dsRNA.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\FWHM_modelB_dsRNA.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\RNA_exp_FWHM_Tikhonov.mat')

F1=figure(1);
plot(8:15,R_mean_all,'Linewidth',3,'Marker','o')
hold on 
plot(8:15,expDistances.RNA.meandist,'Linewidth',3,'Marker','o')
xlabel('Position of 2^{nd} spin label')
ylabel('Distance [nm]')
legend('predicted','exp.')
set(legend,'Location','northwest');
ylim([1.5 4.5])
xlim([8 15])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'linewidth',1.5) 


F2=figure(2);
plot(8:15,FWHM_all,'Linewidth',3,'Marker','o')
hold on 
plot(8:15,FWHM_Tikhonov,'Linewidth',3,'Marker','o')
xlabel('Position of 2^{nd} spin label')
ylabel('FWHM [nm]')
legend('predicted','exp.')
set(legend,'Location','northwest');
ylim([0 1])
xlim([8 15])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'linewidth',1.5) 
title('FWHM')