for nr=8:15
% clear;

%%plot PELDOR or distance distribution
% str_2='PELDOR';
str_2='Distance';
% model A or B or C
model='A';
% model='B';
% model='C';
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsRNA.Xband.EP.Sequence = EP.Sequence;
dsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
load('D:\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
% nd='Which 2nd position? (8-15):';
% str=input(nd,'s');

    
str=num2str(nr);

switch (str)
 case '8'
zeit = real(RNA.RNA1_8.Xband.time);
Experimental.Sexp = RNA.RNA1_8.Xband.Sexp(:,1:6);
nr=8;
t_max=1100;
 case '9'
zeit = real(RNA.RNA1_9.Xband.time);
Experimental.Sexp = RNA.RNA1_9.Xband.Sexp(:,1:6);
nr=9;
t_max=1100;
 case '10'
zeit = real(RNA.RNA1_10.Xband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp(:,1:6);
nr=10;
t_max=1200;
 case '11'
zeit = real(RNA.RNA1_11.Xband.time);
Experimental.Sexp = RNA.RNA1_11.Xband.Sexp(:,1:6);
nr=11;
t_max=1600;
 case '12'
zeit = real(RNA.RNA1_12.Xband.time);
Experimental.Sexp = RNA.RNA1_12.Xband.Sexp(:,1:6);
nr=12;
t_max=2000;
 case '13'
zeit = real(RNA.RNA1_13.Xband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp(:,1:6);
nr=13;
t_max=2200;
 case '14'
zeit = real(RNA.RNA1_14.Xband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp(:,1:6);
nr=14;
t_max=2200;
 case '15'
zeit = real(RNA.RNA1_15.Xband.time);
Experimental.Sexp = RNA.RNA1_15.Xband.Sexp(:,1:6);
nr=15;
t_max=2500;
end

sigma_y=6;

% 
% % % Model A/B but change sigma_r
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:20
% sigma_r=0.5+i*0.02;
% Simulated = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(1:60,1),Experimental.Sexp(1:60,2)+o,Experimental.Sexp(1:60,3)+2*o,Experimental.Sexp(1:60,4)+3*o,Experimental.Sexp(1:60,5)+4*o,Experimental.Sexp(1:60,6)+5*o];
% SC = [SC(1:60,1),SC(1:60,2)+o,SC(1:60,3)+o*2,SC(1:60,4)+o*3,SC(1:60,5)+o*4,SC(1:60,6)+o*5];
% RMSD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RMSD,2)./6;
% % RSMD_alloffset=RMSD(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_r_min=0.5+n*0.02;
% end 

% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r_min,nr,EP,zeit);

% Model A/B but change sigma_h
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:15
% sigma_h=2.5+i*0.1;
% Simulated = AM_PELDOR_RNA(sigma_h,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(1:45,1),Experimental.Sexp(1:45,2)+o,Experimental.Sexp(1:45,3)+2*o,Experimental.Sexp(1:45,4)+3*o,Experimental.Sexp(1:45,5)+4*o,Experimental.Sexp(1:45,6)+5*o];
% SC = [SC(1:45,1),SC(1:45,2)+o,SC(1:45,3)+o*2,SC(1:45,4)+o*3,SC(1:45,5)+o*4,SC(1:45,6)+o*5];
% RMSD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RMSD,2)./6;
% % RSMD_alloffset=RMSD(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_h_min=2.5+n*0.1;
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_h_min,nr,EP,zeit);

% % % 
% sigma_y=0;
% [Simulated,R_mean,FWHM] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,'B');
% [Simulated,R_mean,FWHM] = AM_PELDOR_ApriRNA(sigma_y,nr,EP,zeit);
% [Simulated,R_mean,FWHM,ymax] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,'B');
% [Simulated,R_mean,FWHM,ymax] = AM_PELDOR_ApriRNA(sigma_y,nr,EP,zeit);

%Model C
[Simulated,R_mean,FWHM] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,model);
% [Simulated,R_mean,FWHM] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
R_mean_all(nr-7,:)=R_mean;
FWHM_all(nr-7,:)=FWHM;
end

% sigma_r=0.81;
% % [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);
% 
% 
% % 
% % %PLOT
% switch (str_2) 
%     case 'PELDOR'
% % [Simulated,R_mean,FWHM] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,model);
% [Simulated,R_mean,FWHM] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% o = 0.1; 
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% 
% 
% F=figure(2);
% plot(zeit*1000,SC,'r','LineWidth',2)
% plot(zeit*1000,SC,'Color',[1 0.411764705882353 0.16078431372549],'LineWidth',2)
% plot(zeit*1000,SC,'Color',[1 0.509803921568627 0.811764705882353],'LineWidth',2)
% 
% hold on
% plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% xlabel('Time [ns]');
% ylabel('Signal intensity')
% str2=num2str(sigma_y);
% title(['Çm RNA1-',str,' (',str2,'^o)']);
% axis([0 t_max 0.3 1.5]);
% set(gca,'FontSize',14,'FontWeight','bold','XTick',...
%     [0 500 1000 1500 2000 2500 3000]);
% set(gca,'linewidth',1.5) 
% 
% savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\Xband_CmARNAModel',model,'1_',str,'.fig'])
% savePDF(['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\'],['Xband_CmARNAModel',model,'1_',str,'.pdf'])
% 
% %CmRNA
%     case 'Distance'
% % [Simulated,R_mean,FWHM,ymax] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,model);      
% [Simulated,R_mean,FWHM,ymax] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
% Exp=importdata(['SumX_CmRNA1_',str,'_distr.dat']);  %%X-band 
% 
% hold on 
% plot(Exp(:,1),ymax/max(Exp(:,2)).*Exp(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])
% 
% 
% F=figure(1);
% set(gca,'FontSize',14,'FontWeight','bold','XTick',...
%     [1 2 3 4 5 6]);
% set(gca,'linewidth',1.5) 
% xlabel('Distance [nm]')
% ylabel('no. of occurences')
% if str2num(str)<13
%     xlim([1 5])
% else 
%     xlim([2 6])
% end 
% str2=num2str(sigma_y);
% title(['Çm RNA1-',str,' (',str2,'^o)']);
% set(gca,'linewidth',1.5) 
% 
% 
% savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\Distr_Xband_CmARNAModel',model,'1_',str,'.fig'])
% savePDF(['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\'],['Distr_Xband_CmARNAModel',model,'1_',str,'.pdf'])
% end 
% 
% 
% % % % 
% % % % %CmRNAmodelA
% % savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_ModelA_Xband\Xband_CmARNAModelA1_',str,'.fig'])
% % savePDF('D:\ChSun\Result_24.08\CmARNA\CmARNA_ModelA_Xband\',['Xband_CmARNAModelA1_',str,'.pdf'])
% % % % % 
% % % % 
% clear;
% clf;
% end 
