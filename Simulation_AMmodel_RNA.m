% for nr=8:15
clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsRNA.Xband.EP.Sequence = EP.Sequence;
dsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
nd='Which 2nd position? (8-15):';
str=input(nd,'s');

%     
% str=num2str(nr);

switch (str)
 case '8'
zeit = real(RNA.RNA1_8.Xband.time);
Experimental.Sexp = RNA.RNA1_8.Xband.Sexp(:,1:6);
nr=8;
sigma_y=8;
t_max=1100;
 case '9'
zeit = real(RNA.RNA1_9.Xband.time);
Experimental.Sexp = RNA.RNA1_9.Xband.Sexp(:,1:6);
nr=9;
sigma_y=6;
t_max=1100;
 case '10'
zeit = real(RNA.RNA1_10.Xband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp(:,1:6);
nr=10;
sigma_y=6;
t_max=1200;
 case '11'
zeit = real(RNA.RNA1_11.Xband.time);
Experimental.Sexp = RNA.RNA1_11.Xband.Sexp(:,1:6);
nr=11;
sigma_y=0;
t_max=1600;
 case '12'
zeit = real(RNA.RNA1_12.Xband.time);
Experimental.Sexp = RNA.RNA1_12.Xband.Sexp(:,1:6);
nr=12;
sigma_y=2;
t_max=2000;
 case '13'
zeit = real(RNA.RNA1_13.Xband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp(:,1:6);
nr=13;
sigma_y=6;
t_max=2500;
 case '14'
zeit = real(RNA.RNA1_14.Xband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp(:,1:6);
nr=14;
sigma_y=0;
t_max=2200;
 case '15'
zeit = real(RNA.RNA1_15.Xband.time);
Experimental.Sexp = RNA.RNA1_15.Xband.Sexp(:,1:6);
nr=15;
sigma_y=0;
t_max=2200;
end

% [Simulated,R_mean] = AM_PELDOR_RNA(nr,EP,zeit);

% Model A/B but change sigma_r
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:15
% sigma_r=0.3+i*0.05;
% Simulated = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% RMSD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% % RSMD_alloffset=sum(RMSD,2)./6;
% RSMD_alloffset=RMSD(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_r_min=0.3+n*0.05;
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);

% % % 
% sigma_y=0;
[Simulated,R_mean,FWHM] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit);
% R_mean_all(nr-7,:)=R_mean;
% FWHM_all(nr-7,:)=FWHM;
% % end

% sigma_r=0.81;
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);

% % Model B change sigma_y
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:10
% sigma_y=0+i*2;
% Simulated = AM_PELDOR_RNA(sigma_y,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% RSMD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RSMD,2)./6;
% % RSMD_alloffset=RSMD(:,4);
% [m,n]=min(RSMD_alloffset);
% sigma_y_best=n*2;
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_y_best,nr,EP,zeit);


%PLOT
% 
[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];


F=figure(2);
% F=open('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\PELDOR_result\1_9_B.fig')
plot(zeit*1000,SC,'r','LineWidth',2)
% plot(zeit*1000,SC,'b','LineWidth',2)
% 
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% % 
xlabel('Time [ns]');
ylabel('Signal intensity')
str2=num2str(sigma_y);
title(['Çm RNA1-',str,' (',str2,'^o)']);
axis([0 t_max 0.3 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

%CmRNAModelB
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\CmRNAModelA1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\',['CmRNAModelA1_',str,'.pdf'])
% 
%CmRNAmodelA
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\CmRNAModelB1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelB\',['CmRNAModelB1_',str,'.pdf'])
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\CmRNAModelA1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\',['CmRNAModelA1_',str,'.pdf'])


% clear;
% clf;
% end 
