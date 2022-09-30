for nr=8:15
% clear;

%%plot PELDOR or distance distribution
% str_2='PELDOR';
str_2='Distance';


model='C';
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

[Simulated,R_mean,FWHM,ymax] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
% [Simulated,R_mean,FWHM] = AM_C_PELDOR_ApriRNA(sigma_y,nr,EP,zeit);
% R_mean_all(nr-7,:)=R_mean;
% FWHM_all(nr-7,:)=FWHM;
% end

%PLOT
switch (str_2) 
    case 'PELDOR'
[Simulated,R_mean,FWHM] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];


F=figure(2);
plot(zeit*1000,SC,'r','LineWidth',2)
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
xlabel('Time [ns]');
ylabel('Signal intensity')
str2=num2str(sigma_y);
title(['Çm RNA1-',str,' (',str2,'^o)']);
axis([0 t_max 0.3 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\Xband_CmARNAModel',model,'1_',str,'.fig'])
savePDF(['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\'],['Xband_CmARNAModel',model,'1_',str,'.pdf'])

%%CmRNA
    case 'Distance'
[Simulated,R_mean,FWHM,ymax] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);        
Exp=importdata(['SumX_CmRNA1_',str,'_distr.dat']);  %%X-band 

hold on 
plot(Exp(:,1),ymax/max(Exp(:,2)).*Exp(:,2),'k','Linewidth',2)
xlim([1 5])
xlabel('distance [nm]')
ylabel('no. of occurences')
xticks([1 2 3 4 5 6])


F=figure(1);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [1 2 3 4 5 6]);
set(gca,'linewidth',1.5) 
xlabel('Distance [nm]')
ylabel('no. of occurences')
if str2num(str)<13
    xlim([1 5])
else 
    xlim([2 6])
end 
str2=num2str(sigma_y);
title(['Çm RNA1-',str,' (',str2,'^o)']);
set(gca,'linewidth',1.5) 


savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\Distr_Xband_CmARNAModel',model,'1_',str,'.fig'])
savePDF(['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Xband\'],['Distr_Xband_CmARNAModel',model,'1_',str,'.pdf'])
end 

% 
% % % % 
% % % % %CmRNAmodelA
% % savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_ModelA_Xband\Xband_CmARNAModelA1_',str,'.fig'])
% % savePDF('D:\ChSun\Result_24.08\CmARNA\CmARNA_ModelA_Xband\',['Xband_CmARNAModelA1_',str,'.pdf'])
% % % % % 
% % % % 
clear;
clf;
end 