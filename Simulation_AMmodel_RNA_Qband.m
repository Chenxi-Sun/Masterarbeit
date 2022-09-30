    for nr=8:15
% clear;

%%plot PELDOR or distance distribution
str_2='Distance';
%model A or B or C
% model='A';
% model='B';
model='C';

Band = "Qband";
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
zeit = real(RNA.RNA1_8.Qband.time);
Experimental.Sexp = RNA.RNA1_8.Xband.Sexp;
nr=8;
t_max=1100;
Exp(:,1)=RNA.RNA1_8.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_8.Qband.distr(:,2); 
 case '9'
zeit = real(RNA.RNA1_9.Qband.time);
Experimental.Sexp = RNA.RNA1_9.Xband.Sexp;
nr=9;
t_max=1100;
Exp(:,1)=RNA.RNA1_9.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_9.Qband.distr(:,2); 
 case '10'
zeit = real(RNA.RNA1_10.Qband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp;
nr=10;
t_max=1200;
Exp(:,1)=RNA.RNA1_10.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_10.Qband.distr(:,2); 
 case '11'
zeit = real(RNA.RNA1_11.Qband.time);
Experimental.Sexp = RNA.RNA1_11.Xband.Sexp;
nr=11;
t_max=1600;
Exp(:,1)=RNA.RNA1_11.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_11.Qband.distr(:,2); 
 case '12'
zeit = real(RNA.RNA1_12.Qband.time);
Experimental.Sexp = RNA.RNA1_12.Xband.Sexp;
nr=12;
t_max=2000;
Exp(:,1)=RNA.RNA1_12.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_12.Qband.distr(:,2); 
 case '13'
zeit = real(RNA.RNA1_13.Qband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp;
nr=13;
t_max=2200;
Exp(:,1)=RNA.RNA1_13.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_13.Qband.distr(:,2); 
 case '14'
zeit = real(RNA.RNA1_14.Qband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp;
nr=14;
t_max=2200;
Exp(:,1)=RNA.RNA1_14.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_14.Qband.distr(:,2); 
 case '15'
zeit = real(RNA.RNA1_15.Qband.time);
Experimental.Sexp = RNA.RNA1_15.Xband.Sexp;
nr=15;
t_max=2500;
Exp(:,1)=RNA.RNA1_15.Qband.distr(:,1); 
Exp(:,2)=RNA.RNA1_15.Qband.distr(:,2); 
end

sigma_y=6;

%%CmRNA

% [Simulated,R_mean,FWHM,ymax] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,model);        
[Simulated,R_mean,FWHM,ymax] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit); 
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


savefig(F,['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Qband\Distr_Qband_CmARNAModel',model,'1_',str,'.fig'])
savePDF(['D:\ChSun\Result_24.08\CmARNA\CmARNA_Model',model,'_Qband\'],['Distr_Qband_CmARNAModel',model,'1_',str,'.pdf'])


clear;
clf;
end 