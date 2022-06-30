% for nr=8:14
% clear;
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
% for nr=8:15
%     
% str=num2str(nd);

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
sigma_y=8;
t_max=1100;
 case '10'
zeit = real(RNA.RNA1_10.Xband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp(:,1:6);
nr=10;
sigma_y=8;
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
sigma_y=8;
t_max=2000;
 case '13'
zeit = real(RNA.RNA1_13.Xband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp(:,1:6);
nr=13;
sigma_y=2;
t_max=2500;
 case '14'
zeit = real(RNA.RNA1_14.Xband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp(:,1:6);
nr=14;
sigma_y=0;
t_max=2200;
 case '15'
zeit = real(RNA.RNA1_14.Xband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp(:,1:6);
nr=15;
sigma_y=0;
t_max=2200;
end
str2=num2str(sigma_y);

F=figure(1)
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [1 2 3 4 5 6]);
set(gca,'linewidth',1.5) 
xlabel('Distance [nm]')
ylabel('no. of occurences')
xlim([1 5])

% title(['Çm RNA1-',str,' (',str2,'^o)']);
% title(['Cm RNA1-',str,' (',str2,'^o)']);
set(gca,'linewidth',1.5) 


% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\CmRNAModelA1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\',['CmRNAModelA1_',str,'.pdf'])
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\CmdotRNAModelA1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\preAMvsexp_CmRNA_ModelA\',['CmdotRNAModelA1_',str,'.pdf'])