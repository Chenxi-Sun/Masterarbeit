% for nd=5:14
% clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsDNA.Xband.EP.Sequence = EP.Sequence;
dsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\dsDNA_Cspin.mat');
% load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_DNA\dsDNA_Cspin.mat');

nd='Which 2nd position? (5-14):';
str=input(nd,'s');

% str=num2str(nd);
t_max=2000;

switch (str)
 case '5'
zeit = real(DNApeldor.S0105.T)./1000;
Experimental.Sexp = DNApeldor.S0105.Sexp;
nr=5;
sigma_y=8;
 case '6'
zeit = real(DNApeldor.S0106.T)./1000;
Experimental.Sexp = DNApeldor.S0106.Sexp;
Experimental.Sexp=Experimental.Sexp(1:185,:);
nr=6;
sigma_y=18;
 case '7'
zeit = real(DNApeldor.S0107.T)./1000;
Experimental.Sexp = DNApeldor.S0107.Sexp;
nr=7;
sigma_y=10;
 case '8'
zeit = real(DNApeldor.S0108.T)./1000;
Experimental.Sexp = DNApeldor.S0108.Sexp;
nr=8;
sigma_y=0;
 case '9'
zeit = real(DNApeldor.S0109.T)./1000;
Experimental.Sexp = DNApeldor.S0109.Sexp;
nr=9;
sigma_y=0;
 case '10'
zeit = real(DNApeldor.S0110.T)./1000;
Experimental.Sexp = DNApeldor.S0110.Sexp;
nr=10;
sigma_y=0;
 case '11'
zeit = real(DNApeldor.S0111.T)./1000;
Experimental.Sexp = DNApeldor.S0111.Sexp;
nr=11;
sigma_y=0;
t_max=1700;
 case '12'
zeit = real(DNApeldor.S0112.T)./1000;
Experimental.Sexp = DNApeldor.S0112.Sexp;
nr=12;
sigma_y=4;
 case '13'
zeit = real(DNApeldor.S0113.T)./1000;
Experimental.Sexp = DNApeldor.S0113.Sexp;
nr=13;
sigma_y=10;
t_max=2000;
 case '14'
zeit = real(DNApeldor.S0114.T)./1000;
Experimental.Sexp = DNApeldor.S0114.Sexp;
nr=14;
sigma_y=12;
t_max=3000;
end

sigma_y=6;

F=figure(1)
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [1 2 3 4 5 6]);
set(gca,'linewidth',1.5) 
xlabel('Distance [nm]')
ylabel('no. of occurences')
xlim([1 5])
% xlim([2 6])
str2=num2str(sigma_y);
% title(['?? DNA1-',str,' (',str2,'^o)']);
title(['C DNA1-',str,' (',str2,'^o)']);
set(gca,'linewidth',1.5) 

% 
savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\Distr_CDNA_ModelB1_',str,'.fig'])
savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\',['Distr_CDNA_ModelB1_',str,'.pdf'])
% savefig(F,['Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\DNA_ModelB_6grad\Distr_CdotModelB1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\DNA_ModelB_6grad\',['Distr_CdotModelB1_',str,'.pdf'])