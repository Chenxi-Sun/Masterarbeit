clear;
% Band = "Xband";
Band = "Gband";

switch Band
    case "Gband"
EP = importdata('ExpParG_RNA1_10.mat');
    case "Xband"
EP = importdata('ExperimentalParameters.mat');
end 

dsDNA.Xband.EP.Sequence = EP.Sequence;
dsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;
sigma_y=6;
switch Band
    case "Gband"
load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\DNA1_9dot.mat')
Experimental.Sexp(:,1)=DNA1_9dot.Gband.Sexp(:,1);
Experimental.Sexp(:,2)=DNA1_9dot.Gband.Sexp(:,2);
Experimental.Sexp(:,3)=DNA1_9dot.Gband.Sexp(:,3);
Experimental.Sexp(:,4)=DNA1_9dot.Gband.Sexp(:,4);
Experimental.Sexp(:,5)=DNA1_9dot.Gband.Sexp(:,5);
zeit=DNA1_9dot.Gband.time;

Simulated = AMPELDOR_CdotDNA(sigma_y,9,EP,zeit,'A')

[devn,SC] = ScaleModdev('single',Experimental.Sexp,Simulated.Sexp);
o = 0.02; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+2*o,Experimental.Sexp(:,3)+4.5*o,Experimental.Sexp(:,4)+5.5*o,Experimental.Sexp(:,5)+6*o];
SC = [SC(:,1),SC(:,2)+2*o,SC(:,3)+o*4.5,SC(:,4)+o*5.5,SC(:,5)+o*6];
    case "Xband"
load('CdotDNA1_9_Xband_data.mat')
Experimental.Sexp(:,1)=CdotDNA1_9.Xband40MHz(1:152);
Experimental.Sexp(:,2)=CdotDNA1_9.Xband50MHz(1:152);
Experimental.Sexp(:,3)=CdotDNA1_9.Xband60MHz(1:152);
Experimental.Sexp(:,4)=CdotDNA1_9.Xband70MHz(1:152);
Experimental.Sexp(:,5)=CdotDNA1_9.Xband80MHz(1:152);
Experimental.Sexp(:,6)=CdotDNA1_9.Xband90MHz(1:152);
zeit=CdotDNA1_9.Time(1:152);

Simulated = AMPELDOR_CdotDNA(sigma_y,9,EP,zeit,'B')

[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
end 

F=figure(2)
plot(zeit*1000,SC,'b','LineWidth',2)
% 
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% % 
xlabel('Time [ns]');

switch Band
    case "Gband"
axis([0 1792 0.99 1.12]);
    case "Xband"
axis([0 1792 0.3 1.5])
end 

str=num2str(9);
str2=num2str(sigma_y);
title(['C DNA1-',str,' (',str2,'^o)']);

set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

% saveas(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\ModelA1_',str,'.png'])


savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA\CDNA_ModelA\Gband_CdotDNA_ModelA1_',str,'.fig'])
savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA\CDNA_ModelA\',['Gband_CdotDNA_ModelA1_',str,'.pdf'])