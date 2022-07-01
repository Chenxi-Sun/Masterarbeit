clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsDNA.Xband.EP.Sequence = EP.Sequence;
dsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;

load('CdotDNA1_9_Xband_data.mat')

Experimental.Sexp(:,1)=CdotDNA1_9.Xband40MHz(1:152);
Experimental.Sexp(:,2)=CdotDNA1_9.Xband50MHz(1:152);
Experimental.Sexp(:,3)=CdotDNA1_9.Xband60MHz(1:152);
Experimental.Sexp(:,4)=CdotDNA1_9.Xband70MHz(1:152);
Experimental.Sexp(:,5)=CdotDNA1_9.Xband80MHz(1:152);
Experimental.Sexp(:,6)=CdotDNA1_9.Xband90MHz(1:152);

zeit=CdotDNA1_9.Time(1:152);
sigma_y=6;
Simulated = AMPELDOR_CdotDNA(sigma_y,9,EP,zeit)

% % Model B change sigma_y
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:10
% sigma_y=0+i*2;
% Simulated = AMPELDOR_CdotDNA(sigma_y,9,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% sigma(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% % RSMD_alloffset=sum(sigma,2)./6;
% RSMD_alloffset=sigma(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_y_best=n*2;
% [Simulated,R_mean,sigma] = AMPELDOR_CdotDNA(sigma_y_best,9,EP,zeit);

[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];

F=figure(2)
% plot(zeit*1000,SC,'r','LineWidth',3)
plot(zeit*1000,SC,'b','LineWidth',2)
% 
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% % 
xlabel('Time [ns]');
axis([0 1700 0.3 1.5]);
str=num2str(9);
str2=num2str(sigma_y);
title(['C DNA1-',str,' (',str2,'^o)']);

% axis([0 2000 0.3 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

% saveas(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\ModelA1_',str,'.png'])


savefig(F,['Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\DNA_ModelB_6grad\CdotModelB1_',str,'.fig'])
savePDF('Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\DNA_ModelB_6grad\',['CdotModelB1_',str,'.pdf'])