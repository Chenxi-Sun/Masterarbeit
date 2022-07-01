clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsRNA.Xband.EP.Sequence = EP.Sequence;
dsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;

% load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Cmdot_RNA_Xband\RNA1_12_new.mat')
load('Z:\Students\ChSun\Cmdot_RNA_Xband\RNA1_12_new.mat')

Experimental.Sexp(:,1)=CmdotdsRNA.Xband.Spektrum40(1:151,2);
Experimental.Sexp(:,2)=CmdotdsRNA.Xband.Spektrum50(1:151,2);
Experimental.Sexp(:,3)=CmdotdsRNA.Xband.Spektrum60(1:151,2);
Experimental.Sexp(:,4)=CmdotdsRNA.Xband.Spektrum70(1:151,2);
Experimental.Sexp(:,5)=CmdotdsRNA.Xband.Spektrum80(1:151,2);
Experimental.Sexp(:,6)=CmdotdsRNA.Xband.Spektrum90(1:151,2);

zeit=CmdotdsRNA.Xband.Spektrum40(:,1);
% 

sigma_y=4;
[Simulated,R_mean] = AM_PELDOR_CmdotRNA(sigma_y,12,EP,zeit)


% % Model B change sigma_y
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:10
% sigma_y=0+i*2;
% Simulated = AM_PELDOR_CmdotRNA(sigma_y,12,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% sigma(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(sigma,2)./6;  %compare all offset
% % RSMD_alloffset=sigma(:,1);  %only compare 40Mhz
% [m,n]=min(RSMD_alloffset);
% sigma_y_best=n*2;
% [Simulated,R_mean] = AM_PELDOR_CmdotRNA(sigma_y_best,12,EP,zeit);



%%plot
[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];

F=figure(2)
plot(zeit*1000,SC,'r','LineWidth',2)
% plot(zeit*1000,SC,'b','LineWidth',2)
% 
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% % 
xlabel('Time [ns]');


xlabel('Time [ns]');
ylabel('Signal intensity')
str2=num2str(sigma_y);
title(['Cm RNA1-12 (',str2,'^o)']);
axis([0 1800 0.3 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

%CmRNAModelA
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\CmdotRNAModelA1_12.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\AMmodelA\',['CmdotRNAModelA1_12.pdf'])

%CmRNAModelB
savefig(F,['Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\RNA_ModelB_4grad\CmdotRNAModelB1_12.fig'])
savePDF('Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\RNA_ModelB_4grad\',['CmdotRNAModelB1_12.pdf'])