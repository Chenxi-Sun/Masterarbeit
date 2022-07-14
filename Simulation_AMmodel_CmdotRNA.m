clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');

% Band = "Gband";
% EP = importdata('ExpParG_RNA1_10.mat');
dsRNA.Xband.EP.Sequence = EP.Sequence;
dsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;

% load('E:\Vorlesungen\EPR\Master\ChSun\Cmdot_RNA_Xband\RNA1_12_new.mat')
load('Z:\Students\ChSun\Cmdot_RNA_Xband\RNA1_12_new.mat')
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\RNA1_12dot.mat')

%Xband
Experimental.Sexp(:,1)=CmdotdsRNA.Xband.Spektrum40(1:151,2);
Experimental.Sexp(:,2)=CmdotdsRNA.Xband.Spektrum50(1:151,2);
Experimental.Sexp(:,3)=CmdotdsRNA.Xband.Spektrum60(1:151,2);
Experimental.Sexp(:,4)=CmdotdsRNA.Xband.Spektrum70(1:151,2);
Experimental.Sexp(:,5)=CmdotdsRNA.Xband.Spektrum80(1:151,2);
Experimental.Sexp(:,6)=CmdotdsRNA.Xband.Spektrum90(1:151,2);
zeit=CmdotdsRNA.Xband.Spektrum40(:,1);

%Gband
% Experimental.Sexp(:,1)=RNA1_12dot.Gband.Sexp(:,1);
% Experimental.Sexp(:,2)=RNA1_12dot.Gband.Sexp(:,2);
% Experimental.Sexp(:,3)=RNA1_12dot.Gband.Sexp(:,3);
% Experimental.Sexp(:,4)=RNA1_12dot.Gband.Sexp(:,4);
% Experimental.Sexp(:,5)=RNA1_12dot.Gband.Sexp(:,5);
% zeit=RNA1_12dot.Gband.time;


sigma_y=4;
[Simulated,R_mean] = AM_PELDOR_CmdotRNA(sigma_y,12,EP,zeit)



%%plot

[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];



%Gband
% [devn,SC] = ScaleModdev('single',Experimental.Sexp,Simulated.Sexp);
% o = 0.02; 
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+2*o,Experimental.Sexp(:,3)+3.5*o,Experimental.Sexp(:,4)+4.5*o,Experimental.Sexp(:,5)+6*o];
% SC = [SC(:,1),SC(:,2)+2*o,SC(:,3)+o*3.5,SC(:,4)+o*4.5,SC(:,5)+o*6];

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
% axis([0 1792 0.985 1.122]);
axis([0 1792 0.3 1.5])
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

% CmRNAModelA
savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CmRNA_ModelA_Xband\Distr_CmdotRNAModelA1_12.fig'])
savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CmRNA_ModelA_Xband\',['Distr_CmdotRNAModelA1_12.pdf'])

%CmRNAModelB
% savefig(F,['Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\RNA_ModelB_4grad\CmdotRNAModelB1_12.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\Seminar_06.07_CS\RNA_ModelB_4grad\',['CmdotRNAModelB1_12.pdf'])