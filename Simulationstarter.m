clear;
% FF = "DESREF"; %OL3, DESREF, BSC1 (even, odd)
% % FF = "OL3";
% Band = "Qband";

%  FF = "DESREF"; %OL3, DESREF, BSC1 (even, odd)
FF = "BSC1";
Band = "Xband";

%%%%MD simulation von old RNA
% MD1 = sprintf('BP1_RNA_%s_even.txt',FF);  %Data from txt
% MD2 = sprintf('BP12_RNA_%s_even.txt',FF);
% 
% EP = importdata('ExperimentalParameters.mat');
% RNA1_12.Xband.EP.Sequence = EP.Sequence;
% RNA1_12.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
% RNA1_12.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
% RNA1_12.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
% RNA1_12=importdata('RNA1_12_old.mat');
% zeit = real(RNA1_12.Xband.Spektrum40(:,1));
% 
% Experimental.Sexp(:,1) = RNA1_12.Xband.Spektrum40(:,2)-RNA1_12.Xband.Spektrum40(:,3);
% Experimental.Sexp(:,1) = Experimental.Sexp(:,1)+(1-max(Experimental.Sexp(:,1)));
% 
% Experimental.Sexp(:,2) = RNA1_12.Xband.Spektrum50(:,2)-RNA1_12.Xband.Spektrum50(:,3);
% Experimental.Sexp(:,2) = Experimental.Sexp(:,2)+(1-max(Experimental.Sexp(:,2)));
% 
% Experimental.Sexp(:,3) = RNA1_12.Xband.Spektrum60(:,2)-RNA1_12.Xband.Spektrum60(:,3);
% Experimental.Sexp(:,3) = Experimental.Sexp(:,3)+(1-max(Experimental.Sexp(:,3)));
% 
% Experimental.Sexp(:,4) = RNA1_12.Xband.Spektrum70(:,2)-RNA1_12.Xband.Spektrum70(:,3);
% Experimental.Sexp(:,4) = Experimental.Sexp(:,4)+(1-max(Experimental.Sexp(:,4)));
% 
% Experimental.Sexp(:,5) = RNA1_12.Xband.Spektrum80(:,2)-RNA1_12.Xband.Spektrum80(:,3);
% Experimental.Sexp(:,5) = Experimental.Sexp(:,5)+(1-max(Experimental.Sexp(:,5)));
% 
% Experimental.Sexp(:,6) = RNA1_12.Xband.Spektrum90(:,2)-RNA1_12.Xband.Spektrum90(:,3);
% Experimental.Sexp(:,6) = Experimental.Sexp(:,6)+(1-max(Experimental.Sexp(:,6)));

%%%%MD simulation von RNA
% MD1 = sprintf('BP1_RNA_%s_even.txt',FF);  %Data from txt
% MD2 = sprintf('BP12_RNA_%s_even.txt',FF);
% 
% EP = importdata('ExperimentalParameters.mat');
% CmdotdsRNA.Xband.EP.Sequence = EP.Sequence;
% CmdotdsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
% CmdotdsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
% CmdotdsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
CmdotdsRNA.Xband.Spektrum40 = importdata('PELDOR_dsRNA_40MHz_bckg.dat');
CmdotdsRNA.Xband.Spektrum50 = importdata('PELDOR_dsRNA_50MHz_bckg.dat');
CmdotdsRNA.Xband.Spektrum60 = importdata('PELDOR_dsRNA_60MHz_bckg.dat');
CmdotdsRNA.Xband.Spektrum70 = importdata('PELDOR_dsRNA_70MHz_bckg.dat');
CmdotdsRNA.Xband.Spektrum80 = importdata('PELDOR_dsRNA_80MHz_bckg.dat');
CmdotdsRNA.Xband.Spektrum90 = importdata('PELDOR_dsRNA_90MHz_bckg.dat');
zeit = real(CmdotdsRNA.Xband.Spektrum40(:,1));
% 
% Experimental.Sexp(:,1) = CmdotdsRNA.Xband.Spektrum40(:,2)-CmdotdsRNA.Xband.Spektrum40(:,3);
% Experimental.Sexp(:,1) = Experimental.Sexp(:,1)+(1-max(Experimental.Sexp(:,1)));
% 
% Experimental.Sexp(:,2) = CmdotdsRNA.Xband.Spektrum50(:,2)-CmdotdsRNA.Xband.Spektrum50(:,3);
% Experimental.Sexp(:,2) = Experimental.Sexp(:,2)+(1-max(Experimental.Sexp(:,2)));
% 
% Experimental.Sexp(:,3) = CmdotdsRNA.Xband.Spektrum60(:,2)-CmdotdsRNA.Xband.Spektrum60(:,3);
% Experimental.Sexp(:,3) = Experimental.Sexp(:,3)+(1-max(Experimental.Sexp(:,3)));
% 
% Experimental.Sexp(:,4) = CmdotdsRNA.Xband.Spektrum70(:,2)-CmdotdsRNA.Xband.Spektrum70(:,3);
% Experimental.Sexp(:,4) = Experimental.Sexp(:,4)+(1-max(Experimental.Sexp(:,4)));
% 
% Experimental.Sexp(:,5) = CmdotdsRNA.Xband.Spektrum80(:,2)-CmdotdsRNA.Xband.Spektrum80(:,3);
% Experimental.Sexp(:,5) = Experimental.Sexp(:,5)+(1-max(Experimental.Sexp(:,5)));
% 
% Experimental.Sexp(:,6) = CmdotdsRNA.Xband.Spektrum90(:,2)-CmdotdsRNA.Xband.Spektrum90(:,3);
% Experimental.Sexp(:,6) = Experimental.Sexp(:,6)+(1-max(Experimental.Sexp(:,6)));

%%%%MD simulation von old DNA
MD1 = sprintf('BP1_DNA_%s_odd.txt',FF);  %Data from txt
MD2 = sprintf('BP9_DNA_%s_odd.txt',FF);

EP = importdata('ExperimentalParameters.mat');
dsDNA.Xband.EP.Sequence = EP.Sequence;
dsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;

DNA1_9 = importdata('DNA1_9.mat');

zeit = real(DNA1_9.Xband.time);

Experimental.Sexp = DNA1_9.Xband.Sexp;

%%%MD simulation von DNA
% MD1 = sprintf('BP1_DNA_%s_odd.txt',FF);  %Data from txt
% MD2 = sprintf('BP9_DNA_%s_odd.txt',FF);
% 
% EP = importdata('ExperimentalParameters.mat');
% CdotdsDNA.Xband.EP.Sequence = EP.Sequence;
% CdotdsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
% CdotdsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
% CdotdsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
% CdotdsDNA.Xband.Spektrum40 = importdata('PELDOR_40MHz_bckg.dat');
% CdotdsDNA.Xband.Spektrum50 = importdata('PELDOR_50MHz_bckg.dat');
% CdotdsDNA.Xband.Spektrum60 = importdata('PELDOR_60MHz_bckg.dat');
% CdotdsDNA.Xband.Spektrum70 = importdata('PELDOR_70MHz_Wdh2_bckg.dat');
% CdotdsDNA.Xband.Spektrum80 = importdata('PELDOR_80MHz_bckg.dat');
% CdotdsDNA.Xband.Spektrum90 = importdata('PELDOR_90MHz_Wdh_bckg.dat');
% zeit = real(CdotdsDNA.Xband.Spektrum40(:,1));
% 
% Experimental.Sexp(:,1) = CdotdsDNA.Xband.Spektrum40(:,2)-CdotdsDNA.Xband.Spektrum40(:,3);
% Experimental.Sexp(:,1) = Experimental.Sexp(:,1)+(1-max(Experimental.Sexp(:,1)));
% 
% Experimental.Sexp(:,2) = CdotdsDNA.Xband.Spektrum50(:,2)-CdotdsDNA.Xband.Spektrum50(:,3);
% Experimental.Sexp(:,2) = Experimental.Sexp(:,2)+(1-max(Experimental.Sexp(:,2)));
% 
% Experimental.Sexp(:,3) = CdotdsDNA.Xband.Spektrum60(:,2)-CdotdsDNA.Xband.Spektrum60(:,3);
% Experimental.Sexp(:,3) = Experimental.Sexp(:,3)+(1-max(Experimental.Sexp(:,3)));
% 
% Experimental.Sexp(:,4) = CdotdsDNA.Xband.Spektrum70(:,2)-CdotdsDNA.Xband.Spektrum70(:,3);
% Experimental.Sexp(:,4) = Experimental.Sexp(:,4)+(1-max(Experimental.Sexp(:,4)));
% 
% Experimental.Sexp(:,5) = CdotdsDNA.Xband.Spektrum80(:,2)-CdotdsDNA.Xband.Spektrum80(:,3);
% Experimental.Sexp(:,5) = Experimental.Sexp(:,5)+(1-max(Experimental.Sexp(:,5)));
% 
% Experimental.Sexp(:,6) = CdotdsDNA.Xband.Spektrum90(:,2)-CdotdsDNA.Xband.Spektrum90(:,3);
% Experimental.Sexp(:,6) = Experimental.Sexp(:,6)+(1-max(Experimental.Sexp(:,6)));

%Result 
% Simulated = Spinlabelgeometry(MD1,MD2,CdotdsDNA.Xband.EP,zeit);
Simulated = AMPELDOR_DNA(EP,zeit);
% Simulated = autoneurechner(MD1,MD2,RNA1_12.Xband.EP,zeit);
% Simulated = autoneurechner(MD1,MD2,dsDNA.Xband.EP,zeit);
% Simulated = Spinlabelgeometry(MD1,MD2,CmdotdsRNA.Xband.EP,zeit);

figure(1)
plotmitoffset(Simulated.Sexp,Simulated.time,0.1)
ylim([0.5 1.55])
legend('40MHz','50MHz','60MHz','70MHz','80MHz','90MHz')

% figure
% plot(Simulated.time,Simulated.Sexp)

%%
[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
figure(2)
plot(zeit*1000,SC,'r','LineWidth',3)
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
% hold off
% title('X-band PELDOR simulated data of Cmdot\_dsRNA1\_12 molecules');
% title('X-band PELDOR simulated data of Çm\_dsRNA1\_12 molecules');
% title('X-band PELDOR exp. data of Cmdot\_dsRNA1\_12 molecules');
xlabel('Time [ns]');
axis([0 1700 0.3 1.5]);
title('X-band PELDOR simulated data of Ç\_dsDNA1\_9 molecules');
% title('X-band PELDOR exp. data of Cdot\_dsDNA1\_9 molecules');


