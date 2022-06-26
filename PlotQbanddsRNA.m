%%%original data 
RNA1_12 = importdata('RNA1_12_old.mat');
CmdotdsRNA.Qband.Spectrum70MHz = importdata('1_Cmdot_dsRNA_PELDOR_R_bckg.dat');
figure(1)
plot(CmdotdsRNA.Qband.Spectrum70MHz(:,1)*1000,CmdotdsRNA.Qband.Spectrum70MHz(:,2),'b','LineWidth',2)
hold on 
plot(CmdotdsRNA.Qband.Spectrum70MHz(:,1)*1000,CmdotdsRNA.Qband.Spectrum70MHz(:,3),'color',[0.9100    0.4100    0.1700],'LineWidth',2)
plot(RNA1_12.Qband.time*1000,RNA1_12.Qband.raw(:,2),'b','LineWidth',2)
plot(RNA1_12.Qband.time*1000,RNA1_12.Qband.raw(:,3),'color',[0.9100    0.4100    0.1700],'LineWidth',2)
axis([0 1800 0.5 1.05])
xlabel('Time [ns]')
ylabel('Signal intensity ')
title('Q-band original PELDOR data of dsRNA1\_12 molecules')
hold off 

%%%simulation
FF = "DESREF"; %OL3, DESREF, BSC1 (even, odd)
% FF = "OL3";
Band = "Qband";

% Zahl = [12]; def = ["even"];k=1;
CmdotdsRNA.Qband.EP.Sequence = importdata('Sequence.mat'); %Sequence ist eine variable die du oeffnen musst
CmdotdsRNA.Qband.EP.Settings.PumpFrequency=34; %in GHz %entweder so lassen oder mit easyspin simulieren
CmdotdsRNA.Qband.EP.Settings.DetectionFrequency=33.93; %in GHz 
CmdotdsRNA.Qband.EP.Settings.B0=1.20925; %in T 
zeit = real(CmdotdsRNA.Qband.Spectrum70MHz(:,1));
zeit2 = RNA1_12.Qband.time;

Experimental.Sexp(:,1) = CmdotdsRNA.Qband.Spectrum70MHz(:,2)-CmdotdsRNA.Qband.Spectrum70MHz(:,3);
Experimental.Sexp(:,1) = Experimental.Sexp(:,1)+(1-max(Experimental.Sexp(:,1)));
% Experimental.Sexp(:,2) = RNA1_12.Qband.Sexp;
Experimentalold.Sexp = RNA1_12.Qband.Sexp;
%MD simulation von RNA
MD1 = sprintf('BP1_RNA_%s_even.txt',FF);  %Data from txt
MD2 = sprintf('BP12_RNA_%s_even.txt',FF);

Simulated = Spinlabelgeometry(MD1,MD2,CmdotdsRNA.Qband.EP,zeit);
Simulated2 = autoneurechner(MD1,MD2,CmdotdsRNA.Qband.EP,zeit2)

[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
[devn2,SC2] = ScaleModdev('alle',Experimentalold.Sexp,Simulated2.Sexp);

figure(2)
% plot(RNA1_12.Qband.time*1000,RNA1_12.Qband.Sexp,'LineWidth',3)
plot(zeit*1000,Experimental.Sexp,'k','LineWidth',2)
hold on 
plot(zeit2*1000,Experimentalold.Sexp,'k','LineWidth',2)
% plot(zeit*1000,SC,'r','LineWidth',3)
% plot(zeit2*1000,SC2,'r','LineWidth',3)
hold off
% title('Q-band exp. PELDOR data of dsRNA1\_12 molecules')
title('Q-band exp. PELDOR data of dsRNA1\_12 molecules')
% title('PELDOR data recorded at X-Band of Cdot\_dsDNA1\_9 molecules')
xlabel('Time [ns]');
ylabel('Signal intensity')
axis([0 1700 0.75 1.05]);



