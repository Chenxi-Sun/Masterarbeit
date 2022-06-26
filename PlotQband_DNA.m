%%%original data 
DNA1_9 = importdata('DNA1_9.mat');
CdotdsDNA.Qband.Spectrum70MHz = importdata('CdotdsDNA_4PPeldor_bckg.dat');

%%%simulation
FF = "BSC1"; %OL3, DESREF, BSC1 (even, odd)
% FF = "OL3";
Band = "Qband";

% Zahl = [12]; def = ["even"];k=1;
CdotdsDNA.Qband.EP.Sequence = importdata('Sequence.mat'); %Sequence ist eine variable die du oeffnen musst
CdotdsDNA.Qband.EP.Settings.PumpFrequency=34; %in GHz %entweder so lassen oder mit easyspin simulieren
CdotdsDNA.Qband.EP.Settings.DetectionFrequency=33.93; %in GHz 
CdotdsDNA.Qband.EP.Settings.B0=1.20925; %in T 
zeit = real(CdotdsDNA.Qband.Spectrum70MHz(:,1));
zeit2 = DNA1_9.Qband.time;

Experimental.Sexp(:,1) = CdotdsDNA.Qband.Spectrum70MHz(:,2)-CdotdsDNA.Qband.Spectrum70MHz(:,3);
Experimental.Sexp(:,1) = Experimental.Sexp(:,1)+(1-max(Experimental.Sexp(:,1)));
% Experimental.Sexp(:,2) = DNA1_9.Qband.Sexp;
Experimentalold.Sexp = DNA1_9.Qband.Sexp;
%MD simulation von RNA
MD1 = sprintf('BP1_DNA_%s_odd.txt',FF);  %Data from txt
MD2 = sprintf('BP9_DNA_%s_odd.txt',FF);

Simulated = Spinlabelgeometry(MD1,MD2,CdotdsDNA.Qband.EP,zeit);
Simulated2 = autoneurechner(MD1,MD2,CdotdsDNA.Qband.EP,zeit2);

[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
[devn2,SC2] = ScaleModdev('alle',Experimentalold.Sexp,Simulated2.Sexp);

figure(2)
% plot(DNA1_9.Qband.time*1000,DNA1_9.Qband.Sexp,'LineWidth',3)
plot(zeit*1000,Experimental.Sexp,'k','LineWidth',2)
hold on
plot(zeit2*1000,Experimentalold.Sexp,'k','LineWidth',2)
plot(zeit*1000,SC,'b','LineWidth',3)
plot(zeit2*1000,SC2,'b','LineWidth',3)
hold off
title('Q-band exp. PELDOR data of dsDNA1\_9 molecules')
% title('Q-band simulated PELDOR data of dsDNA1\_9 molecules')
% title('PELDOR data recorded at X-Band of Cdot\_dsDNA1\_9 molecules')
xlabel('Time [ns]');
ylabel('Signal intensity')
axis([0 1700 0.5 1.05]);