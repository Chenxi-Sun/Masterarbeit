% dlmwrite('RNCdot.txt',[time, Exp, zeros(length(time))])

load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')


summe= sum(RNA.RNA1_15.Xband.Sexp(:,1:6),2);

dlmwrite('SumX_CmRNA1_15.DTA',[RNA.RNA1_15.Xband.time, summe, zeros(length(summe),1)])


% Spektrum40 = importdata('PELDOR_dsRNA_40MHz_bckg.dat');
% Spektrum50 = importdata('PELDOR_dsRNA_50MHz_bckg.dat');
% Spektrum60 = importdata('PELDOR_dsRNA_60MHz_bckg.dat');
% Spektrum70 = importdata('PELDOR_dsRNA_70MHz_bckg.dat');
% Spektrum80 = importdata('PELDOR_dsRNA_80MHz_bckg.dat');
% Spektrum90 = importdata('PELDOR_dsRNA_90MHz_bckg.dat');
% zeit = real(Spektrum40(:,1));
% 
% CmdotRNA.Xband(:,1) = Spektrum40(:,2);
% CmdotRNA.Xband(:,2) = Spektrum50(:,2);
% CmdotRNA.Xband(:,3) = Spektrum60(:,2);
% CmdotRNA.Xband(:,4) = Spektrum70(:,2);
% CmdotRNA.Xband(:,5) = Spektrum80(:,2);
% CmdotRNA.Xband(:,6) = Spektrum90(:,2);
% 
% summe= sum(CmdotRNA.Xband,2);
% 
% dlmwrite('SumX_CmdotRNA1_12.DTA',[Spektrum40(:,1), summe, zeros(length(summe),1)])

% C_DNA1_9 = importdata('DNA1_9.mat')
% 
% C_DNA.Xband(:,1) = C_DNA1_9.Xband.Sexp(:,1);
% C_DNA.Xband(:,2) = C_DNA1_9.Xband.Sexp(:,2);
% C_DNA.Xband(:,3) = C_DNA1_9.Xband.Sexp(:,3);
% C_DNA.Xband(:,4) = C_DNA1_9.Xband.Sexp(:,4);
% C_DNA.Xband(:,5) = C_DNA1_9.Xband.Sexp(:,5);
% C_DNA.Xband(:,6) = C_DNA1_9.Xband.Sexp(:,6);
% 
% summe= sum(C_DNA.Xband,2);
% % 
% dlmwrite('SumX_C_DNA1_9.DTA',[C_DNA1_9.Xband.time, summe, zeros(length(summe),1)])

% Spektrum40 = importdata('PELDOR_40MHz_bckg.dat');
% Spektrum50 = importdata('PELDOR_50MHz_bckg.dat');
% Spektrum60 = importdata('PELDOR_60MHz_bckg.dat');
% Spektrum70 = importdata('PELDOR_70MHz_Wdh2_bckg.dat');
% Spektrum80 = importdata('PELDOR_80MHz_bckg.dat');
% Spektrum90 = importdata('PELDOR_90MHz_Wdh_bckg.dat');
% 
% CdotDNA.Xband(:,1) = Spektrum40(:,2);
% CdotDNA.Xband(:,2) = Spektrum50(:,2);
% CdotDNA.Xband(:,3) = Spektrum60(:,2);
% CdotDNA.Xband(:,4) = Spektrum70(:,2);
% CdotDNA.Xband(:,5) = Spektrum80(:,2);
% CdotDNA.Xband(:,6) = Spektrum90(:,2);
% 
% zeit = real(Spektrum40(:,1));
% 
% summe= sum(CdotDNA.Xband,2);
% 
% dlmwrite('sumX_Cdot_DNA1_9.DTA',[zeit, summe, zeros(length(summe),1)])
