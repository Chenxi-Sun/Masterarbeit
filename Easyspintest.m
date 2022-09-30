% clear, clf
% 
% % Parameters
% %-------------------------------------------------------------
% g = 2.005;
% A_N = 40.40 + [70,20,20];
% A_mH = [1 1 1]*3.16;
% lw = 0.01;
% tcorr = 8e-11;
% 
% Sys.g = g;
% Sys.Nucs = '14N,1H';
% Sys.A = [A_N;A_mH];
% Sys.lw = [0 lw];
% Exp.mwFreq = 9.5;
% Exp.nPoints = 1e4;
% 
% % Fast-motion regime spectrum
% %-------------------------------------------------------------
% Sys.tcorr = tcorr;
% garlic(Sys,Exp);

Sys.g = 2;
Sys.Nucs = '1H,14N';
Sys.n = [1 1];
Sys.A = [10 11];
Sys.lwpp = [0 0.01];
Exp.mwFreq = 9.7;
Exp.CenterSweep = [346.5 2.8];
garlic(Sys,Exp);