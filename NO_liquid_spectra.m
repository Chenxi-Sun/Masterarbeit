clear;

Nitroxide.g = 2.0053;
Nitroxide.Nucs = '14N';
Nitroxide.n = [1];
aisoNitroxide=1.4;                % hyperfine coupling in mT
Nitroxide.A = mt2mhz(aisoNitroxide,gfree);    % conversion from mT to MHz units
Nitroxide.lwpp = [0,0.02];         % linewidth in mT   
Exp.mwFreq = 9.5;           % mw frequency in GHz
Exp.nPoints = 8192;         % resolution of cw-EPR spectra  
Exp.Harmonic = 1;             %absorption spectrum
[B,spec]=garlic(Nitroxide,Exp);     % routine to generate spectra 

plot(B,spec,'b','LineWidth',2,'Color',[0 0.447058823529412 0.741176470588235])
xlabel('magnetic field [mT]','FontSize',15)
ylabel('first derivative cw-EPR signal','FontSize',15)
box off
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [336 337 338 339 340 341]);
