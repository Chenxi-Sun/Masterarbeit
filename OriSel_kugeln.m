Sys = struct('S',1/2,'g',[2.0088 2.0065 2.0027]);
Sys = nucspinadd(Sys,'14N',[20 20 98]);

% X-band
% deltaf = [0.09 0.08 0.07 0.06 0.05 0.04 -0.11];
% i=6;
% Exp = struct('mwFreq',9.519,'Field',339.1,'ExciteWidth',16);
% Exp = struct('mwFreq',9.519+deltaf(1,i),'Field',339.1,'ExciteWidth',16);


% Q-band
% Exp = struct('mwFreq',34-0.070,'Field',1209.24,'ExciteWidth',16);
% Exp = struct('mwFreq',34,'Field',1209.24,'ExciteWidth',16);

% G-band
fld = [6403.8 6406.6 6409.1 6415.1 6421.5];
i = 5;
Exp = struct('mwFreq',180,'Field',fld(1,i),'ExciteWidth',7); %excitation width ist geschaetzt


Weights = orisel(Sys,Exp);

figure
Opt = struct('Symmetry','C1');
Opt.nKnots = 200; %geht schneller wenn man weniger knots hat (sieht aber auch haesslicher aus :D ) 
orisel(Sys,Exp,Opt);
set(gca,'visible','off')
set(colorbar,'visible','off')
