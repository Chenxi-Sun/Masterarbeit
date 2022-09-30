clear all
%Nitroxide solid state spectra
band='G'; %choose X,Q,G

Sys.S = 1/2;
Sys.g =[2.0088 2.0065 2.0027];
Sys.lw = 1;  % einfach geschaetzt damit es aussieht wie experimental spectra
Sys.Nucs = '14N';
Sys.A = [20 20 98];%[5 5 34.5]; 

figure
switch band
    case 'X'
    Exp.mwFreq = 9.5012;
    center = (planck*Exp.mwFreq*10^9/(bmagn*mean(Sys.g(1,:))))*1000;
    Exp.CenterSweep = [center 14]; Exp.nPoints = 1400;
    hohe = [50,30,10];
    
    % Experimental spektrum of RNA 1-8
%     [Bexp,Iexp] = eprload('Z:\Members\Max\dsRNA_Cspin\A_data\Xband\2020_08_10_RNA8\EDFS.DSC');
    
    case 'Q'
    Exp.mwFreq = 34;
    center = (planck*Exp.mwFreq*10^9/(bmagn*mean(Sys.g)))*1000;
    Exp.CenterSweep = [center 15];  Exp.nPoints = 1400;
    hohe = [150,145,140];
    
%     [Bexp,Iexp] = eprload('Z:\Members\Max\newCspin\data\Qband\Cmdot_dsRNA1_12\Cmdot_sdRNA_FS.DSC');
%     Iexp=Iexp-Iexp(1,1);
    
    case 'G'
    Exp.mwFreq = 180;
    center = (planck*Exp.mwFreq*10^9/(bmagn*mean(Sys.g)))*1000;
    Exp.CenterSweep = [center 40];
    hohe = [25,52.5,25];Exp.nPoints = 1400;
        
%     Gspec = importdata('Z:\Members\Max\newCspin\data\G\2021_10_12_Cdot_DNA_specfixed\FS.dat');
%     Bexp = Gspec(:,1); Iexp = Gspec(:,2);
    
end

Exp.Harmonic = 0;
Opt.Output='separate';
% Opt.Transitions = [1 6; 2 5; 3 4]; 
[B,I1,p] = pepper(Sys,Exp,Opt); 

plot(B,I1(1,:)+I1(3,:)+I1(9,:),'LineStyle','-.','Linewidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333])
hold on
plot(B,I1(5,:)+I1(7,:)+I1(8,:),'LineStyle','-.','Linewidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333])
plot(B,I1(2,:)+I1(4,:)+I1(6,:),'LineStyle','-.','Linewidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333])
axis tight
xlabel('magnetic field / mT')
hold on 
Opt.Output='summed';
[B,I1] = pepper(Sys,Exp); 
plot(B,I1,'Color',[0 0.447058823529412 0.741176470588235],'Linewidth',3)
axis tight
xlabel('magnetic field / mT')
set(gca,'FontSize',14,'FontWeight','bold')
set(gca,'linewidth',1.5) 

%gx,Ax
plot((1000*planck*Exp.mwFreq*10^9/(Sys.g(:,1)*bmagn)),hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,1)*bmagn)*1000)+Sys.A(:,1)/2.8/10,hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,1)*bmagn)*1000)-Sys.A(:,1)/2.8/10,hohe(1,1),'+');hold on
%gy,Ay
plot((1000*planck*Exp.mwFreq*10^9/(Sys.g(:,2)*bmagn)),hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,2)*bmagn)*1000)+Sys.A(:,2)/2.8/10,hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,2)*bmagn)*1000)-Sys.A(:,2)/2.8/10,hohe(1,1),'+');hold on
%gz,Az
plot((1000*planck*Exp.mwFreq*10^9/(Sys.g(:,3)*bmagn)),hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,3)*bmagn)*1000)+Sys.A(:,3)/2.8/10,hohe(1,1),'+');hold on
plot((planck*Exp.mwFreq*10^9/(Sys.g(:,3)*bmagn)*1000)-Sys.A(:,3)/2.8/10,hohe(1,1),'+');hold on
legend('m_I=+1','m_I=0','m_I=-1','summed')

% set(gca,'FontSize',14,'FontWeight','bold','XTick',...
%     [1205 1207 1209 1211 1213 1215 1217]);