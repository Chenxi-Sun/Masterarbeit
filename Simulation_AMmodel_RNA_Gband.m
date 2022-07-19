% for nr=8:15
% clear;
Band = "Gband";
EP = importdata('ExpParG_RNA1_10.mat');
dsRNA.Xband.EP.Sequence = EP.Sequence;
dsRNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsRNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsRNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
% load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
nd='Which 2nd position? (8-15):';
str=input(nd,'s');
% % 
%     
% str=num2str(nr);

switch (str)
 case '8'
zeit = real(RNA.RNA1_8.Gband.time);
Experimental.Sexp = RNA.RNA1_8.Gband.Sexp(:,1:5);
nr=8;
t_max=1100;
 case '9'
zeit = real(RNA.RNA1_9.Gband.time);
Experimental.Sexp = RNA.RNA1_9.Gband.Sexp(:,1:5);
nr=9;
t_max=1100;
 case '10'
zeit = real(RNA.RNA1_10.Gband.time);
Experimental.Sexp = RNA.RNA1_10.Gband.Sexp(:,1:5);
nr=10;
t_max=1200;
 case '11'
zeit = real(RNA.RNA1_11.Gband.time);
Experimental.Sexp = RNA.RNA1_11.Gband.Sexp(:,1:5);
nr=11;
t_max=1600;
 case '12'
zeit = real(RNA.RNA1_12.Gband.time);
Experimental.Sexp = RNA.RNA1_12.Gband.Sexp(:,1:5);
nr=12;
t_max=2000;
 case '13'
zeit = real(RNA.RNA1_13.Gband.time);
Experimental.Sexp = RNA.RNA1_13.Gband.Sexp(:,1:5);
nr=13;
t_max=2500;
 case '14'
zeit = real(RNA.RNA1_14.Gband.time);
Experimental.Sexp = RNA.RNA1_14.Gband.Sexp(:,1:5);
nr=14;
t_max=2200;
 case '15'
zeit = real(RNA.RNA1_15.Gband.time);
Experimental.Sexp = RNA.RNA1_15.Gband.Sexp(:,1:5);
nr=15;
t_max=2200;
end
% Model A/B but change sigma_r
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:20
% sigma_r=0.6+i*0.02;
% Simulated = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);
% [devn,SC] = ScaleModdev('single',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(1:25,1),Experimental.Sexp(1:25,2)+o*2,Experimental.Sexp(1:25,3)+3*o,Experimental.Sexp(1:25,4)+4*o,Experimental.Sexp(1:25,5)+5*o];
% SC = [SC(1:25,1),SC(1:25,2)+o*2,SC(1:25,3)+o*3,SC(1:25,4)+o*4,SC(1:25,5)+o*5];
% RMSD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RMSD,2)./5;
% % RSMD_alloffset=RMSD(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_r_min(nr-7)=0.6+n*0.02;
% % [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r_min,nr,EP,zeit);
% end 


% Model A/B but change sigma_h
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:10
% sigma_h=4+i*0.1;
% Simulated = AM_PELDOR_RNA(sigma_h,nr,EP,zeit);
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(1:25,1),Experimental.Sexp(1:25,2)+o*2,Experimental.Sexp(1:25,3)+3*o,Experimental.Sexp(1:25,4)+4*o,Experimental.Sexp(1:25,5)+5*o];
% SC = [SC(1:25,1),SC(1:25,2)+o*2,SC(1:25,3)+o*3,SC(1:25,4)+o*4,SC(1:25,5)+o*5];
% RMSD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RMSD,2)./6;
% % RSMD_alloffset=RMSD(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_h_min(nr-7)=4+n*0.1;
% end
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_h_min,nr,EP,zeit);

% 
sigma_y=0;
% [Simulated,R_mean,FWHM] = AM_PELDOR_RNA(sigma_y,nr,EP,zeit,'B');
% [Simulated,R_mean,FWHM] = AM_PELDOR_ApriRNA(sigma_y,nr,EP,zeit,'B');

[Simulated,R_mean,FWHM] = AM_C_PELDOR_RNA(sigma_y,nr,EP,zeit);
% R_mean_all(nr-7,:)=R_mean;
% FWHM_all(nr-7,:)=FWHM;
% end

% sigma_r=0.81;
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_r,nr,EP,zeit);

% Model B change sigma_y

% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:6
% sigma_y=0+i*1;
% Simulated = AM_PELDOR_RNA(sigma_y,nr,EP,zeit);
% [devn,SC] = ScaleModdev('single',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(1:25,1),Experimental.Sexp(1:25,2)+o*2,Experimental.Sexp(1:25,3)+3*o,Experimental.Sexp(1:25,4)+4*o,Experimental.Sexp(1:25,5)+5*o];
% SC = [SC(1:25,1),SC(1:25,2)+o*2,SC(1:25,3)+o*3,SC(1:25,4)+o*4,SC(1:25,5)+o*5];
% RSMD(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(RSMD,2)./6;
% % RSMD_alloffset=RSMD(:,4);
% [m,n]=min(RSMD_alloffset);
% sigma_y_best(nr-7)=n*1;
% end 
% [Simulated,R_mean] = AM_PELDOR_RNA(sigma_y_best,nr,EP,zeit);


% % PLOT
% 
[devn,SC] = ScaleModdev('single',Experimental.Sexp,Simulated.Sexp);
o = 0.02; 

switch (str)
 case '8'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*2,Experimental.Sexp(:,3)+3*o,Experimental.Sexp(:,4)+4*o,Experimental.Sexp(:,5)+5*o];
SC = [SC(:,1),SC(:,2)+o*2,SC(:,3)+o*3,SC(:,4)+o*4,SC(:,5)+o*5];
 case '9'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*1,Experimental.Sexp(:,3)+4*o,Experimental.Sexp(:,4)+5*o,Experimental.Sexp(:,5)+6*o];
SC = [SC(:,1),SC(:,2)+o*1,SC(:,3)+o*4,SC(:,4)+o*5,SC(:,5)+o*6];
 case '10'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*3,Experimental.Sexp(:,3)+4*o,Experimental.Sexp(:,4)+5*o,Experimental.Sexp(:,5)+6*o];
SC = [SC(:,1),SC(:,2)+o*3,SC(:,3)+o*4,SC(:,4)+o*5,SC(:,5)+o*6];
 case '11'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*1,Experimental.Sexp(:,3)+3*o,Experimental.Sexp(:,4)+4*o,Experimental.Sexp(:,5)+5*o];
SC = [SC(:,1),SC(:,2)+o*1,SC(:,3)+o*3,SC(:,4)+o*4,SC(:,5)+o*5];
 case '12'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*1,Experimental.Sexp(:,3)+3*o,Experimental.Sexp(:,4)+4*o,Experimental.Sexp(:,5)+5*o];
SC = [SC(:,1),SC(:,2)+o*1,SC(:,3)+o*3,SC(:,4)+o*4,SC(:,5)+o*5];
 case '13'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*1,Experimental.Sexp(:,3)+3*o,Experimental.Sexp(:,4)+4*o,Experimental.Sexp(:,5)+6*o];
SC = [SC(:,1),SC(:,2)+o*1,SC(:,3)+o*3,SC(:,4)+o*4,SC(:,5)+o*6];
 case '14'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*1,Experimental.Sexp(:,3)+3*o,Experimental.Sexp(:,4)+4*o,Experimental.Sexp(:,5)+5*o];
SC = [SC(:,1),SC(:,2)+o*1,SC(:,3)+o*3,SC(:,4)+o*4,SC(:,5)+o*5];
 case '15'
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o*3,Experimental.Sexp(:,3)+4*o,Experimental.Sexp(:,4)+5*o,Experimental.Sexp(:,5)+6*o];
SC = [SC(:,1),SC(:,2)+o*3,SC(:,3)+o*4,SC(:,4)+o*5,SC(:,5)+o*6];
end


F=figure(2);
plot(zeit*1000,SC,'r','LineWidth',2)
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
xlabel('Time [ns]');
ylabel('Signal intensity')
str2=num2str(sigma_y);
title(['Çm RNA1-',str,' (',str2,'^o)']);
% title(['Çm RNA1-',str,' ']);

switch (str)
 case '8'
axis([0 1008 0.98 1.11]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '9'
axis([0 1200 0.98 1.13]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1 1.15]);
 case '10'
axis([0 1600 0.96 1.13]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '11'
axis([0 1200 0.95 1.11]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '12'
axis([0 1600 0.925 1.11]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '13'
axis([0 1584 0.935 1.125]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '14'
axis([0 1968 0.965 1.112]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
 case '15'
axis([0 1992 0.985 1.122]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000],'YTick',...
    [1 1.05 1.1]);
end
set(gca,'linewidth',1.5) 
% 
% 
% % % % %CmRNAModelB
% savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CmARNA_falscheRichtung\CmRNA_ModelB_Gband_y=0\Gband_CmARNA_ModelB1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CmARNA_falscheRichtung\CmRNA_ModelB_Gband_y=0\',['Gband_CmARNA_ModelB1_',str,'.pdf'])
% % % % 
% % % %CmRNAmodelA
% % savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CmRNA_ModelB_Gband_y=0\Gband_CmRNA_ModelB1_',str,'.fig'])
% % savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CmRNA_ModelB_Gband_y=0\',['Gband_CmRNA_ModelB1_',str,'.pdf'])
% % % 
% % % % 
% clear;
% clf;
% end 
