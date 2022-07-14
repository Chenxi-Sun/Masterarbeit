
% for nd=5:14
% clear;
Band = "Xband";
EP = importdata('ExperimentalParameters.mat');
dsDNA.Xband.EP.Sequence = EP.Sequence;
dsDNA.Xband.EP.Settings.PumpFrequency = EP.Settings.PumpFrequency;
dsDNA.Xband.EP.Settings.DetectionFrequency = EP.Settings.DetectionFrequency;
dsDNA.Xband.EP.Settings.B0 = EP.Settings.B0;
% 
% importdata('dsDNA_Cspin.mat');
load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\dsDNA_Cspin.mat');
% load('E:\Vorlesungen\EPR\Masterarbeit\ChSun\Masterarbeit\AMmodel_DNA\dsDNA_Cspin.mat');

nd='Which 2nd position? (5-14):';
str=input(nd,'s');

% str=num2str(nd);

t_max=2000;

switch (str)
 case '5'
zeit = real(DNApeldor.S0105.T)./1000;
Experimental.Sexp = DNApeldor.S0105.Sexp;
nr=5;
sigma_y=8;
 case '6'
zeit = real(DNApeldor.S0106.T)./1000;
Experimental.Sexp = DNApeldor.S0106.Sexp;
Experimental.Sexp=Experimental.Sexp(1:185,:);
nr=6;
sigma_y=18;
 case '7'
zeit = real(DNApeldor.S0107.T)./1000;
Experimental.Sexp = DNApeldor.S0107.Sexp;
nr=7;
sigma_y=6;
 case '8'
zeit = real(DNApeldor.S0108.T)./1000;
Experimental.Sexp = DNApeldor.S0108.Sexp;
nr=8;
sigma_y=0;
 case '9'
zeit = real(DNApeldor.S0109.T)./1000;
Experimental.Sexp = DNApeldor.S0109.Sexp;
nr=9;
sigma_y=0;
 case '10'
zeit = real(DNApeldor.S0110.T)./1000;
Experimental.Sexp = DNApeldor.S0110.Sexp;
nr=10;
sigma_y=0;
 case '11'
zeit = real(DNApeldor.S0111.T)./1000;
Experimental.Sexp = DNApeldor.S0111.Sexp;
nr=11;
sigma_y=0;
t_max=1700;
 case '12'
zeit = real(DNApeldor.S0112.T)./1000;
Experimental.Sexp = DNApeldor.S0112.Sexp;
nr=12;
sigma_y=4;
t_max=2000;
 case '13'
zeit = real(DNApeldor.S0113.T)./1000;
Experimental.Sexp = DNApeldor.S0113.Sexp;
nr=13;
sigma_y=10;
t_max=2000;
 case '14'
zeit = real(DNApeldor.S0114.T)./1000;
Experimental.Sexp = DNApeldor.S0114.Sexp;
nr=14;
sigma_y=12;
t_max=3000;

end

sigma_y=6;

%%Model A/B
% 
% for i=1:100
% alpha=76.2/360*2*pi;               %first rotation angle alpha
% beta(i)=beta_range(i)/360*2*pi;                %second rotation angle beta
% % [Simulated,R_mean,sigma] = AMPELDOR_DNA(nr,EP,zeit,'A');
% sigma_y=0;
% [Simulated,R_mean,FWHM,ymax] = AMPELDOR_DNA(sigma_y,nr,EP,zeit,'A');
[Simulated,R_mean,FWHM] = AMPELDOR_DNA(sigma_y,nr,EP,zeit,'B');
% R_mean_all(nd-4,:)=R_mean;
% FWHM_all(nd-4,:)=FWHM;
% end

%%Model C
% [Simulated] = AM_C_PELDOR_DNA(nr,EP,zeit,'A');

% % Model A/B change sigma_y
% N=length(Experimental.Sexp(:,1));
% o=0.1;
% for i=1:11
% sigma_y=0+(i-1)*2;
% Simulated = AMPELDOR_DNA(sigma_y,nr,EP,zeit,'A');
% [devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
% Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
% SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];
% sigma(i,:)=sqrt(sum((SC-Experimental.stack).^2,1)/N);
% end 
% RSMD_alloffset=sum(sigma,2)./6;
% % RSMD_alloffset=sigma(:,1);
% [m,n]=min(RSMD_alloffset);
% sigma_y_best=(n-1)*2;
% [Simulated,R_mean,sigma] = AMPELDOR_DNA(sigma_y_best,nr,EP,zeit,'A');
% % % 


% PLOT
[devn,SC] = ScaleModdev('alle',Experimental.Sexp,Simulated.Sexp);
o = 0.1; 
Experimental.stack = [Experimental.Sexp(:,1),Experimental.Sexp(:,2)+o,Experimental.Sexp(:,3)+2*o,Experimental.Sexp(:,4)+3*o,Experimental.Sexp(:,5)+4*o,Experimental.Sexp(:,6)+5*o];
SC = [SC(:,1),SC(:,2)+o,SC(:,3)+o*2,SC(:,4)+o*3,SC(:,5)+o*4,SC(:,6)+o*5];

F=figure(2);

% plot(zeit*1000,SC,'r','LineWidth',2)
plot(zeit*1000,SC,'b','LineWidth',2)
hold on
plot(zeit*1000,Experimental.stack,'k','LineWidth',2)
xlabel('Time [ns]');
ylabel('Signal intensity')
str2=num2str(sigma_y);
title(['Ç DNA1-',str,' (',str2,'^o)']);
axis([0 t_max 0.3 1.5]);
% axis([0 2000 0.3 1.5]);
set(gca,'FontSize',14,'FontWeight','bold','XTick',...
    [0 500 1000 1500 2000 2500 3000]);
set(gca,'linewidth',1.5) 

% saveas(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\helixdy_SL\ModelA1_',str,'.png'])


% savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\CDNAModelB1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelB\',['CDNAModelB1_',str,'.pdf'])

% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\helixdy_SL\ModelB1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\helixdy_SL\',['ModelB1_',str,'.pdf'])

% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\only_helixdynamics\ModelB1_',str2,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelB\only_helixdynamics\',['ModelB1_',str2,'.pdf'])
% 
% savefig(F,['Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\only_dynamics\ModelA1_',str2,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\AMmodel_result\AMmodelA\only_dynamics\',['ModelA1_',str2,'.pdf'])
% clear;
% clf;
% end 


% importdata(['C_DNA_1_',str,'_distr.dat']);
% hold on 
% plot(ans(:,1),ymax/max(ans(:,2)).*ans(:,2),'k','Linewidth',2)
% xlim([1 5])
% xlabel('distance [nm]')
% ylabel('no. of occurences')
% xticks([1 2 3 4 5 6])
% xlim([1 5])
% 
% 
% F=figure(1);
% set(gca,'FontSize',14,'FontWeight','bold','XTick',...
%     [1 2 3 4 5 6]);
% set(gca,'linewidth',1.5) 
% xlabel('Distance [nm]')
% ylabel('no. of occurences')
% if str2num(str)<12
%     xlim([1 5])
% else 
%     xlim([2 6])
% end 
% str2=num2str(sigma_y);
% title(['Ç DNA1-',str,' (',str2,'^o)']);
% % title(['C DNA1-',str,' (',str2,'^o)']);
% set(gca,'linewidth',1.5) 
% 
% % 
% savefig(F,['Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelA\CDNA_ModelA1_',str,'.fig'])
% savePDF('Z:\Students\ChSun\Masterarbeit\11.07_Result\CDNA_ModelA\',['CDNA_ModelA1_',str,'.pdf'])
% % 
% clear;
% clf;
% end 

