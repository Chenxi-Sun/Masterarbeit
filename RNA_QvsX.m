% for nr=8:15
% 
clear;
% load('Z:\Students\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
load('D:\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
nd='Which 2nd position? (8-15):';
str=input(nd,'s');

    
% str=num2str(nr);

switch (str)
 case '8'
zeit = real(RNA.RNA1_8.Qband.time);
Experimental.Sexp = RNA.RNA1_8.Xband.Sexp;
nr=8;
t_max=1100;
Exp(:,1)=RNA.RNA1_8.Qband.time; 
Exp(:,2)=RNA.RNA1_8.Qband.Sexp; 
[n_1,m_1]=min(Exp(:,2));
X_Exp(:,1)=RNA.RNA1_8.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_8.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2),'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:78,2));
plot(X_Exp(1:78,1).*1000,SC,'b','Linewidth',2)
 case '9'
zeit = real(RNA.RNA1_9.Qband.time);
Experimental.Sexp = RNA.RNA1_9.Xband.Sexp;
nr=9;
t_max=1100;
Exp(:,1)=RNA.RNA1_9.Qband.time; 
Exp(:,2)=RNA.RNA1_9.Qband.Sexp; 
X_Exp(:,1)=RNA.RNA1_9.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_9.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2)+0.2,'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:78,2));
plot(X_Exp(1:78,1).*1000,SC+0.2,'b','Linewidth',2)
 case '10'
zeit = real(RNA.RNA1_10.Qband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp;
nr=10;
t_max=1200;
Exp(:,1)=RNA.RNA1_10.Qband.time; 
Exp(:,2)=RNA.RNA1_10.Qband.Sexp; 

X_Exp(:,1)=RNA.RNA1_10.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_10.Xband.Sexp(:,1:6),2)/6; 

plot(Exp(:,1).*1000,Exp(:,2),'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:96,2));
plot(X_Exp(1:96,1).*1000,SC,'b','Linewidth',2)
 case '11'
zeit = real(RNA.RNA1_11.Qband.time);
Experimental.Sexp = RNA.RNA1_11.Xband.Sexp;
nr=11;
t_max=1600;
Exp(:,1)=RNA.RNA1_11.Qband.time; 
Exp(:,2)=RNA.RNA1_11.Qband.Sexp; 
X_Exp(:,1)=RNA.RNA1_11.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_11.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2)+0.2,'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:96,2));
plot(X_Exp(1:96,1).*1000,SC+0.2,'b','Linewidth',2)
 case '12'
zeit = real(RNA.RNA1_12.Qband.time);
Experimental.Sexp = RNA.RNA1_12.Xband.Sexp;
nr=12;
t_max=2000;
Exp(:,1)=RNA.RNA1_12.Qband.time; 
Exp(:,2)=RNA.RNA1_12.Qband.Sexp; 

X_Exp(:,1)=RNA.RNA1_12.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_12.Xband.Sexp(:,1:6),2)/6; 

plot(Exp(:,1).*1000,Exp(:,2),'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:96,2));
plot(X_Exp(1:96,1).*1000,SC,'b','Linewidth',2)
 case '13'
zeit = real(RNA.RNA1_13.Qband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp;
nr=13;
t_max=2200;
Exp(:,1)=RNA.RNA1_13.Qband.time; 
Exp(:,2)=RNA.RNA1_13.Qband.Sexp; 

X_Exp(:,1)=RNA.RNA1_13.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_13.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2)+0.2,'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:122,2));
plot(X_Exp(1:122,1).*1000,SC+0.2,'b','Linewidth',2)
 case '14'
zeit = real(RNA.RNA1_14.Qband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp;
nr=14;
t_max=2200;
Exp(:,1)=RNA.RNA1_14.Qband.time; 
Exp(:,2)=RNA.RNA1_14.Qband.Sexp; 

X_Exp(:,1)=RNA.RNA1_14.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_14.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2),'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:146,2));
plot(X_Exp(1:146,1).*1000,SC,'b','Linewidth',2)
 case '15'
zeit = real(RNA.RNA1_15.Qband.time);
Experimental.Sexp = RNA.RNA1_15.Xband.Sexp;
nr=15;
t_max=2500;
Exp(:,1)=RNA.RNA1_15.Qband.time; 
Exp(:,2)=RNA.RNA1_15.Qband.Sexp; 

X_Exp(:,1)=RNA.RNA1_15.Xband.time; 
X_Exp(:,2)=sum(RNA.RNA1_15.Xband.Sexp(:,1:6),2)/6; 
plot(Exp(:,1).*1000,Exp(:,2)+0.2,'k','Linewidth',2)
hold on 
[devn,SC] = ScaleModdev('alle',Exp(:,2),X_Exp(1:146,2));
plot(X_Exp(1:146,1).*1000,SC+0.2,'b','Linewidth',2)
end
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'linewidth',1.5) 
xlim([0 2000])
legend('Q-band','X-band')
xlabel('Time [ns]');
ylabel('Normalised intensity')
% end 