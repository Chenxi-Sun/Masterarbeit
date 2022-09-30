% % %dsDNA1_5-1_14
% % for num=5:14
% % str=num2str(num);
% % importdata(['C_DNA_1_',str,'_distr.dat']);
% % Y=ans(:,1);
% % f=ans(:,2);
% % [max_f,ind]=max(f);
% % N=length(Y);
% % HM=0.5*max(f);
% % for i=(Y(ind)-1):0.00001:Y(ind)
% %     vq=interp1(Y,f,i);
% %     if abs(vq-HM)<0.00001
% %         break
% %     end
% % end 
% % y1=i;
% % 
% % for i=Y(ind):0.0001:(Y(ind)+1)
% %     vq=interp1(Y,f,i);
% %     if abs(vq-HM)<0.0001
% %         break
% %     end
% % end 
% % y2=i;
% % FWHM(num-4)=y2-y1;
% % end 
load('D:\ChSun\Masterarbeit\AMmodel_RNA\allRNAdata.mat')
%%dsRNA1_8-1_15
for num=8:15
str=num2str(num);

switch (str)
 case '8'
zeit = real(RNA.RNA1_8.Qband.time);
Experimental.Sexp = RNA.RNA1_8.Xband.Sexp;
nr=8;
t_max=1100;
Y=RNA.RNA1_8.Qband.distr(:,1); 
F=RNA.RNA1_8.Qband.distr(:,2); 
 case '9'
zeit = real(RNA.RNA1_9.Qband.time);
Experimental.Sexp = RNA.RNA1_9.Xband.Sexp;
nr=9;
t_max=1100;
Y=RNA.RNA1_9.Qband.distr(:,1); 
F=RNA.RNA1_9.Qband.distr(:,2); 
 case '10'
zeit = real(RNA.RNA1_10.Qband.time);
Experimental.Sexp = RNA.RNA1_10.Xband.Sexp;
nr=10;
t_max=1200;
Y=RNA.RNA1_10.Qband.distr(:,1); 
F=RNA.RNA1_10.Qband.distr(:,2); 
 case '11'
zeit = real(RNA.RNA1_11.Qband.time);
Experimental.Sexp = RNA.RNA1_11.Xband.Sexp;
nr=11;
t_max=1600;
Y=RNA.RNA1_11.Qband.distr(:,1); 
F=RNA.RNA1_11.Qband.distr(:,2); 
 case '12'
zeit = real(RNA.RNA1_12.Qband.time);
Experimental.Sexp = RNA.RNA1_12.Xband.Sexp;
nr=12;
t_max=2000;
Y=RNA.RNA1_12.Qband.distr(:,1); 
F=RNA.RNA1_12.Qband.distr(:,2); 
 case '13'
zeit = real(RNA.RNA1_13.Qband.time);
Experimental.Sexp = RNA.RNA1_13.Xband.Sexp;
nr=13;
t_max=2200;
Y=RNA.RNA1_13.Qband.distr(:,1); 
F=RNA.RNA1_13.Qband.distr(:,2); 
 case '14'
zeit = real(RNA.RNA1_14.Qband.time);
Experimental.Sexp = RNA.RNA1_14.Xband.Sexp;
nr=14;
t_max=2200;
Y=RNA.RNA1_14.Qband.distr(:,1); 
F=RNA.RNA1_14.Qband.distr(:,2); 
 case '15'
zeit = real(RNA.RNA1_15.Qband.time);
Experimental.Sexp = RNA.RNA1_15.Xband.Sexp;
nr=15;
t_max=2500;
Y=RNA.RNA1_15.Qband.distr(:,1); 
F=RNA.RNA1_15.Qband.distr(:,2); 
end


Y=Y;
f=F;
[max_f,ind]=max(f);
N=length(Y);
HM=0.5*max(f);
for i=(Y(ind)-1):0.00001:Y(ind)
    vq=interp1(Y,f,i);
    if abs(vq-HM)<0.00001
        break
    end
end 
y1=i;

for i=Y(ind):0.0001:(Y(ind)+1)
    vq=interp1(Y,f,i);
    if abs(vq-HM)<0.0001
        break
    end
end 
y2=i;
FWHM(num-7)=y2-y1;
end 







