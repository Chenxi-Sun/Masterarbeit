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

%%dsRNA1_8-1_15
for num=8:15
str=num2str(num);
importdata(['SumX_CmRNA1_',str,'_distr.dat']);
Y=ans(:,1);
f=ans(:,2);
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







