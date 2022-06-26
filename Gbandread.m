% written by Grytz

% sum up single time traces measured for example at G-band and check for
% artifacts

clear all
figure

% datei='PELDOR_xx_0G';
% a=39;b=183;

% datei='PELDOR-64xy';
% a=50;b=87;

datei='PELDOR_yy_51G';
a=1;b=22;
% % a=149;b=249;

% datei='PELDOR+5yz';
% a=98;b=142;

% datei='PELDOR_zz_195G';
% a=184;b=252;

dumm = 'add'; %choose oder add


if a<=9
    erste=['00',num2str(a)];
elseif a<=99
    erste=['0',num2str(a)];
else
    erste=[num2str(a)];
end


    
[t,re1,im1]=textread([erste,datei,'.dat'],'%f%f%f');
if abs(sum(re1))>=abs(sum(im1))
if sum(re1)>=0
    re=re1;
    im=im1;
else
    re=-re1;
    im=-im1;
end
else
    if sum(im1)>=0
        re=im1;
        im=re1;
    else
        re=-im1;
        im=-re1;
    end
end


good.ren=zeros(1,length(re))'; %length of real signal and imagin???r signal vector
good.imn=zeros(1,length(re))';
bad.ren=zeros(1,length(re))'; %for bad timetraces
bad.imn=zeros(1,length(re))';
all.ren=zeros(1,length(re))'; %for all timetraces
all.imn=zeros(1,length(re))';

good.list=[];
bad.list=[];
a1=a+1;
for k=a1:b%loop over all files in the folder
    if k<=9
        name=sprintf(['00',num2str(k),datei,'.dat']);  %experimental files should have an index
    elseif k<=99
        name=sprintf(['0',num2str(k),datei,'.dat']) ; %experimental files should have an index
    else
        name=sprintf([num2str(k),datei,'.dat']) ; %experimental files should have an index
    end
    
    [t,re1,im1]=textread(name,'%f%f%f');
    if abs(sum(re1))>=abs(sum(im1))
        if sum(re1)>=0
            re=re1;
            im=im1;
        else
            re=-re1;
            im=-im1;
        end
    else
        if sum(im1)>=0
            re=im1;
            im=re1;
        else
            re=-im1;
            im=-re1;
        end
    end
    all.ren=all.ren+re;all.imn=all.imn+im; %add all timetraces
    
    subplot(2,1,1)
    plot(t,re,'b')
    xlabel('time'); ylabel('signal')
    subplot(2,1,2)
    plot(t,im,'r')
    % plot sum timetrace
    xlabel('time'); ylabel('signal')
    title('sum timetrace selected')
    
    switch dumm
        case 'choose'    
            x=input('Is timetrace ok, then typ 0 for yes and 1 for no:')
            if x==0%add selected timetraces
                good.ren=good.ren+re;good.imn=good.imn+im;
                good.list=cat(3,good.list,name);
            else
                bad.ren=bad.ren+re;bad.imn=bad.imn+im;%add not selected timetraces

                bad.list=cat(3,bad.list,name);
            end
            fprintf('%d done, ende %d\n',k,b)
        case 'add'
            good.ren=good.ren+re;good.imn=good.imn+im;
            good.list=cat(3,good.list,name);
    end
  
    
end


subplot(2,1,1)
plot(t,good.ren,'b') 
xlabel('time'); ylabel('signal')
subplot(2,1,2)
plot(t,good.imn,'r')
% plot sum timetrace
xlabel('time'); ylabel('signal')
title('sum timetrace selected')

% 
save ([datei,'.mat'],'good','bad','all','t')
dlmwrite([datei,'.dat'],[t*1000000,good.ren,good.imn])
