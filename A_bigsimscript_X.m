clear all

FF = "DESREF"; %OL3, DESREF, BSC1 (even, odd)
% FF = "OL3";
% FF = "BSC1";
Band = "Xband";

if FF == "OL3"
%     Zahl = [8:15]; def = ["even","odd","even","odd","even","odd","even","odd"];m=2; n=4;
    Zahl = [13]; def = ["odd"]; m=1; n=1;
    
elseif FF == "DESREF"
%     Zahl = [8,10,12,14]; def = ["even","even","even","even"];m=2; n=2;
%     Zahl = [8:15]; def = ["even","odd","even","odd","even","odd","even","odd"];m=2; n=4;
    Zahl = [9,11,13,15]; def = ["odd","odd","odd","odd"]; m=2; n=2;
%     Zahl = [8]; def = ["even"];m=1; n=1;
%     Zahl = [8,10]; def = ["even","even"];m=2; n=2;
    
elseif FF == "BSC1"
    Zahl = [9,11,13,15]; def = ["Odd","Odd","Odd","Odd"]; m=2; n=2;    
end

figure

for k = 1:length(Zahl)
    load(sprintf('RNA1_%d_new.mat',Zahl(1,k)));
    probe = sprintf('RNA1_%d',Zahl(1,k));
    probe = eval(probe);
%     
%     load('RNA1_8_rect.mat') 
%     forrect.Xband=forrect.Xbandrect;
%     probe=forrect;
% %     
%     probe.Xband.EP.Sequence.FirstPulse.Length=46;
%     probe.Xband.EP.Sequence.SecondPulse.Length=72;
%     probe.Xband.EP.Sequence.ThirdPulse.Length=72;
%     probe.Xband.EP.Sequence.PumpPulse.Length=72;

%     lang=100;
%     probe.Xband.EP.Sequence.FirstPulse.Length=lang;
%     probe.Xband.EP.Sequence.SecondPulse.Length=lang;
%     probe.Xband.EP.Sequence.ThirdPulse.Length=lang;
%     probe.Xband.EP.Sequence.PumpPulse.Length=lang;
     
    
    MD1 = sprintf('BP1_RNA_%s_%s.txt',FF,def(1,k));
    MD2 = sprintf('BP%d_RNA_%s_%s.txt',Zahl(1,k),FF,def(1,k));
    Param = sprintf('probe.%s.EP',Band);
    zeit = sprintf('probe.%s.time',Band);
    name = sprintf("RNA1_%d",Zahl(1,k));
    
    RNA.(name) = autoneurechner(MD1,MD2,Param,zeit,probe);
    SIM.(Band).(FF) = RNA.(name);
%     save(sprintf('RNA1_%d_new.mat',Zahl(1,k)),'SIM','-append')
    
    %changes Mod depth, stacks graphs
    [devn,SC] = ScaleModdev('alle',probe.(Band).Sexp,RNA.(name).Sexp);

    o = 0;
    for j=1:size(probe.(Band).Sexp,2)
        probe.(Band).stack(:,j) = probe.(Band).Sexp(:,j) + o;
        RNA.(name).SC(:,j) = SC(:,j) + o;
        o = o+0.2 ;
    end 
        
    subplot(m,n,k)
    
    plot(probe.(Band).time,RNA.(name).SC,'r',probe.(Band).time,probe.(Band).stack,'k')
    title(sprintf("RNA 1-%d",Zahl(1,k)))
    xlabel('time [\mus]');
    ylim([min(probe.(Band).stack(:,1))-0.05 max(probe.(Band).stack(:,size(probe.(Band).stack,2)))+0.05]);
    xlim([0 max(probe.(Band).time)]);
    xticks([0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0]);
    hold on
end

