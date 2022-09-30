
%%CdotDNA
Cdot_Q=importdata('CdotdsDNA_4PPeldor_fit.dat');
plot(Cdot_Q(:,1)*1000,Cdot_Q(:,2),'k','LineWidth',2);
hold on 
Cdot_X=importdata('sumX_Cdot_DNA1_9_fit.dat');
[devn,SC] = ScaleModdev('alle',Cdot_Q(1:157,2),Cdot_X(:,2));
diff=min(Cdot_Q(1:192,2))-min(SC)
plot(Cdot_X(:,1)*1000,SC+diff,'b','LineWidth',2);

%%CmdotRNA
% Cmdot_Q=importdata('Cmdot_sdRNA_4PPeldor_fit.dat');
% plot(Cmdot_Q(:,1).*1000,Cmdot_Q(:,2),'k','LineWidth',2);
% hold on 
% Cmdot_X=importdata('SumX_CmdotRNA1_12_fit.dat');
% [devn,SC] = ScaleModdev('alle',Cmdot_Q(1:156,2),Cmdot_X(:,2));
% diff=min(Cmdot_Q(1:156,2))-min(SC)
% plot(Cmdot_X(:,1).*1000,SC+diff,'b','LineWidth',2);


set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'linewidth',1.5) 

legend('Q-band','X-band')
xlabel('Time [ns]');
ylabel('Normalised intensity')