clear
DFT = importdata('Cspin_DFT_Edwards2011_import_noH.xyz');
turn_file = string(zeros(100,4));
turn_file(1,1)=string(14);
turn_file(1,2:4)=missing;
turn_file(2,1)="SpinlabelTurn";
turn_file(2,2:4)=missing;
turn_file(3:3+13,1)=['C';'C';'C';'C';'O';'N';'C';'C';'C';'C';'C';'C';'C';'C'];

turn_file(3:3+13,2:4)=DFT.data;
xyz=DFT.data;
anfang=coord_smallSL(xyz);

tnrangl = [0:10:90];
alle=[];
for d=1:length(tnrangl)
   alle=[alle;repmat(tnrangl(d),30,1)]; 
end

% for mdl=1:90
for mdl=1:length(alle)
    turn_file(16*mdl+1,1)=string(14);    turn_file(16*mdl+1,2:4)=missing;
    turn_file(16*mdl+2,1)="SpinlabelTurn";    turn_file(16*mdl+2,2:4)=missing;
    
    for d=1:length(xyz)
%         str90(d,:)=drehenRodrigues(anfang.xyz(d,:)-anfang.M,anfang.xaxis1,mdl)+anfang.M;
            str90(d,:)=drehenRodrigues(anfang.xyz(d,:)-anfang.M,-anfang.zaxis1,alle(mdl))+anfang.M;
%             str90(d,:)=drehenRodrigues(anfang.xyz(d,:)-anfang.M,anfang.xaxis1,alle(mdl))+anfang.M;
    end
    
    str90coords = coord_smallSL(str90);
    str90 = str90+ (anfang.M-str90coords.M);
    
    turn_file(16*mdl+3:16*mdl+16,2:4)=str90;
    turn_file(16*mdl+3:16*mdl+16,1)=['C';'C';'C';'C';'O';'N';'C';'C';'C';'C';'C';'C';'C';'C'];

end
    

writematrix(turn_file, 'alphaturn.txt')

% file_ID = fopen('testfile.txt','w');
% fprintf(file_ID,'%s %s %s\n',turn_file);
% 
% fclose(file_ID);
%%
function all=coord_smallSL(xyz)
all.N  = xyz(10,:);
all.C2 = xyz(1,:);
all.C4 = xyz(3,:);
all.C5 = xyz(8,:);

all.NNO = xyz(6,:);
all.ONO = xyz(5,:);


NObond= all.ONO-all.NNO;
NC5   = all.C5-all.N;
NC2   = all.C2-all.N;

for ra=1:1
    all.xaxis1(ra,:) = NObond(ra,:)/norm(NObond(ra,:));
    NC2(ra,:) = NC2(ra,:)/norm(NC2(ra,:));
    all.zaxis1(ra,:) = cross(all.xaxis1(ra,:),NC2(ra,:));
    all.zaxis1(ra,:) = all.zaxis1(ra,:)/norm(all.zaxis1(ra,:));
    all.yaxis1(ra,:) = cross(all.zaxis1(ra,:),all.xaxis1(ra,:));
end
% M = C5 + 0.5*(C4-C5) + xaxis1.*(sqrt(7^2-0.5*vecnorm((C4-C5)')'.^2)+sqrt(8.1^2-0.5*vecnorm((C4-C5)')'.^2))/2;
% all.M = all.C5 + 0.5*(all.C4-all.C5) + all.xaxis1*(sqrt(7^2-0.75^2)+sqrt(8.2^2-0.75^2))/2;
all.M =0.5*(all.NNO+all.ONO);
all.xyz=xyz;
end


