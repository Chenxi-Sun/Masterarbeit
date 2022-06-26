load('Z:\Students\ChSun\Masterarbeit\AMmodel_DNA\dsDNA_Cspin.mat');

zeit = real(DNApeldor.S0115.T(:,1));

summe= sum(DNApeldor.S0115.Sexp,2);

dlmwrite('C_DNA_1_15.DTA',[zeit, summe, zeros(length(summe),1)])