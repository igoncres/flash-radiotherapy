function [D50c, D50f, D90c, D90f] = getD50_D90(dose, TCPc, TCPf)

%INFO: function to calculate D50 and D90 from TCP curves

%INPUT
%dose [vector] doses in Gy
%TCPc: [matrix] TCP curves with CONV-RT (tumors x doses)
%TCPc: [matrix] TCP curves with FLASH-RT (tumors x doses)

%OUTPUT
%D50c, D50f, D90c, D90f: [real]

    TCP = mean(TCPc);

    ind=find(TCP>0.1 & TCP<0.9);
    pp = spline(TCP(ind),dose(ind));
    D50c = ppval(pp,0.5);

    ind=find(TCP>0.8 & TCP<0.99);
    pp = spline(TCP(ind),dose(ind));
    D90c = ppval(pp,0.9);

    TCP = mean(TCPf);

    ind=find(TCP>0.1 & TCP<0.9);
    pp = spline(TCP(ind),dose(ind));
    D50f = ppval(pp,0.5);

    ind=find(TCP>0.8 & TCP<0.99);
    pp = spline(TCP(ind),dose(ind));
    D90f = ppval(pp,0.9);

end