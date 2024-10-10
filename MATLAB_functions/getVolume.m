function [time, C, Cd, V] = getVolume(V0, SF, lambda, phi, K, tspan, RTtime, dt_vol)

%INFO: function to obtain volume curves

%INPUT
%V0: [real] Initial volume
%SF: [real] Surviving fraction (=1 for control)
%lambda: [real] Proliferation rate
%phi: [real] Elimination rate
%K: [real] Carrying capacity
%tspan: [vector] Initial and final time of simulation
%RTtime: [vector] Irradiation times
%dt_vol: [real] Time step

%OUTPUT
%time: [vector] Time points
%C: [vector] Non-damaged tumor cells
%Cd: [vector] Damaged tumor cells
%V: [vector] Volume

time = (0:dt_vol:tspan(end));
nt = length(time);
C = zeros(1,nt);
Cd = zeros(1,nt);

if RTtime(1)==0
    C(1) = SF*V0;
    Cd(1) = (1-SF)*V0;
end

for i = 1:nt-1
    % One RT fraction
    C(i+1) = C(i)+dt_vol*lambda*C(i)*(1-(C(i)+Cd(i))/K);
    Cd(i+1) = max(0, Cd(i)-dt_vol*(lambda*Cd(i)*(C(i)+Cd(i))/K+phi*Cd(i)));

    %Multiple RT fractions
    %RT: 1 yes; 0 no
%     RT = sum(abs(RTtime-time(i+1))<dt_vol/2);
%     if RT
%         Cd(i+1) = Cd(i+1) + (1-SF)*C(i+1);
%         C(i+1) = SF*C(i+1);
%     end

end

V = C + Cd;

end