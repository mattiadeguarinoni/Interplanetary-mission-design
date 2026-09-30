function [deltat] = TOF(a, e, th1, th2, mu)
% legge oraria: noti il semiasse maggiore, l'eccentricità, l'anomalia vera
% iniziale, l'anomalia vera finale e la costante gravitazionale possiamo
% ricavare l'intervallo di tempo necessario per andare da th1 a th2
% NB: questo algoritmo è valido soltanto per orbite ellittiche

E1 = mod(2*atan(sqrt((1-e)/(1+e))*tan(th1/2)), 2*pi);
E2 = mod(2*atan(sqrt((1-e)/(1+e))*tan(th2/2)), 2*pi);

if E2 < E1
    E2 = E2 + 2*pi;
end

deltat = sqrt(a^3/mu)*((E2-E1)-e*(sin(E2)-sin(E1)));

end