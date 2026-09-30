function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)
% trasferimento biellittico bitangente: imponiamo la geometria dell'orbita
% iniziale, la geometria dell'orbita finale, il raggio di apocentro delle
% due orbite di trasferimento e la costante gravitazionale --> troviamo il
% costo della manovra in termini di velocità (dobbiamo fare tre manovre,
% quindi avremo tre DeltaV) e il tempo di manovra (abbiamo due orbite di
% trasferimento, quindi avremo due Deltat)

rp_i = ai*(1 - ei);
rp_f = af*(1 - ef);

rp_t1 = ai*(1 - ei);
rp_t2 = af*(1 - ef);

at1 = (rp_t1 + ra_t)/2;
at2 = (rp_t2 + ra_t)/2;

DeltaV1 = sqrt(mu)*sqrt((2/rp_t1)-(1/at1)) - sqrt(mu)*sqrt((2/rp_i)-(1/ai));
DeltaV2 = sqrt(mu)*sqrt((2/ra_t)-(1/at2)) - sqrt(mu)*sqrt((2/ra_t)-(1/at1));
DeltaV3 = sqrt(mu)*sqrt((2/rp_f)-(1/af)) - sqrt(mu)*sqrt((2/rp_t2)-(1/at2));

Deltat1 = pi*sqrt((at1^3)/mu);
Deltat2 = pi*sqrt((at2^3)/mu);

end