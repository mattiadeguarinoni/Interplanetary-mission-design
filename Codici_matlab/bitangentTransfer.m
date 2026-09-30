function [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)
% trasferimento bitangente: imponiamo la geometria dell'orbita iniziale, la
% geometria dell'orbita finale, il tipo di manovra che vogliamo effettuare
% e la costante gravitazionale --> troviamo il costo della manovra in
% termini di velocità (dobbiamo fare due manovre, quindi avremo due DeltaV)
% e il tempo di manovra

rp_i = a_i*(1 - e_i);
ra_i = a_i*(1 + e_i);
rp_f = a_f*(1 - e_f);
ra_f = a_f*(1 + e_f);

switch type

    case 'pa'
        rp_t = rp_i;
        ra_t = ra_f;
        a_t = (rp_t + ra_t)/2;
        DeltaV1 = sqrt(mu)*sqrt((2/rp_t)-(1/a_t)) - sqrt(mu)*sqrt((2/rp_i)-(1/a_i));
        DeltaV2 = sqrt(mu)*sqrt((2/ra_f)-(1/a_f)) - sqrt(mu)*sqrt((2/ra_t)-(1/a_t));
        Deltat = pi*sqrt((a_t^3)/mu);

    case 'ap'
        rp_t = rp_f;
        ra_t = ra_i;
        a_t = (rp_t + ra_t)/2;
        DeltaV1 = sqrt(mu)*sqrt((2/ra_t)-(1/a_t)) - sqrt(mu)*sqrt((2/ra_i)-(1/a_i));
        DeltaV2 = sqrt(mu)*sqrt((2/rp_f)-(1/a_f)) - sqrt(mu)*sqrt((2/rp_t)-(1/a_t));
        Deltat = pi*sqrt((a_t^3)/mu);

    case 'pp'
        rp_t = rp_i;
        ra_t = rp_f;
        a_t = (rp_t + ra_t)/2;
        DeltaV1 = sqrt(mu)*sqrt((2/rp_t)-(1/a_t)) - sqrt(mu)*sqrt((2/rp_i)-(1/a_i));
        DeltaV2 = sqrt(mu)*sqrt((2/rp_f)-(1/a_f)) - sqrt(mu)*sqrt((2/ra_t)-(1/a_t));
        Deltat = pi*sqrt((a_t^3)/mu);

    case 'aa'
        rp_t = ra_i;
        ra_t = ra_f;
        a_t = (rp_t + ra_t)/2;
        DeltaV1 = sqrt(mu)*sqrt((2/rp_t)-(1/a_t)) - sqrt(mu)*sqrt((2/ra_i)-(1/a_i));
        DeltaV2 = sqrt(mu)*sqrt((2/ra_f)-(1/a_f)) - sqrt(mu)*sqrt((2/ra_t)-(1/a_t));
        Deltat = pi*sqrt((a_t^3)/mu);

end

end