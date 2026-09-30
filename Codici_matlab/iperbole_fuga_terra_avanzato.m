%% iperbole di fuga dalla Terra avanzato

clear
clc

v_inf_eclip = [1.4130 2.0178 -2.0936]'; % trovata da scenario 2
a_p = 2.2119*10^(4);
e_p = 0.5897;
minimo = inf;
mu_t = 398600.4;
R_terra = 6371;

options = optimoptions('fsolve', ...
    'Display', 'none', ...
    'TolFun', 1e-10, ...
    'TolX', 1e-10, ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 1000);

for theta = 0:0.1:2*pi

    p = a_p*(1 - e_p^2);
    r_h = p/(1 + e_p*cos(theta));
    eps = deg2rad(23.45);
    T_eci_eclip = [1 0 0; 0 cos(eps) sin(eps); 0 -sin(eps) cos(eps)];
    v_inf_eci = T_eci_eclip'*v_inf_eclip;
    
    a_h = -mu_t/(norm(v_inf_eci)^2);

    [r_h_vec,v_op] = par2car(a_p,e_p,0.5560,1.9093,1.8783,theta,mu_t);
    r_h_vers = r_h_vec/norm(r_h_vec);
    r_inf = v_inf_eci/norm(v_inf_eci);

    % angolo alpha tra r_h_vers e r_inf
    alpha = acos(dot(r_h_vers, r_inf));
    
    % x = [e_h theta_h theta_inf] vettore contenente le incognite
    fun = @(x) [r_h-(a_h*(1-x(1)^2)/(1+x(1)*cos(x(2))));...
                cos(x(3))+1/x(1); ...
                x(3)-alpha-x(2)];

    x0 = [1.5, deg2rad(80), deg2rad(120)]'; % condizioni iniziali
    sol = fsolve(fun,x0,options);

    r_p_h = -a_h*(sol(1) - 1);

    if r_p_h > R_terra && sol(3) > pi/2 && sol(3) < pi
        % sol_v = [sol_v sol];
        h_h = cross(r_h_vec,r_inf)/norm(cross(r_h_vec,r_inf));
        i_h = acos(h_h(3));
        k = [0 0 1]';
        N = cross(k,h_h)/norm(cross(k,h_h));

        if N(2) >= 0
            OM_h = acos(N(1));
        else
            OM_h = 2*pi - acos(N(1));
        end

        cos_b = dot(N,r_h_vers);
        sin_b = dot(cross(N,r_h_vers),h_h);
        beta = atan2(sin_b,cos_b);
        om_h = beta - sol(2);
       
        [r_ip,v_ip] = par2car(a_h,sol(1),i_h,OM_h,om_h,sol(2),mu_t);
        
        Delta_v = norm(v_ip - v_op);
        
        if Delta_v < minimo
            minimo = Delta_v;
            e_h_best = sol(1);
            i_h_best = i_h;
            om_h_best = om_h;
            OM_h_best = OM_h;
            theta_best = sol(2);
        end
    end

end

% tempo di trasferimento
m_terra = 5.972*10^(24);
m_sole = 1.989*10^(30);
a_terra = 1.4946 * 10^(8);
r_SOI = a_terra*((m_terra/m_sole)^(2/5));
p = abs(a_h)*(e_h_best^2 - 1);
theta1 = theta_best;
theta2 = acos((p/r_SOI -1)/e_h_best);

F1 = 2*atanh(sqrt((e_h_best-1)/(e_h_best+1))*tan(theta1/2));
F2 = 2*atanh(sqrt((e_h_best-1)/(e_h_best+1))*tan(theta2/2));

delta_t = sqrt(abs(a_h)^3/mu_t)*((e_h_best*sinh(F2)-F2)-(e_h_best*sinh(F1)-F1));