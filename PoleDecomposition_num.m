% Base function to calculate numerically the response of each pole of the BOLD signal
% k and w can either be single values or matrix values
% Revised May 8, 2015 -> made structure arrays
% Revised May 15, 2015 -> used params class params.m
% James Pang

function [w_pole, a, T, kz] = PoleDecomposition_num(p, v_b, Gamma, k, w)
    
    % additional parameters
    D  = p.rho_f*(2*Gamma - p.beta*p.Cz/p.tau);    
    kz = sqrt((p.k_0)^2 + 1/(v_b)^2*p.Cz*(p.beta/p.tau)*(D/p.rho_f));


    P = -p.Cz*(p.k2 - p.k3 - p.V_0*(p.k1 + p.k2));
    Q = p.Cz*((p.k2 - p.k3)*(p.eta + p.tau^(-1)) - ...
        (p.k1 + p.k2)*p.Cz*(p.eta - (p.tau^(-1))*(p.beta - 2)) + ...
        (D/p.rho_f)*(p.k2 - p.k3 - p.V_0*(p.k1 + p.k2)));
    R = p.Cz*(D/p.rho_f)*((p.k2 - p.k3)*(p.eta + p.tau^(-1)) - ...
        (p.k1 + p.k2)*p.Cz*(p.eta - (p.tau^(-1))*(p.beta - 2)));


    w1 = -1i*Gamma - sqrt(k.^2*v_b^2 + kz^2*v_b^2 - Gamma^2);
    w2 = -1i*Gamma + sqrt(k.^2*v_b^2 + kz^2*v_b^2 - Gamma^2);
    w3 = (-0.5*1i*p.kappa - p.w_f)*ones(size(k));
    w4 = (-0.5*1i*p.kappa + p.w_f)*ones(size(k));
    w5 = (-1i*p.eta - 1i*(p.tau^(-1)))*ones(size(k));

    a1 = (w1.^2*P*1i + w1*Q + R*1i)./((w1-w2).*(w1-w3).*(w1-w4).*(w1-w5)+eps);
    a2 = (w2.^2*P*1i + w2*Q + R*1i)./((w2-w1).*(w2-w3).*(w2-w4).*(w2-w5)+eps);
    a3 = (w3.^2*P*1i + w3*Q + R*1i)./((w3-w1).*(w3-w2).*(w3-w4).*(w3-w5)+eps);
    a4 = (w4.^2*P*1i + w4*Q + R*1i)./((w4-w1).*(w4-w2).*(w4-w3).*(w4-w5)+eps);
    a5 = (w5.^2*P*1i + w5*Q + R*1i)./((w5-w1).*(w5-w2).*(w5-w3).*(w5-w4)+eps);

    T1 = a1./(w - w1);
    T2 = a2./(w - w2);
    T3 = a3./(w - w3);
    T4 = a4./(w - w4);
    T5 = a5./(w - w5);
    Ttotal = T1 + T2 + T3 + T4 + T5;

    w_pole = struct('w1',w1,'w2',w2,'w3',w3,'w4',w4,'w5',w5);
    a = struct('a1',a1,'a2',a2,'a3',a3,'a4',a4,'a5',a5);
    T = struct('T1',T1,'T2',T2,'T3',T3,'T4',T4,'T5',T5,'Ttotal',Ttotal);

end
    
    