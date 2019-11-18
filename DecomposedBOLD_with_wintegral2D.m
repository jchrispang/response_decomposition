%%
% This file calculates the 2D hemodynamic decomposedresponse function, 
% as a function of radius and time with the w integral left for ctfft
 
% Dependencies:
% 1. params.m from J.C. Pang
% 2. PoleDecomposition_num.m from J.C. Pang
% 3. Hankel transform algorithm based on work provided by
%    M. Guizar-Sicairos and J. C. Gutierrez-Vega, Computation of 
%    quasi-discrete Hankel transforms of integer order for propagating optical 
%    wave fields, J. Opt. Soc. Am. A 21, 53-58 (2004).
% 4. c.mat from mathworks upload of (3)
% 5. ctfft.m from K. Aquino

%%
function [Y, kz] = DecomposedBOLD_with_wintegral2D(p, v_b, Gamma, stimulus)

    % INPUT
    % 1. parameters: p
    % 2. wave properties: speed v_b, damping rate Gamma
    % 3. stimulus: stimulus in k, w space (Note that w is in x and k is in y)
    %    --> sample form [f, fk] = meshgrid(linspace(-freqMax, freqMax, Nw), ...
    %                    linspace(-spatialFreqMax, spatialFreqMax, Nk)); 
    %                    stimulus = f^2 + fk^2 
    %
    % OUTPUT
    % 1. structure Y containing BOLD response in r vs t of each pole contribution
    %
    % EXAMPLE: Y = DecomposedBOLD_with_wintegral(params, 1e-3, 1, stimulus)

    
    Nw = p.Nw;
    freqMax = p.freqMax;
    Nk = p.Nk;
    spatialFreqMax = p.spatialFreqMax;
    
    % predefine vectors for temporal frequency
    jvec = 0:Nw-1;
    f = (jvec - Nw/2)/Nw*freqMax*2;
    w = 2*pi*f;                                 %% w = \omega
    
    t = 1/(f(2)-f(1))*1/Nw*(jvec - Nw/2);       %% time vector


    % The Hankel transform: Algorithm by M. Guizar-Sicairos and J. C.
    % Gutierrez-Vega, using a fourier-bessel expansion to calculate the Hankel
    % transform.

    % Input parameters  

    kMax = spatialFreqMax*2*pi;     %% maximum spatial frequency k 
    N = Nk;                         %% Number of sampling points
    ord = 0;                        %% Transformation order

    % Matrix and vectors computing  
    % This operations may only be performed once for iterative algorithms

    load c.mat;
    c = c(ord+1, 1:N+1);

    rMax = c(N+1)/(2*pi*kMax);      %% maximum radius r
    k = c(1:N)'*kMax/c(N+1);        %% k vector
    r = c(1:N)'/(2*pi*kMax);        %% radius vector

    r = 2*pi*r;     % Hankel transform is with 2pi, ours is with 1/2pi

    [Jn,Jm] = meshgrid(c(1:N),c(1:N));
    C = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
    % C is the transformation matrix

    m1 = (abs(besselj(ord+1,c(1:N)))/kMax)';    %% m1 prepares input vector for transformation
    m2 = m1*kMax/rMax;                          %% m2 prepares output vector for display
    clear Jn
    clear Jm


    % Hankel transform over k iterated over each time element,
    % this then means we will have the following response:
    % y(r,\omega), then all we have to do after the hankel transform is done is
    % to invert the \omega component to retrieve it all back in (r,t) space.

    %     FWHMx = 0.25;
    %     sigmax = 1e-3*FWHMx/sqrt(log(2));
    %     FWHMt = 0.25;
    %     sigmat = FWHMt/sqrt(log(2));
    %     zw = exp(-(sigmat^2/4)*w.^2);       % prop to Fourier trans of exp(-t.^2/sigmat^2)

    Y1 = zeros(Nk, Nw);                 % initialize BOLD response of pole
    Y2 = zeros(Nk, Nw);
    Y3 = zeros(Nk, Nw);
    Y4 = zeros(Nk, Nw);
    Y5 = zeros(Nk, Nw);

    for j=1:length(w)

        [~, ~, T, kz] = PoleDecomposition_num(p, v_b, Gamma, k, w(j));

        if isstruct(stimulus)
            F1 = T.T1.*stimulus.T1(:,j); 
            F2 = T.T2.*stimulus.T2(:,j);
            F3 = T.T3.*stimulus.T3(:,j);
            F4 = T.T4.*stimulus.T4(:,j);
            F5 = T.T5.*stimulus.T5(:,j);
        else
            F1 = T.T1.*stimulus(:,j); 
            F2 = T.T2.*stimulus(:,j);
            F3 = T.T3.*stimulus(:,j);
            F4 = T.T4.*stimulus(:,j);
            F5 = T.T5.*stimulus(:,j);
        end

        FT1 = F1./m1;                   %% Prepare vector for transformation
        transformT1 = C*FT1;            %% Obtain the Hankel transform    
        fT1(:,j) = transformT1.*m2;     %% Prepare vector for display

        FT2 = F2./m1;                   %% Prepare vector for transformation
        transformT2 = C*FT2;            %% Obtain the Hankel transform    
        fT2(:,j) = transformT2.*m2;     %% Prepare vector for display

        FT3 = F3./m1;                   %% Prepare vector for transformation
        transformT3 = C*FT3;            %% Obtain the Hankel transform    
        fT3(:,j) = transformT3.*m2;     %% Prepare vector for display

        FT4 = F4./m1;                   %% Prepare vector for transformation
        transformT4 = C*FT4;            %% Obtain the Hankel transform    
        fT4(:,j) = transformT4.*m2;    %% Prepare vector for display

        FT5 = F5./m1;                   %% Prepare vector for transformation
        transformT5 = C*FT5;            %% Obtain the Hankel transform    
        fT5(:,j) = transformT5.*m2;    %% Prepare vector for display
    end

    for j=1:length(k)
        [Y1(j,:), ~] = ctfft(squeeze(fT1(j,:)),w);  % response of single pole 
        [Y2(j,:), ~] = ctfft(squeeze(fT2(j,:)),w);
        [Y3(j,:), ~] = ctfft(squeeze(fT3(j,:)),w);
        [Y4(j,:), ~] = ctfft(squeeze(fT4(j,:)),w);
        [Y5(j,:), ~] = ctfft(squeeze(fT5(j,:)),w);
    end

    Ytotal = Y1 + Y2 + Y3 + Y4 + Y5;        % total BOLD response Y(r,t)

    clear fT1 fT2 fT3 fT4 fT5 fTtotal

    Y = struct('Y1',Y1', 'Y2',Y2', 'Y3',Y3', 'Y4',Y4', ...
        'Y5',Y5', 'Ytotal',Ytotal', 'time', t, 'position', r);

end
