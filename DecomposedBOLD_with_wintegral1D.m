%%
% This file calculates the 1D hemodynamic decomposedresponse function, 
% as a function of radius and time with the w integral left for ctfft
 
% Dependencies:
% 1. params.m from J.C. Pang
% 2. PoleDecomposition_num.m from J.C. Pang

%%
function [Y, kz] = DecomposedBOLD_with_wintegral1D(p, v_b, Gamma, stimulus)

    % INPUT
    % 1. parameters: p
    % 2. wave properties: speed v_b, damping rate Gamma
    % 3. stimulus: stimulus in k, w space (Note that w is in x and k is in y)
    %    --> sample form [f, fk] = meshgrid(linspace(-freqMax, freqMax, Nw), ...
    %                    linspace(-spatialFreqMax, spatialFreqMax, Nk)); 
    %                    stimulus = f^2 + fk^2 
    %
    % OUTPUT
    % 1. structure Y containing BOLD response in x vs t of each pole contribution
    %
    % EXAMPLE: Y = DecomposedBOLD_with_wintegral1D(params, 1e-3, 1, stimulus)

    
    Nw = p.Nw;
    freqMax = p.freqMax;
    Nk = p.Nk;
    spatialFreqMax = p.spatialFreqMax;
    
    % predefine vectors for temporal frequency
    jvec = 0:Nw-1;
    f = (jvec - Nw/2)/Nw*freqMax*2;
    w = 2*pi*f;                                 %% w = \omega
    
    jveck = 0:Nk-1;
    fk = (jveck - Nk/2)/Nk*spatialFreqMax*2;
    k = 2*pi*fk;                                 %% kx
    
    t = 1/(f(2)-f(1))*1/Nw*(jvec - Nw/2);           %% time vector
    x = 1/(fk(2)-fk(1))*1/Nk*(jveck - Nk/2);        %% position vector
    
    [wmat, kmat] = meshgrid(w, k);
    
    wvals = (-1).^(1:length(w));
    kvals = (-1).^(1:length(k));
    
    wvals_mat = repmat(wvals, length(k), 1);        %% prepare matrix for Fourier transform
    kvals_mat = repmat(kvals.', 1, length(w));
    
    [~, ~, T, kz] = PoleDecomposition_num(p, v_b, Gamma, kmat, wmat);

    if isstruct(stimulus)
        F1 = T.T1.*stimulus.T1; 
        F2 = T.T2.*stimulus.T2;
        F3 = T.T3.*stimulus.T3;
        F4 = T.T4.*stimulus.T4;
        F5 = T.T5.*stimulus.T5;
    else
        F1 = T.T1.*stimulus; 
        F2 = T.T2.*stimulus;
        F3 = T.T3.*stimulus;
        F4 = T.T4.*stimulus;
        F5 = T.T5.*stimulus;
    end
    Y1 = kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*F1,[],2),[],1);
    Y2 = kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*F2,[],2),[],1);
    Y3 = kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*F3,[],2),[],1);
    Y4 = kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*F4,[],2),[],1);
    Y5 = kvals_mat.*ifft(kvals_mat.*wvals_mat.*fft(wvals_mat.*F5,[],2),[],1);

    Ytotal = Y1 + Y2 + Y3 + Y4 + Y5;        % total BOLD response Y(r,t)

    Y = struct('Y1',Y1', 'Y2',Y2', 'Y3',Y3', 'Y4',Y4', 'Y5',Y5',...
        'Ytotal',Ytotal', 'time', t, 'position', x);

end
