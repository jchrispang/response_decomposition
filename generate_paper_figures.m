%% FIGURES FOR MANUSCRIPT 1

%% |Tj|^2 in +-k and +- \omega space

FigureNum = 2;      % figure number for manuscript
p = params;         % load nominal model parameters;

v_b = 1e-3;         % wave speed
Gamma = 1;          % wave damping rate

freqMax = 1;                % maximum temporal frequency
spatialFreqMax = 1000;      % maximum spatial frequency
freqMin = -freqMax;         % minimum temporal frequency
spatialFreqMin = -spatialFreqMax;   % minimum spatial frequency
 
% forming matrices of frequencies
[f, fk] = meshgrid(linspace(freqMin, freqMax, p.Nw/4), ...
                   linspace(spatialFreqMin, spatialFreqMax, p.Nk/4));   

kval = 2*pi*fk;
wval = 2*pi*f;

% calculating the decomposed transfer functions
[~, ~, T, ~] = PoleDecomposition_num(p, v_b, Gamma, kval, wval);

order = [1,3,5,2,4,6];
positions = {[0.10 0.55 0.22 0.33], [0.40 0.55 0.22 0.33], [0.70 0.55 0.27 0.33], ...
             [0.10 0.10 0.22 0.33], [0.40 0.10 0.22 0.33], [0.70 0.10 0.27 0.33]};

cmap = colormap('bone');
cmap_new = cmap(size(cmap,1)/4+1:end,:);

fig =figure('Position', [200, 200, 600, 400], 'Visible', 'on');
for i=1:6
    if i==6
        data = eval('T.Ttotal');
        plot_title = '$|T_{Y\zeta}|^2$';
    else
        data = eval(['T.T',num2str(order(i))]);
        if i==1
            plot_title = '$|T_1|^2$';
        elseif i==2
            plot_title = '$|T_2|^2$';
        elseif i==3
            plot_title = '$|T_3|^2$';
        elseif i==4
            plot_title = '$|T_4|^2$';
        elseif i==5
            plot_title = '$|T_5|^2$';
        end
    end
    
    ylim_max = 6000;
    
    subplot(2,3,i, 'Parent', fig, 'Position', positions{i})
    contourf(f, kval, log10(abs(data).^2))
    caxis([-7 1])
    colormap(flipud(cmap_new))
    xlim([freqMin, freqMax])
    ylim([-6000, 6000])
    set(gca, 'FontSize', 11, 'XTick',-1.0:0.5:1.0, 'YTick',-ylim_max:ylim_max/2:ylim_max)
    title(plot_title,'fontsize',15,'interpreter', 'latex');
    
    if i==1 || i==4
        ylabel('$k$ (m$^{-1}$)','fontsize',15,'interpreter', 'latex')
    end
    if i==4 || i==5 || i==6
        xlabel('$f$ (Hz)','fontsize',15,'interpreter', 'latex')
    end
    if i==3 || i==6
        colorbar('Ytick', -7:1)
    end
    
    % resonance curves
    if i==1
        line([-1 0], [ylim_max, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
        line([-1 0], [-ylim_max, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    elseif i==2
        line([-p.w_f/(2*pi) -p.w_f/(2*pi)], [ylim_max, -ylim_max], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    elseif i==3
        line([0 0], [ylim_max, -ylim_max], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    elseif i==4
        line([1 0], [ylim_max, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
        line([1 0], [-ylim_max, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    elseif i==5
        line([p.w_f/(2*pi) p.w_f/(2*pi)], [ylim_max, -ylim_max], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    end
    
    % subplot letter
    annotation('textbox', [0.018,0.86,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.018+0.3,0.86,0.1,0.1], 'String', '(c)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.018+2*0.3,0.86,0.1,0.1], 'String', '(e)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.018,0.86-0.45,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.018+0.3,0.86-0.45,0.1,0.1], 'String', '(d)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.018+2*0.3,0.86-0.45,0.1,0.1], 'String', '(f)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])


%% 2D Fourier transform of |Tj|^2 in r-t space

FigureNum = 3;      % figure number for manuscript
p = params;         % load nominal model parameters;

v_b = 1e-3;         % wave speed
Gamma = 1;          % wave damping rate

ord = 0;            % Bessel function order

load c.mat;         % pre-calculated Bessel function zeros
c = c(ord+1, 1:p.Nk+1);
k = c(1:p.Nk)'*2*pi*p.spatialFreqMax/c(p.Nk+1);       % k vector
r = c(1:p.Nk)/(2*pi*p.spatialFreqMax);                % radius vector
clear c

jvec = 0:p.Nw-1;
f = (jvec - p.Nw/2)/p.Nw*p.freqMax*2;
w = 2*pi*f;                                     % \omega

t = 1/(f(2)-f(1))*1/p.Nw*(jvec - p.Nw/2);       % time vector

[wmat, kmat] = meshgrid(w, k);          % forming matrices of frequencies

% calculating the decomposed transfer functions
[~, ~, T, ~] = PoleDecomposition_num(p, v_b, Gamma, kmat, wmat);

% Fourier transform of |T|^2
FT1 = FourierTransform2D(p, abs(T.T1).^2);
FT2 = FourierTransform2D(p, abs(T.T2).^2);
FT3 = FourierTransform2D(p, abs(T.T3).^2);
FT4 = FourierTransform2D(p, abs(T.T4).^2);
FT5 = FourierTransform2D(p, abs(T.T5).^2);
FTtotal = FourierTransform2D(p, abs(T.Ttotal).^2);

minTime = 0;
maxTime = 8;
maxSpace = 8;
clim_min = -0.5;
clim_max = 1;
% cmap1 = colormap('bone');
% cmap2 = colormap('hot');
% cmap_new = [cmap1(3*size(cmap1,1)/4+1:end,:); cmap1(end,:); flipud(cmap2(size(cmap2,1)/4+1:end,:))];
% cmap = colormap_helper('temp*');

% normalize with respect to maximum of FTtotal
basis_norm = real(FTtotal);     
normalization = max(basis_norm(:));
contours = 10;
tol = 1e-15;

positions_top = {[0.05 0.55 0.12 0.33], [0.2 0.55 0.12 0.33], ...
                 [0.35 0.55 0.12 0.33], [0.5 0.55 0.12 0.33], ...
                 [0.65 0.55 0.12 0.33], [0.8 0.55 0.18 0.33]};
letter_top = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
letter_bot = {'(g)', '(h)', '(i)', '(j)', '(k)', '(l)'};
             
fig = figure('Position', [200, 200, 600, 400], 'Visible', 'on');
for i=1:6
    if i==6
        data = eval('FTtotal');
        plot_title_top = 'Re[$F_{Y\zeta}$ ]';
        plot_title_bot = 'Im[$F_{Y\zeta}$ ]';
    else
        data = eval(['FT',num2str(i)]);
        if i==1
            plot_title_top = 'Re[$F_1$]';
            plot_title_bot = 'Im[$F_1$]';
        elseif i==2
            plot_title_top = 'Re[$F_2$]';
            plot_title_bot = 'Im[$F_2$]';
        elseif i==3
            plot_title_top = 'Re[$F_3$]';
            plot_title_bot = 'Im[$F_3$]';
        elseif i==4
            plot_title_top = 'Re[$F_4$]';
            plot_title_bot = 'Im[$F_4$]';
        elseif i==5
            plot_title_top = 'Re[$F_5$]';
            plot_title_bot = 'Im[$F_5$]';
        end
    end
    
    subplot(2,6,i, 'Parent', fig, 'Position', positions_top{i})
    
    contourf(r/1e-3, t, ...
        (real(data)/normalization).*((real(data)/normalization>=tol)+(real(data)/normalization<=-tol)),...
        contours); 
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [r(1)/1e-3, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', [r(1)/1e-3, 2:2:maxSpace], 'ytick', minTime:2:maxTime, ...
        'xticklabel', [0:2:maxSpace])
    title(plot_title_top,'fontsize',15, 'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87,0.1,0.1], 'String', letter_top{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
%     colormap(colormap(hot)/2 + colormap(cool)/2) 
    colormap(ThomasColorbar(clim_min, clim_max))    % Thomas colorbar
%     colormap(cmap);
    if i==1
        ylabel('$t$ (s)','fontsize',15, 'interpreter', 'latex')
    elseif i==6
        colorbar
    end
    
    subplot(2,6,6+i, 'Parent', fig, 'Position', positions_top{i}-[0 0.46 0 0])
    contourf(r/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)),...
        contours);  
    if strcmp(lastwarn, 'Contour not rendered for constant ZData')
        imagesc(r/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)));
        set(gca, 'YDir', 'normal')
    end 
        
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [r(1)/1e-3, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', [r(1)/1e-3, 2:2:maxSpace], 'ytick', minTime:2:maxTime, ...
        'xticklabel', [0:2:maxSpace])
    title(plot_title_bot, 'fontsize',15, 'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87-0.46,0.1,0.1], 'String', letter_bot{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    xlabel('$r$ (mm)','fontsize',15, 'interpreter', 'latex')
%     colormap(colormap(hot)/2 + colormap(cool)/2) 
    colormap(ThomasColorbar(clim_min, clim_max))
%     colormap(cmap)
    
    if i==1
        ylabel('$t$ (s)','fontsize',15, 'interpreter', 'latex')
    elseif i==6
        colorbar
    end
end
% 
set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])


%% Temporal power P^t_j(\omega) as a function of \omega

FigureNum = 4;      % figure number for manuscript
p = params;         % load nominal model parameters;

p.Nw = 4500;        % number of omega points
p.Nk = 3500;        % number of k points
p.freqMax = 1;      % maximum temporal frequency
p.spatialFreqMax = 4000;    % maximum spatial frequency

jvec = 0:p.Nw-1;
f = (jvec - p.Nw/2)/p.Nw*p.freqMax*2;
w = 2*pi*f;                             % \omega

kvec = 0:p.Nk-1;
fk = kvec/p.Nk*p.spatialFreqMax;
k = 2*pi*fk;                            % k

[wmat, kmat] = meshgrid(w, k);          % forming matrices of frequencies

v_b = [0.001, 0.0025, 0.005, 0.001, 0.0025, 0.005];
Gamma = [1, 1, 1, 0.6, 0.6, 0.6];

fig = figure('Position', [200, 200, 650, 400], 'Visible','on');
for i=1:length(v_b)
    [~, ~, T, ~] = PoleDecomposition_num(p, v_b(i), Gamma(i), kmat, wmat);

%     Pw_func = @(F) trapz(k, kmat.*abs(F).^2/(2*pi))/...
%         trapz(w, trapz(k, kmat.*abs(F).^2/(2*pi))/(2*pi)); % normalized
    Pw_func = @(F) trapz(k, kmat.*abs(F).^2/(2*pi)); % not normalized
    
    Pw = struct('T1', Pw_func(T.T1), 'T2', Pw_func(T.T2), ...
                'T3', Pw_func(T.T3), 'T4', Pw_func(T.T4), ...
                'T5', Pw_func(T.T5), ...
                'Ttotal', Pw_func(T.Ttotal));
            
    H = subplot(2,3,i, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    HAx = axes('Position', Hpos+[-0.065, -0.01, 0.02, 0.01]);
    
    color = bone;   
    Htot = semilogy(HAx, f, Pw.Ttotal, 'Color', color(1,:), 'linewidth', 1.5);
    hold on;
    H1 = semilogy(HAx, f, Pw.T1, 'Color', color(58,:), 'linewidth', 1.5);
    H2 = semilogy(HAx, f, Pw.T2, 'Color', color(48,:), 'linewidth', 1.5);
    H3 = semilogy(HAx, f, Pw.T3, 'Color', color(38,:), 'linewidth', 1.5);
    H4 = semilogy(HAx, f, Pw.T4, 'Color', color(28,:), 'linewidth', 1.5);
    H5 = semilogy(HAx, f, Pw.T5, 'Color', color(18,:), 'linewidth', 1.5);
    hold off;
    
    set(HAx, 'FontSize', 11)
    title(HAx, ['\nu_\beta=', sprintf('%0.1f', v_b(i)*1e3),' mm s^{-1}, \Gamma=',...
                sprintf('%0.1f', Gamma(i)),' s^{-1}'], ...
                'fontsize', 12);
    xlim(HAx, [-p.freqMax p.freqMax])
    ylim(HAx, [3*10^0, 3*10^7])
    set(HAx, 'xtick', [-p.freqMax, -p.freqMax/2, 0, p.freqMax/2, p.freqMax])
    set(HAx, 'ytick', [10^0, 10.^(1:7)])
    if i==1
        ylabel(HAx, '$P_j^t(f)$ (arb. units)', 'fontsize', 15, 'interpreter', 'latex')
        leg = legend(HAx, [H1, H2, H3, H4, H5, Htot], '$P^t_{1}(f)$', '$P^t_{2}(f)$', ...
                '$P^t_{3}(f)$', '$P^t_{4}(f)$', '$P^t_{5}(f)$', '$P^t_{Y\zeta}(f)$');
        set(leg, 'FontSize', 10, ...
            'Orientation','Vertical','position',[0.928,0.474,0.01,0.10], ...
            'interpreter', 'latex')
    elseif i==4
        label = ylabel(HAx, '$P_j^t(f)$ (arb. units)', 'fontsize', 15, ...
            'interpreter', 'latex');
    end
    if i>3
        xlabel(HAx, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
    end
end

% subplot letter
annotation('textbox', [0.065,0.83,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.065+0.28,0.83,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.065+2*0.28,0.83,0.1,0.1], 'String', '(c)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.065,0.83-0.48,0.1,0.1], 'String', '(d)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.065+0.28,0.83-0.48,0.1,0.1], 'String', '(e)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.065+2*0.28,0.83-0.48,0.1,0.1], 'String', '(f)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])


%% Spatial power P^s_j(k/2pi) as a function of k/2pi

FigureNum = 5;      % figure number for manuscript
p = params;         % load nominal model parameters;

p.Nw = 4500;        % number of omega points
p.Nk = 3500;        % number of k points
p.freqMax = 50;     % maximum temporal frequency
p.spatialFreqMax = 4000;    % maximum spatial frequency

jvec = 0:p.Nw-1;
f = (jvec - p.Nw/2)/p.Nw*p.freqMax*2;
w = 2*pi*f;                             % \omega

kvec = 0:p.Nk-1;
fk = kvec/p.Nk*p.spatialFreqMax;
k = 2*pi*fk;                            % k

[wmat, kmat] = meshgrid(w, k);          % forming matrices of frequencies

v_b = [0.001, 0.0025, 0.005, 0.001, 0.0025, 0.005];
Gamma = [1, 1, 1, 0.6, 0.6, 0.6];

fig = figure('Position', [200, 200, 650, 400], 'Visible','on');
for i=1:length(v_b)
    [~, ~, T, ~] = PoleDecomposition_num(p, v_b(i), Gamma(i), kmat, wmat);

%     Pk_func = @(F) trapz(w, abs(F).^2/(2*pi), 2)/...
%                 trapz(w, trapz(k, kmat.*abs(F).^2/(2*pi))/(2*pi)); % normalized
    Pk_func = @(F) trapz(w, abs(F).^2/(2*pi), 2); % not normalized
            
    Pk = struct('T1', Pk_func(T.T1), 'T2', Pk_func(T.T2), ...
                'T3', Pk_func(T.T3), 'T4', Pk_func(T.T4), ...
                'T5', Pk_func(T.T5), ...
                'Ttotal', Pk_func(T.Ttotal));
    
    color = bone;
    
    H = subplot(2,3,i, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    HAx = axes('Position', Hpos+[-0.055, -0.01, 0, 0.01]);
    
    Htot = loglog(k, Pk.Ttotal, 'Color', color(1,:), 'linewidth', 1.5);
    hold on;
    H1 = loglog(k, Pk.T1, 'Color', color(58,:), 'linewidth', 1.5);
    H2 = loglog(k, Pk.T2, 'Color', color(48,:), 'linewidth', 1.5);
    H3 = loglog(k, Pk.T3, 'Color', color(38,:), 'linewidth', 1.5);
    H4 = loglog(k, Pk.T4, 'Color', color(28,:), 'linewidth', 1.5);
    H5 = loglog(k, Pk.T5, 'Color', color(18,:), 'linewidth', 1.5);
    hold off;
    
    title(HAx, ['\nu_\beta=', sprintf('%0.1f', v_b(i)*1e3),' mm s^{-1}, \Gamma=',...
                sprintf('%0.1f', Gamma(i)),' s^{-1}'], ...
                'fontsize', 12);
            
    set(HAx, 'FontSize', 11)
    xlim(HAx, [0 10^4])
    ylim(HAx, [10^-8, 10^4])
    set(HAx, 'ytick', [10^-8, 10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4])
    set(HAx, 'xtick', [0, 10^0, 10^1, 10^2, 10^3, 10^4])
    if i==1
        ylabel(HAx, '$P^s_j(k)$ (arb. units)','fontsize',15, 'interpreter', 'latex')
        leg = legend(HAx, [H1, H2, H3, H4, H5, Htot], '$P^s_1(k)$', ...
            '$P^s_2(k)$', '$P^s_3(k)$', '$P^s_4(k)$', ...
            '$P^s_5(k)$', '$P^s_{Y\zeta}(k)$');
        set(leg, 'FontSize', 10, ...
            'Orientation','Vertical','position',[0.925,0.474,0.01,0.10], 'interpreter', 'latex')
    elseif i==4
        ylabel(HAx, '$P^s_j(k)$ (arb. units)','fontsize',15,'interpreter', 'latex')
    end
    if i>3
        xlabel(HAx, '$k$ (m$^{-1}$)','fontsize',15, 'interpreter', 'latex')
    end
end

% subplot letter
annotation('textbox', [0.076,0.83,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.076+0.28,0.83,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.076+2*0.28,0.83,0.1,0.1], 'String', '(c)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.076,0.83-0.48,0.1,0.1], 'String', '(d)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.076+0.28,0.83-0.48,0.1,0.1], 'String', '(e)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.076+2*0.28,0.83-0.48,0.1,0.1], 'String', '(f)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])


%% Behavior of \omega_1(k), \omega_2(k) as a function of k

FigureNum = 6;
p = params;         % Load nominal model parameters;

spatialFreqMin = 0;     
spatialFreqMax = 200;
fk = linspace(spatialFreqMin, spatialFreqMax, p.Nk);    
kval = 2*pi*fk;

v_b_vector = [1e-3, 2.5e-3, 5e-3, 1e-3, 2.5e-3, 5e-3];
Gamma_vector = [1.0, 1.0, 1.0, 0.6, 0.6, 0.6] ;

fig = figure('Position', [200, 200, 600, 400], 'Visible', 'on');

for i=1:6
    [w_pole, ~, ~, ~] = PoleDecomposition_num(p, v_b_vector(i), Gamma_vector(i), kval, 10);

    ylim_max = 5;
    ylim_min = -ylim_max;
    
    subplot(2,3,i)
    plot(kval, real(w_pole.w1),'-', 'Color', [0,0,0], 'linewidth', 2)
    hold on;
    plot(kval, real(w_pole.w2), '-', 'Color', [0.7,0.7,0.7], 'linewidth', 2)
    plot(kval, imag(w_pole.w1), '--', 'Color', [0,0,0], 'linewidth', 2)
    plot(kval, imag(w_pole.w2), '--', 'Color', [0.7,0.7,0.7], 'linewidth', 2)
    hold off;
    
%     plot(kval, real(w_pole.w1),'ko-', 'MarkerSize', 3)
%     hold on;
%     plot(kval, real(w_pole.w2), 'ks-', 'MarkerSize', 3)
%     plot(kval, imag(w_pole.w1), 'ko--', 'MarkerSize', 3)
%     plot(kval, imag(w_pole.w2), 'ks--', 'MarkerSize', 3)
%     hold off;
    xlabel('$k$ (m$^{-1}$)','fontsize',15,'interpreter', 'latex')
    
    set(gca, 'FontSize', 11, ...
        'xlim', [spatialFreqMin, 1000], 'ylim', [ylim_min ylim_max], ...
        'xtick', 0:1000/5:1000, 'ytick', ylim_min:ylim_max)
    text(spatialFreqMax/12, ylim_max*0.75, ...
        ['\nu_\beta=',sprintf('%0.1f', v_b_vector(i)*1e3),' mm s^{-1}, \Gamma=',...
            sprintf('%0.1f', Gamma_vector(i)),' s^{-1}'],'fontsize',10);

    if i==1
        leg = legend('Re[\omega_1]','Re[\omega_2]',...
            'Im[\omega_1]','Im[\omega_2]');
        set(leg, 'FontSize', 10, ...
            'Orientation','Horizontal','position',[0.51,0.96,0.01,0.01])
        ylabel('Re[\omega_j] & Im[\omega_j]','fontsize',15)
    elseif i==4
        ylabel('Re[\omega_j] & Im[\omega_j]','fontsize',15)
    end
end
annotation('textbox', [0.055,0.84,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.055+0.28,0.84,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.055+2*0.28,0.84,0.1,0.1], 'String', '(c)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.055,0.84-0.47,0.1,0.1], 'String', '(d)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.055+0.28,0.84-0.47,0.1,0.1], 'String', '(e)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')
annotation('textbox', [0.055+2*0.28,0.84-0.47,0.1,0.1], 'String', '(f)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])

%% Contour plot of kc

FigureNum = 7;
p = params;         % Load nominal model parameters;

v_b = linspace(1e-3, 20e-3, 2000);
Gamma= linspace(0.1, 1, 2000);
[v_b_mat, Gamma_mat] = meshgrid(v_b, Gamma); 

D  = p.rho_f.*(2.*Gamma_mat - p.beta*p.Cz/p.tau);  
kz = sqrt((p.k_0)^2 + 1./(v_b_mat.^2).*p.Cz.*(p.beta/p.tau).*(D./p.rho_f));

kc = (1./v_b_mat).*sqrt(Gamma_mat.^2 - kz.^2.*v_b_mat.^2);

kc_real = real(kc);
kc(find(~kc_real)) = -1*imag(kc(find(~kc_real)));

kclow = -225;%min(kc(:));
kcmid = 0;
kchigh = 1080;%max(kc(:));

% lowColor = [0 255 0];
% lowColor = [255 140 100];
lowColor = [255 140 100];
midColor = [255 255 255];
highColor = [0 0 0];

map = interp1( [kclow kcmid kchigh], ...
  [lowColor; midColor; highColor]/255, ...
  linspace(kclow, kchigh, 256));

contour_spacing = [linspace(kclow,0,8),linspace(0,kchigh,12)];
v_b1 = v_b(1);
Gamma1 = -p.k_0*v_b1 + p.Cz*p.beta/p.tau;
Gamma2 = Gamma(1);
v_b2 = (Gamma2 - p.Cz*p.beta/p.tau)/(-p.k_0);

fig = figure('Visible', 'on');
[C,h] = contourf(v_b*1e3,Gamma,kc,contour_spacing,'linewidth',0.5,'LineStyle','--');
hold on;
line([v_b1*1e3 v_b2*1e3],[Gamma1 Gamma2],'LineStyle','-', 'LineWidth',2,'color','b')
hold off;
caxis([kclow kchigh])
colormap(map)
cbar = colorbar;
xlim([1, 20])
ylim([0.1, 1])

tick1 = linspace(kclow,kcmid,4);
tick2 = linspace(kcmid,kchigh,13);
label1 = num2str(tick1',3);
label2 = num2str(tick2(2:end)', 3);

set(cbar, 'YTick', cat(2, tick1, tick2(2:end)))
set(cbar, 'YTickLabel', {'225i','150i','75i','0+0i','90','180','270','360','450',...
                         '540','630','720','810','900','990','1080'})
ylabel(cbar, '$k_c$','fontsize',15,'interpreter', 'latex')
xlabel('$\nu_\beta$ (mm s$^{-1}$)','fontsize',15,'interpreter', 'latex')
ylabel('$\Gamma$ (s$^{-1}$)','fontsize',15,'interpreter', 'latex')
set(gca, 'FontSize', 11)

% set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
% print('-depsc', ['Figure',num2str(FigureNum),'.eps'])

%% Modulus and Argument of a_j(k) as a function of k

FigureNum = 8;
p = params;         % Load nominal model parameters;

spatialFreqMin = 0;     
spatialFreqMax = 400;
fk = linspace(spatialFreqMin, spatialFreqMax, p.Nk);    
kval = 2*pi*fk;

v_b_vector = [1e-3, 1e-3];
Gamma_vector = [1.0, 0.6] ;

color = bone;
    
fig = figure('Visible', 'on');
for i=1:2
    [~, a, ~, kz] = PoleDecomposition_num(p, v_b_vector(i), Gamma_vector(i), kval, 1);
    
    subplot(2,2,2*(i-1)+1)
    semilogy(kval, abs(a.a1), 'Color', color(1,:),'linewidth',2)
    hold on;
    semilogy(kval, abs(a.a2), 'Color', color(20,:),'linewidth',2)
    semilogy(kval, abs(a.a3), 'Color', color(30,:),'linewidth',2)
    semilogy(kval, abs(a.a4), 'Color', color(40,:),'linewidth',2)
    semilogy(kval, abs(a.a5), 'Color', color(50,:),'linewidth',2)
    hold off;
    xlabel('$k$ (m$^{-1}$)','fontsize',15,'interpreter', 'latex')
    ylabel('$|a_j|$','fontsize',15,'interpreter', 'latex')
    set(gca, 'FontSize', 11, ...
        'xlim', [spatialFreqMin, 2000], 'ylim', [10^(-2), 10^(2)])
    text(-720, 0.05, ['\Gamma = ',sprintf('%0.1f', Gamma_vector(i)),' s^{-1}'], ...
        'fontsize', 18, 'fontweight', 'b', 'rotation', 90)
    
    subplot(2,2,2*(i-1)+2)
    plot(kval, angle(a.a1), 'Color', color(1,:),'linewidth',2)
    hold on;
    plot(kval, angle(a.a2), 'Color', color(25,:),'linewidth',2)
    plot(kval, angle(a.a3), 'Color', color(35,:),'linewidth',2)
    plot(kval, angle(a.a4), 'Color', color(45,:),'linewidth',2)
    plot(kval, angle(a.a5), 'Color', color(55,:),'linewidth',2)
    hold off;
    xlabel('$k$ (m$^{-1}$)','fontsize',15,'interpreter', 'latex')
    ylabel('Arg[$a_j$]','fontsize',15,'interpreter', 'latex')
    set(gca, 'FontSize', 11, ...
        'xlim', [spatialFreqMin, 2000], 'ylim', [-pi, pi])
    set(gca, 'yticklabel', {'-p', '-p/2', '0', 'p/2', 'p'}, ...
        'ytick', [-pi(), -pi()/2, 0, pi()/2, pi()],'fontname','symbol')
    
    if i==1
        leg = legend('$a_1$','$a_2$','$a_3$','$a_4$','$a_5$');
        set(leg, 'FontSize', 13, ...
            'Orientation','Horizontal','position',[0.51,0.97,0.01,0.01], 'fontname', 'Helvetica','interpreter', 'latex')
        annotation('textbox', [0.04,0.83,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
        annotation('textbox', [0.48,0.83,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
        annotation('textbox', [0.04,0.36,0.1,0.1], 'String', '(c)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
        annotation('textbox', [0.48,0.36,0.1,0.1], 'String', '(d)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    end
end

k12 = (1./v_b_vector(i)).*sqrt(Gamma_vector(i).^2 - kz.^2.*v_b_vector(i).^2);
k34 = (1./v_b_vector(i)).*sqrt(p.w_f^2 - (0.5*p.kappa - Gamma_vector(i))^2 + Gamma_vector(i)^2 - kz^2*v_b_vector(i)^2);
set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])


%% Real and imaginary part of 2D Green function G(r,t)

FigureNum = 9;
p = params;         % Load nominal model parameters;

v_b = 1e-3;
Gamma = 1;

% p.Nk = 2^11;
% p.Nw = 2^11;
% p.freqMax = 20;
% p.spatialFreqMax = 5000;
[Y, ~] = DecomposedBOLD_solution2D(p, v_b, Gamma, 0, 'dirac_delta', {});
r = Y.position;
t = Y.time;

minTime = 0;
maxTime = 8;
minSpace = 0;
maxSpace = 6;
clim_min = -0.5;
clim_max = 1;

basis_norm = real(Y.Ytotal);
normalization = max(basis_norm(:));
contour_spacing = clim_min:0.0909:clim_max;
% tol = 1e-20;
% contours = 10;
tol = 1e-5;
        
positions_top = {[0.05 0.55 0.12 0.33], [0.2 0.55 0.12 0.33], ...
                 [0.35 0.55 0.12 0.33], [0.5 0.55 0.12 0.33], ...
                 [0.65 0.55 0.12 0.33], [0.8 0.55 0.18 0.33]};
letter_top = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
letter_bot = {'(g)', '(h)', '(i)', '(j)', '(k)', '(l)'};
             
fig = figure('Position', [200, 200, 600, 400], 'Visible', 'on');
for i=1:6
    if i==6
        data = eval('Y.Ytotal');
        plot_title_top = 'Re[$G_{\rm sum}$ ]';
        plot_title_bot = 'Im[$G_{\rm sum}$ ]';
    else
        data = eval(['Y.Y',num2str(i)]);
        if i==1
            plot_title_top = 'Re[$G_1$]';
            plot_title_bot = 'Im[$G_1$]';
        elseif i==2
            plot_title_top = 'Re[$G_2$]';
            plot_title_bot = 'Im[$G_2$]';
        elseif i==3
            plot_title_top = 'Re[$G_3$]';
            plot_title_bot = 'Im[$G_3$]';
        elseif i==4
            plot_title_top = 'Re[$G_4$]';
            plot_title_bot = 'Im[$G_4$]';
        elseif i==5
            plot_title_top = 'Re[$G_5$]';
            plot_title_bot = 'Im[$G_5$]';
        end     
    end
        
    subplot(2,6,i, 'Parent', fig, 'Position', positions_top{i}) 
    contourf(r/1e-3, t, ...
        (real(data)/normalization).*((real(data)/normalization>=tol)+(real(data)/normalization<=-tol)),...
        [min(min(real(data)/normalization)), contour_spacing]);  
    if strcmp(lastwarn, 'Contour not rendered for constant ZData')
        imagesc(r/1e-3, t, ...
        (real(data)/normalization).*((real(data)/normalization>=tol)+(real(data)/normalization<=-tol)));
        set(gca, 'YDir', 'normal')
    end 
    colormap(ThomasColorbar(clim_min, clim_max))    % Thomas colorbar
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [r(1)/1e-3, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', [r(1)/1e-3, 2:2:maxSpace], 'ytick', minTime:2:maxTime, ...
        'xticklabel', [0:2:maxSpace])
    title(plot_title_top,'fontsize',15,'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87,0.1,0.1], 'String', letter_top{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    if i==1 || i==2
        A = real(data)/normalization;
        [~, max_ind] = max(A, [], 2);
        rmax = r(max_ind)/1e-3;
        [P, ~] = polyfit(rmax(513:633).', t(513:633), 1);
        hold on;
%         plot(rmax(513:633), t(513:633), 'k.')
%         plot(rmax(513:633).', polyval(P, rmax(513:633).'), 'k-', 'LineWidth', 1)
        hold off;
    end
        
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==6
        colorbar
    end
        
    subplot(2,6,6+i, 'Parent', fig, 'Position', positions_top{i}-[0 0.46 0 0])
    contourf(r/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)),...
        [min(min(imag(data)/normalization)), contour_spacing]);  
    if strcmp(lastwarn, 'Contour not rendered for constant ZData')
        imagesc(r/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)));
        set(gca, 'YDir', 'normal')
    end 
    colormap(ThomasColorbar(clim_min, clim_max))    % Thomas colorbar    
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [r(1)/1e-3, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', [r(1)/1e-3, 2:2:maxSpace], 'ytick', minTime:2:maxTime, ...
        'xticklabel', [0:2:maxSpace])
    title(plot_title_bot,'fontsize',15,'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87-0.46,0.1,0.1], 'String', letter_bot{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    xlabel('$r$ (mm)','fontsize',15,'interpreter', 'latex')
    
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==6
        colorbar
    end
end

set(fig, 'PaperPositionMode','auto')     % WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])

%% Real and imaginary part of 1D Green function G(x,t)

FigureNum = 10;
p = params;         % Load nominal model parameters;

v_b = 1e-3;
Gamma = 1;

[Y, kz] = DecomposedBOLD_solution1D(p, v_b, Gamma, 0, 'dirac_delta', {});
x = Y.position;
t = Y.time;

minTime = 0;
maxTime = 8;
minSpace = -6;
maxSpace = 6;
clim_min = -0.5;
clim_max = 1;

basis_norm = real(Y.Ytotal);
normalization = max(basis_norm(:));
contour_spacing = clim_min:0.0909:clim_max;
tol = 1e-10;
% contours = 10;
% tol = 1e-10;
        
positions_top = {[0.05 0.55 0.12 0.33], [0.2 0.55 0.12 0.33], ...
                 [0.35 0.55 0.12 0.33], [0.5 0.55 0.12 0.33], ...
                 [0.65 0.55 0.12 0.33], [0.8 0.55 0.18 0.33]};
letter_top = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
letter_bot = {'(g)', '(h)', '(i)', '(j)', '(k)', '(l)'};
             
fig = figure('Position', [200, 200, 600, 400], 'Visible', 'on');
for i=1:6
    if i==6
        data = eval('Y.Ytotal');
        plot_title_top = 'Re[$G_{\rm sum}$]';
        plot_title_bot = 'Im[$G_{\rm sum}$]';
    else
        data = eval(['Y.Y',num2str(i)]);
        if i==1
            plot_title_top = 'Re[$G_1$]';
            plot_title_bot = 'Im[$G_1$]';
        elseif i==2
            plot_title_top = 'Re[$G_2$]';
            plot_title_bot = 'Im[$G_2$]';
        elseif i==3
            plot_title_top = 'Re[$G_3$]';
            plot_title_bot = 'Im[$G_3$]';
        elseif i==4
            plot_title_top = 'Re[$G_4$]';
            plot_title_bot = 'Im[$G_4$]';
        elseif i==5
            plot_title_top = 'Re[$G_5$]';
            plot_title_bot = 'Im[$G_5$]';
        end
    end
        
    subplot(2,6,i, 'Parent', fig, 'Position', positions_top{i})        
    contourf(x/1e-3, t, ...
        (real(data)/normalization).*((real(data)/normalization>=tol)+(real(data)/normalization<=-tol)),...
        [min(min(real(data)/normalization)), contour_spacing]);  
    if strcmp(lastwarn, 'Contour not rendered for constant ZData')
        imagesc(x/1e-3, t, ...
        (real(data)/normalization).*((real(data)/normalization>=tol)+(real(data)/normalization<=-tol)));
        set(gca, 'YDir', 'normal')
    end 
    colormap(ThomasColorbar(clim_min, clim_max))    % Thomas colorbar  
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [minSpace, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', minSpace:2:maxSpace, 'ytick', minTime:2:maxTime)
    title(plot_title_top,'fontsize',15,'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87,0.1,0.1], 'String', letter_top{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    
    if i==1 || i==2
        A = real(data)/normalization;
        [~, max_ind_posx] = max(A(:,513:end), [], 2);
        xmax_pos = x(max_ind_posx+512)/1e-3;
        [~, max_ind_negx] = max(A(:,1:513), [], 2);
        xmax_neg = x(max_ind_negx)/1e-3;
        [Ppos, ~] = polyfit(xmax_pos(513:643), t(513:643), 1);
        [Pneg, ~] = polyfit(xmax_neg(513:643), t(513:643), 1);
        hold on;
%         plot(xmax_pos(513:643), t(513:643), 'k.')
%         plot(xmax_neg(513:643), t(513:643), 'k.')
        plot(x(513:end)/1e-3, x(513:end)/1e-3, 'k-', 'LineWidth', 1)
        plot(x(1:513)/1e-3, abs(x(1:513)/1e-3), 'k-', 'LineWidth', 1)
        
%         plot(xmax_pos(513:643), polyval(Ppos, xmax_pos(513:643)), 'k-', 'LineWidth', 1)
%         plot(xmax_neg(513:643), polyval(Pneg, xmax_neg(513:643)), 'k-', 'LineWidth', 1)
        hold off;
    end
    
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==6
        colorbar
    end
        
    subplot(2,6,6+i, 'Parent', fig, 'Position', positions_top{i}-[0 0.46 0 0])
    contourf(x/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)),...
        [min(min(imag(data)/normalization)), contour_spacing]);  
    if strcmp(lastwarn, 'Contour not rendered for constant ZData')
        imagesc(x/1e-3, t, ...
        (imag(data)/normalization).*((imag(data)/normalization>=tol)+(imag(data)/normalization<=-tol)));
        set(gca, 'YDir', 'normal')
    end 
    colormap(ThomasColorbar(clim_min, clim_max))    % Thomas colorbar  
    caxis([clim_min clim_max])
    set(gca, 'FontSize', 11, ...
        'xlim', [minSpace, maxSpace], 'ylim', [minTime, maxTime], ...
        'xtick', minSpace:2:maxSpace, 'ytick', minTime:2:maxTime)
    title(plot_title_bot,'fontsize',15,'interpreter', 'latex');
    annotation('textbox', [0.015+(i-1)*0.15,0.87-0.46,0.1,0.1], 'String', letter_bot{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==6
        colorbar
    end
end

set(fig, 'PaperPositionMode','auto')     % WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])

%% Real and imaginary parts of 1D solution to Gaussian stimuli
%%% Compares with experiment too

FigureNum = 12;
p = params;
subject = 'aLeft';

load singleSourceData/editedFileParams.mat      %% Loading v_b and Gamma parameters

load(cat(2,'singleSourceData/rawTS', subject, '.mat'))  %% Loading experiment results
load(cat(2,'singleSourceData/', subject, '.mat'))

% avgTSER : raw data
% smoothTSER : smooth data
v_b = Vmat(find(strcmp(associatedParams, subject)))*1e-3;     
Gamma = -Gmat(find(strcmp(associatedParams, subject))); 

v_b = 2.0120e-3;
Gamma = 0.7220;

t_exp = 0.25:0.25:20;     % 20 seconds with 250ms interval
x_exp = distanceInterp;
[t_exp_mat, x_exp_mat] = meshgrid(t_exp, x_exp*1e-3);

jvec = 0:p.Nw-1;
f = (jvec - p.Nw/2)/p.Nw*p.freqMax*2;
w = 2*pi*f;                                 % w = \omega

jveck = 0:p.Nk-1;
fk = (jveck - p.Nk/2)/p.Nk*p.spatialFreqMax*2;
k = 2*pi*fk;                                 % kx

t = 1/(f(2)-f(1))*1/p.Nw*(jvec - p.Nw/2);           % time vector
x = 1/(fk(2)-fk(1))*1/p.Nk*(jveck - p.Nk/2);        % position vector

[tt, xx] = meshgrid(t, x);          % mesh matrix of space and time
[wmat, kmat] = meshgrid(w, k);      % mesh matrix of spatial and temporal freq

sigma_x = (1/sqrt(log(2)))*1e-3;
sigma_t = 2/sqrt(log(2));
t0 = 5;

zeta_xt = (1/(pi*sigma_x*sigma_t))*exp(-xx.^2/(2*sigma_x^2)).*exp(-(tt-t0).^2/(2*sigma_t^2));

stimulus = 1/p.Nk * ifftshift(fftshift(ifft(fft(zeta_xt, [], 1), [], 2), 1), 2);

[Y, ~] = DecomposedBOLD_with_wintegral1D(p, v_b, Gamma, stimulus);

Ytotal_new = real(fftshift(Y.Ytotal.'));

tstart_ind = find(t==0);
tend_ind = find(t==20);
tend_ind2 = find(t==20);
xstart_ind = find(x==-5*1e-3);
xend_ind = find(x==5*1e-3);
xstart_ind2 = find(x==-10*1e-3);
xend_ind2 = find(x==10*1e-3);

cmap1 = bone;
cmap2 = hot;
cmap_new = [cmap1(24+1:end,:); flipud(cmap2(24+1:end,:))];

positions_top = {[0.06 0.38 0.14 0.18], [0.24 0.38 0.14 0.18], ...
                 [0.42 0.38 0.14 0.18], [0.60 0.38 0.14 0.18], ...
                 [0.78 0.38 0.22 0.18]};
letter_top = {'(c)', '(d)', '(e)', '(f)', '(g)'};
letter_bot = {'(h)', '(i)', '(j)', '(k)', '(l)'};

fig = figure('Position', [200, 200, 600, 400], 'Visible', 'on');
subplot(3,5,1, 'Parent', fig, 'Position', [0.24 0.73 0.20 0.21])
data_exp = smoothTSER;
contourf(x_exp, t_exp, data_exp.'/max(data_exp(:)), 10);   % NORMALIZED
colormap(cmap_new)
caxis([-1 1])
set(gca, 'FontSize', 11, ...
    'xlim', [-5, 5], 'ylim', [0, 20], ...
    'xtick', -5:2.5:5, 'ytick', 0:4:20)
title('experiment','fontsize',15);
xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
annotation('textbox', [0.16,0.9,0.1,0.1], 'String', '(a)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')

subplot(3,5,3, 'Parent', fig, 'Position', [0.54 0.73 0.28 0.21])
data_pred = Ytotal_new(xstart_ind:xend_ind,tstart_ind:tend_ind);
data_pred = data_pred - mean(data_pred(:));     % Correcting fluctuations around mean
contourf(x(xstart_ind:xend_ind)/1e-3, t(tstart_ind:tend_ind), ...
        (data_pred.'/max(data_pred(:))), 10);    % NORMALIZED
colormap(cmap_new)
caxis([-1 1])
set(gca, 'FontSize', 11, ...
    'xlim', [-5, 5], 'ylim', [0, 20], ...
    'xtick', -5:2.5:5, 'ytick', 0:4:20)
title('simulated','fontsize',15);
xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
colorbar
annotation('textbox', [0.46,0.9,0.1,0.1], 'String', '(b)', 'fontsize', 18, ...
    'fontweight', 'b', 'linestyle', 'none')

for i=1:5
    data = eval(['Y.Y',num2str(i)]);
    data = fftshift(data.');
    data = data(xstart_ind2:xend_ind2,tstart_ind:tend_ind2) - mean(data_pred(:));
    if i==1
        plot_title_top = 'Re[$Y_1$]';
        plot_title_bot = 'Im[$Y_1$]';
    elseif i==2
        plot_title_top = 'Re[$Y_2$]';
        plot_title_bot = 'Im[$Y_2$]';
    elseif i==3
        plot_title_top = 'Re[$Y_3$]';
        plot_title_bot = 'Im[$Y_3$]';
    elseif i==4
        plot_title_top = 'Re[$Y_4$]';
        plot_title_bot = 'Im[$Y_4$]';
    elseif i==5
        plot_title_top = 'Re[$Y_5$]';
        plot_title_bot = 'Im[$Y_5$]';
    end
    
    subplot(3,5,5+i, 'Parent', fig, 'Position', positions_top{i})
    contourf(x(xstart_ind2:xend_ind2)/1e-3, ...
             t(tstart_ind:tend_ind2), ...
             real(data).'/max(data_pred(:)), 10);
    colormap(cmap_new)
    caxis([-1 1])
    set(gca, 'FontSize', 11, ...
        'xlim', [-10, 10], 'ylim', [0, 20], ...
        'xtick', -10:5:10, 'ytick', 0:4:20)
    title(plot_title_top,'fontsize',15,'interpreter', 'latex');
    annotation('textbox', [0.021+(i-1)*0.18,0.54,0.1,0.1], 'String', letter_top{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==5
        colorbar
    end
    
    subplot(3,5,10+i, 'Parent', fig, 'Position', positions_top{i}-[0 0.29 0 0])
    if i==5
        imagesc(x(xstart_ind2:xend_ind2)/1e-3, ...
             t(tstart_ind:tend_ind2), ...
             imag(data).'/max(data_pred(:)));
    else
        contourf(x(xstart_ind2:xend_ind2)/1e-3, ...
             t(tstart_ind:tend_ind2), ...
             imag(data).'/max(data_pred(:)), 10);
    end
    colormap(cmap_new)
    caxis([-1 1])
    set(gca, 'FontSize', 11, ...
        'xlim', [-10, 10], 'ylim', [0, 20], ...
        'xtick', -10:5:10, 'ytick', 0:4:20)
    title(plot_title_bot,'fontsize',15,'interpreter', 'latex');
    xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    annotation('textbox', [0.021+(i-1)*0.18,0.25,0.1,0.1], 'String', letter_bot{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    
    if i==1
        ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    elseif i==5
        colorbar
    end
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['Figure',num2str(FigureNum),'.eps'])

