%%
% This file shows or plots the real and imaginary part of the 
% 2D hemodynamic decomposed response function as contours in a 
% radius and time space to various stimuli
 
% Dependencies:
% 1. params.m from J.C. Pang
% 2. PoleDecomposition_num.m from J.C. Pang
% 3. Hankel transform algorithm based on work provided by
%    M. Guizar-Sicairos and J. C. Gutierrez-Vega, Computation of 
%    quasi-discrete Hankel transforms of integer order for propagating optical 
%    wave fields, J. Opt. Soc. Am. A 21, 53-58 (2004).
% 4. c.mat from mathworks upload of (3)
% 5. ctfft.m from K. Aquino
% 6. DecomposedBOLD_with_wintegral2D.m from J.C. Pang

%%
function [Y, kz] = DecomposedBOLD_solution2D(p, v_b, Gamma, plotting, neuralDrive, vargin)
    
    % plotting = 0 -> no plot
    % plotting = 1 -> show plot but do not save
    % plotting = 2 -> don't show plot but save
    
    % For gaussian: vargin{1} = sigma_x, vargin{2} = sigma_t
    % For gabor: vargin{1} = omega_x, vargin{2} = omega_t, 
    %            vargin{3} = sigma_x, vargin{4} = sigma_t
    % For transfer_x: vargin{1} = which transfer function
    
    Nw = p.Nw;
    freqMax = p.freqMax;
    Nk = p.Nk;
    spatialFreqMax = p.spatialFreqMax;

    ord = 0;

    load c.mat;
    c = c(ord+1, 1:Nk+1);
    k = c(1:Nk)'*2*pi*spatialFreqMax/c(Nk+1);        %% k vector
    clear c

    jvec = 0:Nw-1;
    f = (jvec - Nw/2)/Nw*freqMax*2;
    w = 2*pi*f;                                 %% w = \omega

    [wmat, kmat] = meshgrid(w, k);

    switch neuralDrive
        case 'dirac_delta'
            stimulus = ones(Nk, Nw);
        case 'gaussian'
            switch length(vargin)
                case 0 % default widths of gaussian stimulus
                    FWHMx = 0.25;
                    sigmax = 1e-3*FWHMx/(2*sqrt(log(2)));
                    FWHMt = 0.25;
                    sigmat = FWHMt/(2*sqrt(log(2)));
                case 1
                    sigmax = 1e-3*vargin{1};
                    FWHMt = 0.25;
                    sigmat = FWHMt/(2*sqrt(log(2)));
                case 2
                    sigmax = 1e-3*vargin{1};
                    sigmat = vargin{2};
            end

            stimulus = exp(-(sigmax^2/4)*kmat.^2).*exp(-(sigmat^2/4)*wmat.^2);
        case 'gaussian_block'
            switch length(vargin)
                case 1 
                    FWHMx = 0.25;
                    sigmax = 1e-3*FWHMx/(2*sqrt(log(2)));
                    t_0 = vargin{1};
                case 2
                    sigmax = 1e-3*vargin{1};
                    t_0 = vargin{2};
            end

            stimulus = exp(-(sigmax^2/4)*kmat.^2).*(exp(1i*wmat*t_0)./(1i*wmat) - 1./(1i*wmat));
        case 'gabor'
            switch length(vargin)
                case 2 % default widths of gaussian stimulus
                    omega_x = vargin{1};
                    omega_t = vargin{2};
                    FWHMx = 0.25;
                    sigmax = 1e-3*FWHMx/(2*sqrt(log(2)));
                    FWHMt = 0.25;
                    sigmat = FWHMt/(2*sqrt(log(2)));
                case 3
                    omega_x = vargin{1};
                    omega_t = vargin{2};
                    sigmax = 1e-3*vargin{3};
                    FWHMt = 0.25;
                    sigmat = FWHMt/(2*sqrt(log(2)));
                case 4
                    omega_x = vargin{1};
                    omega_t = vargin{2};
                    sigmax = 1e-3*vargin{3};
                    sigmat = vargin{4};
            end

            stimulus = exp(-(sigmax^2/4)*(-omega_x+kmat).^2).*exp(-(sigmat^2/4)*(-omega_t+wmat).^2);
        case 'gaussian_separated'
            [w_pole, ~, ~, ~] = PoleDecomposition_num(p, v_b, Gamma, kmat, wmat);

            FWHMx = 0.25;
            sigmax = 1e-3*FWHMx/(2*sqrt(log(2)));
            FWHMt = 0.25;
            sigmat = FWHMt/(2*sqrt(log(2)));

            stim = @(wpole) exp(-(sigmax^2/4)*kmat.^2).*exp(-(sigmat^2/4)*wpole.^2);

            stimulus = struct('T1',stim(w_pole.w1), 'T2',stim(w_pole.w2), ...
                              'T3',stim(w_pole.w3), 'T4',stim(w_pole.w4), ...
                              'T5',stim(w_pole.w5));
        case 'transfer_raw'
            [~, ~, transfer, ~] = PoleDecomposition_num(p, v_b, Gamma, kmat, wmat);
            
            switch vargin{1}
                case '1'
                    stimulus = (transfer.T1);
                case '2'
                    stimulus = (transfer.T2);
                case '3'
                    stimulus = (transfer.T3);
                case '4'
                    stimulus = (transfer.T4);
                case '5'
                    stimulus = (transfer.T5);
                case 'all'
                    stimulus = (transfer.Ttotal);
            end
        case 'transfer_abs'
            [~, ~, transfer, ~] = PoleDecomposition_num(p, v_b, Gamma, kmat, wmat);
            
            switch vargin{1}
                case '1'
                    stimulus = abs(transfer.T1).^2;
                case '2'
                    stimulus = abs(transfer.T2).^2;
                case '12'
                    stimulus = abs(transfer.T1).^2 + abs(transfer.T2).^2;
                case '3'
                    stimulus = abs(transfer.T3).^2;
                case '4'
                    stimulus = abs(transfer.T4).^2;
                case '34'
                    stimulus = abs(transfer.T3).^2 + abs(transfer.T4).^2;
                case '5'
                    stimulus = abs(transfer.T5).^2;
                case 'all'
                    stimulus = abs(transfer.Ttotal).^2;
            end
    end

    [Y, kz] = DecomposedBOLD_with_wintegral2D(p, v_b, Gamma, stimulus);
    r = Y.position;
    t = Y.time;
    
    %% PLOTTING
    
    if plotting>0
        maxTime = 15;
        minSpace = r(1)/1e-3;
        maxSpace = 3;
        clim_min = -0.5;
        clim_max = 2;

        basis_norm = real(Y.Ytotal);
        normalization = max(basis_norm(:));
        contour_spacing = clim_min:0.096:clim_max;
        tol = 1e-10;

        fig1 = figure('Position', [200, 200, 800, 500], 'Visible','off');

        subplot(2,3,1, 'Parent', fig1, 'Position',[0.1 0.55 0.22 0.33]);
        F = real(Y.Y1);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        ylabel('t (s)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.03,0.79,0.1,0.1], 'String', '(a)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{1}(r, t) ]','fontsize',15);

        subplot(2,3,4, 'Parent', fig1, 'Position',[0.1 0.1 0.22 0.33]);
        F = real(Y.Y2);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        ylabel('t (s)','fontsize',15)
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.03,0.34,0.1,0.1], 'String', '(b)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{2}(r, t) ]','fontsize',15);

        subplot(2,3,2, 'Parent', fig1, 'Position',[0.4 0.55 0.22 0.33]);
        F = real(Y.Y3);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.335,0.79,0.1,0.1], 'String', '(c)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{3}(r, t) ]','fontsize',15);

        text(-0.75, maxTime*1.25, ...
            ['\nu_\beta=',sprintf('%0.1f', v_b*1e3),'mm/s,  k_z\nu_\beta=',...
                sprintf('%0.2f', kz*v_b),'/s,  \Gamma=',...
                sprintf('%0.1f', Gamma),'/s'],...
                        'fontsize',18); % title for entire figure

        subplot(2,3,5, 'Parent', fig1, 'Position',[0.4 0.1 0.22 0.33]);
        F = real(Y.Y4);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.335,0.34,0.1,0.1], 'String', '(d)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{4}(r, t) ]','fontsize',15);

        subplot(2,3,3, 'Parent', fig1, 'Position',[0.7 0.55 0.27 0.33]);
        F = real(Y.Y5);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.635,0.79,0.1,0.1], 'String', '(e)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{5}(r, t) ]','fontsize',15);
        colorbar

        subplot(2,3,6, 'Parent', fig1, 'Position',[0.7 0.1 0.27 0.33]);
        F = real(Y.Ytotal);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.635,0.34,0.1,0.1], 'String', '(f)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Re[ Y_{total}(r, t) ]','fontsize',15);
        colorbar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fig2 = figure('Position', [200, 200, 800, 500], 'Visible','off');

        subplot(2,3,1, 'Parent', fig2, 'Position',[0.1 0.55 0.22 0.33]);
        F = imag(Y.Y1);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        ylabel('t (s)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.03,0.79,0.1,0.1], 'String', '(a)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{1}(r, t) ]','fontsize',15);

        subplot(2,3,4, 'Parent', fig2, 'Position',[0.1 0.1 0.22 0.33]);
        F = imag(Y.Y2);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        ylabel('t (s)','fontsize',15)
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.03,0.34,0.1,0.1], 'String', '(b)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{2}(r, t) ]','fontsize',15);

        subplot(2,3,2, 'Parent', fig2, 'Position',[0.4 0.55 0.22 0.33]);
        F = imag(Y.Y3);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end  
        caxis([clim_min clim_max])
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.335,0.79,0.1,0.1], 'String', '(c)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{3}(r, t) ]','fontsize',15);

        text(-0.75, maxTime*1.25, ...
            ['\nu_\beta=',sprintf('%0.1f', v_b*1e3),'mm/s,  k_z\nu_\beta=',...
                sprintf('%0.2f', kz*v_b),'/s,  \Gamma=',...
                sprintf('%0.1f', Gamma),'/s'],...
                        'fontsize',18); % title for entire figure

        subplot(2,3,5, 'Parent', fig2, 'Position',[0.4 0.1 0.22 0.33]);
        F = imag(Y.Y4);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end    
        caxis([clim_min clim_max])
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.335,0.34,0.1,0.1], 'String', '(d)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{4}(r, t) ]','fontsize',15);

        subplot(2,3,3, 'Parent', fig2, 'Position',[0.7 0.55 0.27 0.33]);
        F = imag(Y.Y5);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.635,0.79,0.1,0.1], 'String', '(e)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{5}(r, t) ]','fontsize',15);
        colorbar
        set(gca,'YDir','normal')

        subplot(2,3,6, 'Parent', fig2, 'Position',[0.7 0.1 0.27 0.33]);
        F = imag(Y.Ytotal);
        contourf(r/1e-3, t, ...
          ((F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol))), ...
          [min(min(F/normalization)), contour_spacing]); 
        if strcmp(lastwarn, 'Contour not rendered for constant ZData')
            imagesc(r/1e-3, t, ...
            (F/normalization).*((F/normalization>=tol)+(F/normalization<=-tol)));
            set(gca,'YDir','normal')
        end 
        caxis([clim_min clim_max])
        xlabel('r (mm)','fontsize',15)
        xlim([minSpace, maxSpace])
        ylim([0, maxTime])
        set(gca, 'FontSize', 11)
        annotation('textbox', [0.635,0.34,0.1,0.1], 'String', '(f)', 'fontsize', 20, ...
        'fontweight', 'b', 'linestyle', 'none')
        title('Im[ Y_{total}(r, t) ]','fontsize',15);
        colorbar
        set(gca,'YDir','normal')
        
        
        if plotting==1
            figure(fig1)
            figure(fig2)
        elseif plotting==2
        %%saving plot
            folder = ['Y_w=',num2str(freqMax),'_k=',num2str(spatialFreqMax),'/'];
            switch length(vargin)
                case 0
                    filename1 = [folder,'RealY_',neuralDrive,...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
                    filename2 = [folder,'ImagY_',neuralDrive,...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
                case 1
                    filename1 = [folder,'RealY_',neuralDrive,'_vargin1=',vargin{1},...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
                    filename2 = [folder,'ImagY_',neuralDrive,'_vargin1=',vargin{1}...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
                case 2
                    filename1 = [folder,'RealY_',neuralDrive,'_vargin1=',vargin{1},'_vargin2=',vargin{2},...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
                    filename2 = [folder,'ImagY_',neuralDrive,'_vargin1=',vargin{1},'_vargin2=',vargin{2},...
                         '_w=',num2str(freqMax),'_Nw=',num2str(Nw),...
                         '_k=',num2str(spatialFreqMax),'_Nk=',num2str(Nk),...
                         '_vb=',num2str(v_b),'_Gamma=',num2str(Gamma)];
            end
            set(fig1, 'PaperPositionMode','auto')     %# WYSIWYG
            print(fig1, '-depsc', [filename1,'.eps'])
            display('Figure saved')
            
            set(fig2, 'PaperPositionMode','auto')     %# WYSIWYG
            print(fig2, '-depsc', [filename2,'.eps'])
            display('Figure saved')
        end
    end
end