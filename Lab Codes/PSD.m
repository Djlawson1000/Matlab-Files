function [Sxx_mean, LLL] = PSD(x, fsamp, pointsperblock, smoothened, excludeFirstPoint)        
    fontsize = 30; LW = 2;

    % Time vector
    n = (1:length(x))./fsamp; 
    difft = diff(n); dt = difft(1); % Time step

    % Define blocks
    Nblocks = floor(length(x) ./ pointsperblock);           % Number of total blocks
    Blocktime = n(1:pointsperblock);                         % Time per block
    TOT_BlockTime = dt .* pointsperblock;                    % Total time per block

    % WINDOWING AND BLOCK AVERAGING 
    x_mat = reshape(x(1:pointsperblock .* Nblocks)', pointsperblock, Nblocks); % Reshape the signal to matrix
    window = .5*(1 - cos(2*pi*Blocktime/TOT_BlockTime))';                  % Hanning window
    win2 = window .* ones(1, Nblocks);                                      % Make window into matrix
    FT = fft(x_mat .* win2) .* (1 ./ pointsperblock) .* sqrt(8/3);         % FFT .* window .* 1/N .* corrector window
    Sxx = abs(FT).^2 .* TOT_BlockTime ;                                    % Matrix of power
    Sxx_mean = mean(Sxx, 2);                                                % Mean PSD

    % FREQUENCY VECTOR
    LL = (0:pointsperblock-1);  LLL = LL ./ pointsperblock .* fsamp;   
    deltaF = diff(LLL); 

    % Energy outputs of signal
    % sigma_square_EXACT = sum(mean(x_mat.^2, 2)) * 1 / (N);   % Exact
    % sigma_square_APPROX = 2 * deltaF(1) .* sum(Sxx_mean(1:round(length(LLL) / 2)));  % Approximation spectral
    Freq_Resolution = deltaF(10) / 1000;

    % If excludeFirstPoint is true, remove the first point from both Sxx_mean and LLL
    if excludeFirstPoint
        Sxx_mean = Sxx_mean(2:end);  % Exclude first point of PSD
        LLL = LLL(2:end);            % Exclude first frequency point
    end

    % Plot the PSD
    if smoothened == 1
        semilogy(LLL(1:round(length(LLL) / 2)) ./ 1000, (smooth((abs(Sxx_mean(1:round(length(LLL) / 2)))),'moving')),'-','linewidth',LW); hold on; grid on; % PSD
    else
        semilogy(LLL(1:round(length(LLL) / 2)) ./ 1000, (abs(Sxx_mean(1:round(length(LLL) / 2)))),'-','linewidth',LW); hold on; grid on; % PSD
    end
    set(gca, 'Box','on','FontName', 'Times','linewidth',1.9,'fontsize',fontsize); xlabel('$f$ (kHz)','Interpreter','latex'); ylabel('$S_{xx}(f)$ ([Pa]$^2$/Hz)','Interpreter','latex');
    set(gcf,'color','w'); 
    hold on;
end
