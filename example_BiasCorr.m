
% Example run of BiasCorr algorithm
% (here, coupled with the LSAPC linear inverse algorithm)
%
% paper: Source term determination with elastic plume bias correction
% authors: Ondrej Tichy, Vaclav Smidl, Nikolaos Evangeliou
% contact: otichy@utia.cas.cz

clear all;
close all;
warning off all;


%% test data load
    % load example data or create your own
    if 1
        load ./data_example.mat
    else
        p = 200;
        n = 10;
        x_true = [0 0 0 1 1 1 1 0 0 0]';
        hh = zeros(p,1);
        hv = zeros(p,1);
        ht = zeros(p,1);

        eM = rand(p,n+2); 
        eM(eM<0) = 0;

        % simulated gradient SRS matrices and elastic coefficients
        M_h = randn(p,n);
        hh([50:120]) = -1; 

        M_v = randn(p,n);
        hv([10:70]) = +1;

        M_t = eM(:,1:(end-2)) - eM(:,3:end);
        ht(140:180) = -1;

        % simulated SRS matrix
        M = eM(:,2:(end-1));

        % measurement and noise
        y = (M + diag(hh)*M_h + diag(hv)*M_v + diag(ht)*M_t)*x_true;
        e = randn(p,1)*0.1;
        y = y + e;
        
%         save ./data_example.mat M M_h M_v M_t y x_true
    end
    
    % end up with measurement vector y, SRS matrix M and gradient matrices M_h, M_v, and M_t

%% setting and run of BiasCorr algorithm coupled with selected linear inverse algorithm
    
    % create index set \mathcal{I} of neighborhood of each measurement coded in indLx matrix, here simple subdiagonal for simplicity
    p = length(y);
    indLx = zeros(p,p) + diag(ones(p-1,1),-1);
    

    iterations = 10; % select number of iterations

    % main algorithm: BiasCorr method, here coupled with the LSPAC algorithm 
    [hat_x,Mtilde,info] = alg_BiasCorr(y,M,M_h,M_v,M_t,iterations,indLx);
        
%% print estimates and scatter plot

    fig = figure(1);
    set(fig, 'Position', [0, 1000, 900, 250]);
    rows = 1;
    cols = 3;
    subplot(rows,cols,[1:2])
        stairs(hat_x,'blue','Linewidth',2)
        hold on
        stairs(x_true,'red--')
        hold off
        xlim([1 10])
        legend('BiasCorr-LSAPC method','ground truth','Location','southeast')
    subplot(rows,cols,3)
        plot(y,M*x_true,'x','MarkerEdgeColor','red')
        hold on
        plot(y,(Mtilde)*hat_x,'x','MarkerEdgeColor','blue')
        hold on
        plot([0 10],[0 10],'black--')
        hold off
        axis('square')
        lim1 = -0;
        lim2 = max([y;M*x_true;Mtilde*hat_x]);
        xlim([lim1 lim2])
        ylim([lim1 lim2])
        xlabel('measurements')
        ylabel('reconstructions')

