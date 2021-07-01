function [hat_x,Mtilde,info] = alg_BiasCorr(y,M,Sh,Sv,St,iterations,indLx)
% Plume bias correction method

% here, coupled here with the LSAPC algorithm: 
% Tichy, O., Smidl, V., Hofman, R., and Stohl, A.: LS-APC v1.0: a tuning-free method 
% for the linear inverse problem and its application to source-term determination, 
% Geosci. Model Dev., 9, 4297–4311, https://doi.org/10.5194/gmd-9-4297-2016, 2016.


[p n] = size(M);


% check the availability of gradients
if sum(sum(Sh)) == 0
    comp_h = 0;
else
    comp_h = 1;
end
if sum(sum(Sv)) == 0
    comp_v = 0;
else
    comp_v = 1;
end
if sum(sum(St)) == 0
    comp_t = 0;
else
    comp_t = 1;
end


%% settings of LSAPC algorithm
indL = zeros(n,n);
s = 2; 
r = 1; 

if 1
    blok = [ones(1,s) zeros(1,n/r-s)];
    indDiag = kron(ones(1,r),blok);
else
    blok = [ones(1,s) zeros(1,n/r-s)];
    blok2 = [ones(1,s-1) zeros(1,n/r-s+1)];
    indDiag = kron(ones(1,r) - [0 ones(1,r-1)],blok) + kron(ones(1,r) - [1 zeros(1,r-1)],blok2);
end
for j = 1:(n)
    if indDiag(j)
        indL = indL + diag(ones(n-j+1,1),j-1);
    end
end
indL = indL' - eye(n);

% start and prior
tlacitko = 1e0;
hat_upsilon = tlacitko*ones(n,1);
hat_varsigma = 1*indL;
hat_omega = 1*max(max(M'*M))^(-1);
hat_x = ones(n,1);
hat_L = eye(n);

alpha0 = 1e-10;
beta0 = 1e-10;
zeta0 = 1e-2;
eta0 = 1e-2;


koef_el0 = -1;
eL0 = zeros(n,n) + koef_el0*diag(ones(n-1,1),-1);

Dsi_lj = 0*ones(n,n);
pomLtxxtL = 0*eye(n);
pomLULt = 0*eye(n);


%% settings of BiasCorr algorithm
% start and prior
% set L_h, L_v a L_t 
kappa0hvt = 1e-10;
nu0hvt = 1e-10;

koefLstart = 0.0; % start L coefficient
koefLhvt0 = 0.0; % prior L coefficient

% horizontal
indLh = indLx;
hat_Lh = eye(p) + koefLstart*indLx;
eLh0 = zeros(p,p) + koefLhvt0*indLx;
zetah0 = 1e-2;
etah0 = 1e-2;
Dsi_lhi = zeros(p,p);
pomLhthhtLh = zeros(p,p);
pomLhWLht = zeros(p,p);
hat_varsigma_h = 1*indLx;
kappa0h = kappa0hvt;
nu0h = nu0hvt;

% vertical
indLv = indLx;
hat_Lv = eye(p) + koefLstart*indLx;
eLv0 = zeros(p,p) + koefLhvt0*indLx;
zetav0 = 1e-2;
etav0 = 1e-2;
Dsi_lvi = zeros(p,p);
pomLvtvvtLv = zeros(p,p);
pomLvWLvt = zeros(p,p);
hat_varsigma_v = 1*indLx;
kappa0v = kappa0hvt;
nu0v = nu0hvt;

% time
indLt = indLx;
hat_Lt = eye(p) + koefLstart*indLx;
eLt0 = zeros(p,p) + koefLhvt0*indLx;
zetat0 = 1e-2;
etat0 = 1e-2;
Dsi_lti = zeros(p,p);
pomLtttttLt = zeros(p,p);
pomLtWLtt = zeros(p,p);
hat_varsigma_t = 1*indLx;
kappa0t = kappa0hvt;
nu0t = nu0hvt;

hat_hh = zeros(p,1);
hat_hv = zeros(p,1);
hat_ht = zeros(p,1);

hat_wh = ones(p,1);
hat_wv = ones(p,1);
hat_wt = ones(p,1);

% omega
vartheta0 = 1e-10;
rho0 = 1e-10;

rozsah_h = 1.1; % coefficient of \delta s_{h,v,t} : h_{h,v,t} \in [-coef*\delta s_{h,v,t} ; +coef*\delta s_{h,v,t}]



%% auxiliary variables
info.omega = zeros(iterations,1);

HhSh = kron(hat_hh,ones(1,n)).*Sh;
HvSv = kron(hat_hv,ones(1,n)).*Sv;
HtSt = kron(hat_ht,ones(1,n)).*St;

hat_hhh = zeros(p,1);
hat_hvv = zeros(p,1);
hat_htt = zeros(p,1);
var_hh = zeros(p,1);
var_hv = zeros(p,1);
var_ht = zeros(p,1);

%% main algorithm
% LSAPC --> BiasCorr --> LSAPC --> BiasCorr --> LSAPC --> BiasCorr --> ...

MtMtilde = M'*M;
Mtilde = M;

konec = 1;
iter = 0;
while konec
    iter = iter + 1;
    if mod(iter,1) == 0
        display(['BiasCorr-LSAPC: iteration ' num2str(iter) '/' num2str(iterations)])    
    end
    
    %% LSAPC algorithm
    for it_lsapc = 1:100
        % x 
        si_x = ( hat_omega*( MtMtilde ) + hat_L*diag(hat_upsilon)*hat_L' + pomLULt)^(-1);
        mu_x = si_x*hat_omega*( Mtilde )'*y;

        if 1 % pozitivity on x
            [hat_x hat_xx] = momtrun_low(mu_x,sqrt(diag( si_x )),0);
            var_x = hat_xx - hat_x.^2;
            hat_xxt = hat_x*hat_x' + diag(var_x);
            if 1
                old_dSi = diag(si_x);
                new_dSi = var_x;
                nasobic = diag(sqrt(new_dSi./old_dSi));
                new_si_x = nasobic*si_x*nasobic;
                hat_xxt = hat_x*hat_x' + new_si_x;
            end
        else
            hat_x = mu_x;
            hat_xxt = mu_x*mu_x' + si_x;
        end

        % upsilon
        alpha = alpha0 + (1/2)*ones(n,1);
        for j = 1:(n-1)
            ind = find(indL(:,j)>0);
            pomLtxxtL(j,j) = hat_L([j; ind],j)'*hat_xxt([j; ind],[j; ind])*hat_L([j; ind],j) + trace(diag(Dsi_lj(ind,j))*hat_xxt(ind,ind)); 
        end
        pomLtxxtL(n,n) = hat_xxt(n,n);
        beta = beta0 + (1/2)*diag(pomLtxxtL);
        hat_upsilon = alpha./beta;

        si_LLt = zeros(n);
        pomLULt = zeros(n);
        pomLtxxtL = zeros(n);
        % l + varsigma
        for j = 1:(n-1)
            ind = find(indL(:,j)>0); 
            si_lj = (hat_upsilon(j)*hat_xxt(ind,ind) + diag(hat_varsigma(ind,j)) )^(-1);
            mu_lj = si_lj*(-hat_upsilon(j)*(hat_xxt(ind,j) ) + diag(hat_varsigma(ind,j))*eL0(ind,j));

            hat_lj = mu_lj;
            hat_ljljt = mu_lj*mu_lj' + si_lj;
            hat_ljtlj = mu_lj'*mu_lj + sum(diag(si_lj));

            hat_L(ind,j) = hat_lj;
            Dsi_lj(ind,j) = diag(si_lj);

            % varsigma(j)
            zeta = zeta0 + (1/2);
            eta = eta0 + (1/2)*diag(hat_ljljt) - diag(eL0(ind,j)*hat_lj') + (1/2)*diag(eL0(ind,j)*eL0(ind,j)');
            hat_varsigma(ind,j) = zeta./eta; 

            % LLt
            si_LLt(ind,ind) = si_lj;
            pomLULt(ind,ind) = pomLULt(ind,ind) + si_lj*hat_upsilon(j); 

        end

        % omega - inner for LSAPC
        vartheta = vartheta0 + p/2;
        rho = rho0 + (1/2)*y'*(y - 2*Mtilde*hat_x) + (1/2)*trace(hat_xxt*( MtMtilde ));
        hat_omega = vartheta/rho;
    end

    
    %% Plume bias correction method: BiasCorr
    
    % omega
    vartheta = vartheta0 + p/2;
    rho = rho0 + (1/2)*y'*(y - 2*Mtilde*hat_x) + (1/2)*trace(hat_xxt*( MtMtilde ));
    hat_omega = vartheta/rho;
    info.omega(iter) = hat_omega;
    
    % Hh a wh
    if comp_h
        
        si_hh = (hat_omega*diag(diag(Sh*hat_xxt*Sh')) + hat_Lh*diag(hat_wh)*hat_Lh' + pomLhWLht)^(-1); 
        mu_hh = si_hh*hat_omega * (diag(Sh*hat_x*y') - diag(Sh*hat_xxt*M') - diag(HvSv*hat_xxt*Sh') - diag(HtSt*hat_xxt*Sh'));
        
        [hat_hh,hat_hhh] = momtrun_tN(mu_hh,sqrt(diag(si_hh)),-rozsah_h,+rozsah_h);
        hat_hh = hat_hh';
        hat_hhh = hat_hhh';
        var_hh = hat_hhh - hat_hh.^2;
        hat_hhhht = hat_hh*hat_hh' + diag(var_hh);

        % w
        kappa_h = kappa0h + (1/2)*ones(p,1);
        for i = 1:(p-1)
            ind = find(indLh(:,i)>0);
            pomLhthhtLh(i,i) = hat_Lh([i; ind],i)'*hat_hhhht([i; ind],[i; ind])*hat_Lh([i; ind],i) + trace(diag(Dsi_lhi(ind,i))*hat_hhhht(ind,ind)); 
        end
        pomLhthhtLh(p,p) = hat_hhhht(p,p);
        
        nu_h = nu0h + (1/2)*diag(pomLhthhtLh);
        hat_wh = kappa_h./nu_h;
        
        % L + varsigma
        pomLhWLht = zeros(p,p);
        for i = 1:(p-1)
            ind = find(indLh(:,i)>0); 
            
            si_lhi = (hat_wh(i)*hat_hhhht(ind,ind) + diag(hat_varsigma_h(ind,i)))^(-1); % hat_hhhht(ind,i)
            mu_lhi = si_lhi*( -hat_wh(i)*hat_hhhht(ind,i) + diag(hat_varsigma_h(ind,i))*eLh0(ind,i) ); % overit rozmery!!!
            
            hat_lhi = mu_lhi;
            hat_lhilhit = mu_lhi*mu_lhi' + si_lhi;
            
            hat_Lh(ind,i) = hat_lhi;
            Dsi_lhi(ind,i) = diag(si_lhi);
            
            % varsigma
            zeta_h = zetah0 + 1/2;
            eta_h = etah0 + (1/2)*diag(hat_lhilhit);
            hat_varsigma_h(ind,i) = zeta_h./eta_h;
            
            pomLhWLht(ind,ind) = pomLhWLht(ind,ind) + si_lhi*hat_wh(i); 
        end
    else
        hat_hh = zeros(p,1);
        hat_wh = zeros(p,1);
    end
    HhSh = kron(hat_hh,ones(1,n)).*Sh;
    
    
    % Hv a wv
    if comp_v

        si_hv = (hat_omega*diag(diag(Sv*hat_xxt*Sv')) + hat_Lv*diag(hat_wv)*hat_Lv' + pomLvWLvt)^(-1);
        mu_hv = si_hv*hat_omega * (diag(Sv*hat_x*y') - diag(Sv*hat_xxt*M') - diag(HhSh*hat_xxt*Sv') - diag(HtSt*hat_xxt*Sv'));
        
        % h
        [hat_hv,hat_hvv] = momtrun_tN(mu_hv,sqrt(diag(si_hv)),-rozsah_h,+rozsah_h);
        hat_hv = hat_hv';
        hat_hvv = hat_hvv';
        var_hv = hat_hvv - hat_hv.^2;
        hat_hvhvt = hat_hv*hat_hv' + diag(var_hv);

        % w
        kappa_v = kappa0v + (1/2)*ones(p,1);
        for i = 1:(p-1)
            ind = find(indLv(:,i)>0);
            pomLvtvvtLv(i,i) = hat_Lv([i; ind],i)'*hat_hvhvt([i; ind],[i; ind])*hat_Lv([i; ind],i) + trace(diag(Dsi_lvi(ind,i))*hat_hvhvt(ind,ind)); 
        end
        pomLvtvvtLv(p,p) = hat_hvhvt(p,p);
        
        nu_v = nu0v + (1/2)*diag(pomLvtvvtLv);
        hat_wv = kappa_v./nu_v;
        
        % L + varsigma
        pomLvWLvt = zeros(p,p);
        for i = 1:(p-1)
            ind = find(indLv(:,i)>0);
            
            si_lvi = (hat_wv(i)*hat_hvhvt(ind,ind) + diag(hat_varsigma_v(ind,i)))^(-1);
            mu_lvi = si_lvi*( -hat_wv(i)*hat_hvhvt(ind,i) + diag(hat_varsigma_v(ind,i))*eLv0(ind,i) );
            
            hat_lvi = mu_lvi;
            hat_lvilvit = mu_lvi*mu_lvi' + si_lvi;
            
            hat_Lv(ind,i) = hat_lvi;
            Dsi_lvi(ind,i) = diag(si_lvi);
            
            % varsigma
            zeta_v = zetav0 + 1/2;
            eta_v = etav0 + (1/2)*diag(hat_lvilvit);
            hat_varsigma_v(ind,i) = zeta_v./eta_v;
            
            pomLvWLvt(ind,ind) = pomLvWLvt(ind,ind) + si_lvi*hat_wv(i); 
            
        end
    else
        hat_hv = zeros(p,1);
        hat_wv = zeros(p,1);
    end
    HvSv = kron(hat_hv,ones(1,n)).*Sv;    
    
    
    % Ht a wt
    if comp_t

        si_ht = (hat_omega*diag(diag(St*hat_xxt*St')) + hat_Lt*diag(hat_wt)*hat_Lt' + pomLtWLtt)^(-1);
        mu_ht = si_ht*hat_omega * (diag(St*hat_x*y') - diag(St*hat_xxt*M') - diag(HhSh*hat_xxt*St') - diag(HvSv*hat_xxt*St'));
        
        % h
        [hat_ht,hat_htt] = momtrun_tN(mu_ht,sqrt(diag(si_ht)),-rozsah_h,+rozsah_h);
        hat_ht = hat_ht';
        hat_htt = hat_htt';
        var_ht = hat_htt - hat_ht.^2;
        hat_hthtt = hat_ht*hat_ht' + diag(var_ht);
        
        % w
        kappa_t = kappa0t + (1/2)*ones(p,1);
        for i = 1:(p-1)
            ind = find(indLt(:,i)>0);
            pomLtttttLt(i,i) = hat_Lt([i; ind],i)'*hat_hthtt([i; ind],[i; ind])*hat_Lt([i; ind],i) + trace(diag(Dsi_lti(ind,i))*hat_hthtt(ind,ind)); 
        end
        pomLtttttLt(p,p) = hat_hthtt(p,p);
        
        nu_t = nu0t + (1/2)*diag(pomLtttttLt);
        hat_wt = kappa_t./nu_t;
        
        % L + varsigma
        pomLtWLtt = zeros(p,p);
        for i = 1:(p-1)
            ind = find(indLt(:,i)>0);
            
            si_lti = (hat_wt(i)*hat_hthtt(ind,ind) + diag(hat_varsigma_t(ind,i)))^(-1);
            
            mu_lti = si_lti*( -hat_wt(i)*hat_hthtt(ind,i) + diag(hat_varsigma_t(ind,i))*eLt0(ind,i) );
            
            hat_lti = mu_lti;
            hat_ltiltit = mu_lti*mu_lti' + si_lti;
            
            hat_Lt(ind,i) = hat_lti;
            Dsi_lti(ind,i) = diag(si_lti);
            
            % varsigma
            zeta_t = zetat0 + 1/2;
            eta_t = etat0 + (1/2)*diag(hat_ltiltit);
            hat_varsigma_t(ind,i) = zeta_t./eta_t;
            
            pomLtWLtt(ind,ind) = pomLtWLtt(ind,ind) + si_lti*hat_wt(i);
            
        end
    else
        hat_ht = zeros(p,1);
        hat_wt = zeros(p,1);
    end
    HtSt = kron(hat_ht,ones(1,n)).*St;        
    
    % output
    Mtilde = M + HhSh + HvSv + HtSt;
    MtMtilde = M'*M + (HhSh'*HhSh + Sh'*diag(var_hh)*Sh) + (HvSv'*HvSv + Sv'*diag(var_hv)*Sv) + (HtSt'*HtSt + St'*diag(var_ht)*St) + 2*M'*(HhSh+HvSv+HtSt) + 2*HhSh'*HvSv + 2*HhSh'*HtSt + 2*HvSv'*HtSt;

    if iter == iterations
        konec = 0;
    end
end

% bias correction field
info.hh = hat_hh;
info.hv = hat_hv;
info.ht = hat_ht;


end