clear; clc;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
rng(12345);

%% Script purpose:
% This script shows how we can use external instruments to identify the
% responses of macro-variables to a shock of interest. The idea is that if
% we can find a time series zt which is correlated with a structural shock
% et, but uncorrelated with all other shocks in our VAR system, then we can
% use this instrument to construct impulse responses etc.

% In particular, we look at the case of so-called high-frequency
% instruments, as proposed for example, in Gertler and Karadi (2015) or
% Miranda-Agrippino and Ricco (2021) for the case of monetary policy
% shocks. Recall that shocks as we understand them are unobserved. However,
% one could argue that the movements in prices of certain assets right
% around the time when a monetary policy announcement is made can be seen
% as purely driven by the shock component of the policy. All expected
% aspects of the policy would have been priced in already, so if prices
% move, it must be because something unexpected -- a shock -- is happening.

% To guarantee that, besides RELEVANCE of the instrument (i.e. that the
% movements measured in the asset around the announcement reflect monetary
% policy shocks), we also have EXOGENEITY with respect to the other shocks
% (economic and financial) of the system, one needs to buy the assumption
% that policy announcements on a given day do not factor in economic or
% financial news of that specific day, only of the days before.

% A common shortcoming of the Cholesky approach in a VAR lies in the fact
% that the ordering matters. We typically assume that real variables
% respond with a lag to policy actions but policy can take innovations in
% real variables into account. However, financial variables can react
% immediately to the policy shock, while the financial innovations are
% likely to enter the policy decision as well. Hence, there is no place for
% financial variables in such a Cholesky VAR.

%% Loading the data

% Excess Bond Premium: 
% https://www.atlantafed.org/research/publications/policy-hub/2021/09/24/12--term-structure-of-excess-bond-premium

% Industrial Production (INDPRO), Short Rate (GS1),  CPI (CPIAUCSL) from
% FRED of St. Louis FED (+ some other series for trials)

% Wu and Xia (2016) shadow rate
% https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate

% Miranda-Agrippino and Ricco (2021) Instrument
% https://www.openicpsr.org/openicpsr/project/116841/version/V1/view

% Gertler and Karadi (2015) Instrument
% https://www.openicpsr.org/openicpsr/project/114082/version/V1/view

% Load the data
[data, ~] = xlsread('Proxydata.xlsx','Monthly','B2:E445');

% data goes from January 1973 to December 2009

IP = log(data(:,1))*100;
CPI = log(data(:,2))*100;
EBP = data(:,3);
SR = data(:,4);

finaldata = [IP, CPI, SR, EBP];
[T, N] = size(finaldata);

%% First we compute a Cholesky VAR for reference
c = 1;
%[p1, p2] = aicbic(finaldata,  24, c);
p = 12;
[beta, residuals] = VAR(finaldata, p, c);

% Compute Wold IRFs
hor = 48;
wold = woldirf(beta, c, p, hor);

% Compute Cholesky IRFs
sigma = (residuals' * residuals) ./ (T - 1 - p - N*p);
S = chol(sigma, 'lower');
cholirf = choleskyIRF(wold, S);

% Runkle (1987) Bootstrap for the Cholesky IRF
nboot = 1000;
prc = 68;
[bootchol, ~, ~, boot_beta] = bootstrapChol(finaldata,p,c,beta,residuals,nboot,hor,prc,'residual');

% Compute the bias-adjusted bootstrap of Kilian (1998)
shock = 3;
varnames = {'IP', 'CPI', 'Shadow Rate', 'EBP'};
shockname = "Monetary Policy Shock";
[upper_corrected, lower_corrected, corrections] = bootstrapChol_corrected(beta,boot_beta,p,c,finaldata,residuals,nboot,hor,prc,'residual');
plotirfchol(cholirf,upper_corrected,lower_corrected,shock,varnames,shockname, prc);


%% Now implement the IV-procedure using first the instrument of GK (2015)
% The authors preffered instrument is the surprise in the three month ahead
% futures rate (FF4)
ivGK = csvread('instruments_GK2015.csv',1);
ff4 = ivGK(127:end,4); % January 1990 to June 2012

%% First, we obtain an estimate of the reduced form residuals from the VAR
% that are related to the policy indicator. In our case it is the shadow
% rate which is ordered third.

u_p = residuals(:,3);

% Recall that we estimated these for 1973 to 2009 data with 12 lags,
% so we have to adapt both the instrument and the residuals 

ff4_final = ff4(1:end-30);
u_p_final = u_p(193:end);
u_q_final = residuals(193:end,[1,2,4]);

%% Second, we regress the reduced form estimate u_p on our instrument of the
% monetary policy shock (and a constant)

[~, uhat, ~] = myols(u_p_final, ff4_final,0);

% Now uhat can be interpreted as the variation in the VAR residual that
% comes from the structural policy shock.

%% Third, we regress the other residuals from the VAR onto uhat to get sq/sp

sq_sp = (uhat'*uhat)\uhat'*u_q_final;

% The ratio sq_sp tells us how strongly the other variables in the system
% correspond to the monetary policy shock relative to the response of the
% policy variable. Since s measures the responses ON IMPACT, it is
% convenient to require s_p = 1. which means that s_q is simply equal to
% the parameter estimated as sq_sp.

%% We normalize the initial impulse of the policy variable to 1

s =[sq_sp(1) sq_sp(2) 1 sq_sp(3)]';

% Now we can compute the structural IRFs as

ivirf = zeros(N, hor+1);

for h=1:hor+1
    ivirf(:,h) = wold(:,:,h)*s;
end

%% We can bootstrap confidence bands for this procedure.
% For this we need to bootstrap the VAR (we will use a Wild bootstrap for
% this as recommended in Gertler and Karadi 2015) and we need to bootstrap
% the instrument as well

[upper, lower, meanirf, medianirf] = ...
    bootstrapIV_corrected(finaldata,p,c,beta,residuals,ff4, 1000, 2000, [1,240], [193,432], 3, hor, 68);

%% Plot the IRFs
plotirf_partial(ivirf,upper,lower,varnames,shockname, prc)

%% FEVD Gertler and Karadi (2015)
% We can scale the variance of the shock and compute the sequence of the
% shock as follows:

scaler = s./norm(S\s);

% Unit variance shock IRF
ivirf_scaled = zeros(N, hor+1);

for h=1:hor+1
    ivirf_scaled(:,h) = wold(:,:,h)*scaler;
end

% Unit variance shock 
mon_pol_unit_var = scaler'*inv(sigma)*residuals';

% Check if this series has unit variance
check_uv = mon_pol_unit_var*mon_pol_unit_var'/(T - 1 - p - N*p);

%% Now we can compute the FEVD using this shock

% Get the total variation of the forecast errors
temp = 0;
denom = zeros(N,hor+1);
for h = 1:hor+1
    temp =  temp + wold(:,:,h) * sigma * wold(:,:,h)';
    denom(:,h) = diag(temp); % we only want the variances, not the covariances
end
    
% Get the variation in forecast errors due to a specific shock
irf_sq = ivirf_scaled.^2;
num = cumsum(irf_sq,2);

% Compute the share of the total
fevd_iv = num./denom;

% Plot the FEVD
plot_vardec(fevd_iv, varnames, shockname)
%% Alternative way: using Stock and Watson's (2012) method 

% Project: zt = delta'residuals + vt
[delta, zhat, vt] = myols(ff4_final,residuals(193:end,:),0);

% Compute the IRF
T1 = size(residuals(193:end,:),1);
sigma = residuals(193:end,:)'*residuals(193:end,:)/(T1 - 1 - p - N*p);
ivirf_scaled_sw = zeros(N, hor+1);

% Compute the unit variance structural shock
mp_uv_sw = delta'./sqrt(delta'*sigma*delta)*residuals(193:end,:)';
check_uv_sw = (mp_uv_sw*mp_uv_sw')/(T1 - 1 - p - N*p);


for h=1:hor+1
    ivirf_scaled_sw(:,h) = wold(:,:,h)*sigma*(delta./sqrt(delta'*sigma*delta));
end

for i = 1:N
    subplot(3,2,i)
    plot(ivirf_scaled_sw(i,:))
end

%% Now we compare this to the instrument of Miranda-Agrippino and Ricco (2021)
% This instrument goes from January 1991 to December 2009 (Hence the
% constraint chosen in 'finaldata')

ivMAR = xlsread('Proxydata.xlsx','Monthly','F218:F445');

% Plot the two instruments against each other
figure;
t1=datetime(1991,1,1);
t2=datetime(2009,12,31);
dates=t1:calmonths(1):t2;
plot(dates,ivMAR,'-','LineWidth',1.5,'Color','b'); hold on;
plot(dates,ff4_final(13:end),'-','LineWidth',1.5,'Color','r'); hold on;
line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;
ylabel('Surprise', 'FontSize', 16);   
title('Monetary Policy Instruments', 'FontSize', 16);  
set(gca,'FontSize',16)
axis tight
legend({'MAR2021','GK2015'})

%% First, we obtain an estimate of the reduced form residuals from the VAR
% that are related to the policy indicator. In our case it is the shadow
% rate which is ordered third.

u_p = residuals(:,3);

% Recall that we estimated these for 1973 to 2009 data with 12 lags,
% so we have to adapt the residuals 

u_p_final = u_p(205:end);
u_q_final = residuals(205:end,[1,2,4]);

%% Second, we regress the reduced form estimate u_p pn our instrument of the
% monetary policy shock (and a constant)

[~, uhat, ~] = myols(u_p_final, ivMAR,0);

% Now uhat can be interpreted as the variation in the VAR residual that
% comes from the structural policy shock.

%% Third, we regress the other residuals from the VAR onto uhat to get sq/sp

sq_sp = (uhat'*uhat)\uhat'*u_q_final;

% The ratio sq_sp tells us how strongly the other variables in the system
% correspond to the monetary policy shock relative to the response of the
% policy variable. Since s measures the respones ON IMPACT, it is
% convenient to require s_p = 1. which means that s_q is simply equal to
% the parameter estimated as sq_sp.

%% We normalize the initial impulse of the policy variable to 1

s =[sq_sp(1) sq_sp(2) 1 sq_sp(3)]';

% Now we can compute the structural IRFs as

ivirf = zeros(N, hor+1);

for h=1:hor+1
    ivirf(:,h) = wold(:,:,h)*s;
end

%% We can bootstrap confidence bands for this procedure.
% For this we need to bootstrap the VAR (we will use a Wild bootstrap for
% this as recommended in Gertler and Karadi 2015) and we need to bootstrap
% the instrument as well

[upper, lower, meanirf, medianirf] = ...
    bootstrapIV_corrected(finaldata,p,c,beta,residuals,ivMAR, 1000, 2000, [1,228], [205,432], 3, hor, prc);

%% Plot the IRFs
plotirf_partial(ivirf,upper,lower,varnames,shockname, prc)

%% Compute the unit variance shock sequence and the corresponding IRF

% If we want to conduct historical decompositions or FEVD we need a unit
% variance shock. Given our normalization above, (1) the shock itself is
% unknown and (2) the variance is not determined and fixed at 1. FOr HD we
% need the time series of the shock, for FEVD we need the IRF to the unit
% variance shock.

% We can scale the variance of the shock and compute the sequence of the
% shock as follows:

scaler = s./norm(S\s);

% Unit variance shock IRF
ivirf_scaled = zeros(N, hor+1);

for h=1:hor+1
    ivirf_scaled(:,h) = wold(:,:,h)*scaler;
end

% Unit variance shock 
mon_pol_unit_var = scaler'*inv(sigma)*residuals';

% Check if this series has unit variance
check_uv = mon_pol_unit_var*mon_pol_unit_var'/(T - 1 - p - N*p);

%% Now we can compute the FEVD using this shock

% Get the total variation of the forecast errors
temp = 0;
denom = zeros(N,hor+1);
for h = 1:hor+1
    temp =  temp + wold(:,:,h) * sigma * wold(:,:,h)';
    denom(:,h) = diag(temp); % we only want the variances, not the covariances
end
    
% Get the variation in forecast errors due to a specific shock
irf_sq = ivirf_scaled.^2;
num = cumsum(irf_sq,2);

% Compute the share of the total
fevd_iv = num./denom;

plot_vardec(fevd_iv, varnames, shockname)

%% Alternative way: using Stock and Watson's (2012) method 

% Project: zt = delta'residuals + vt
[delta, zhat, vt] = myols(ivMAR,residuals(205:end,:),0);

% Compute the IRF
T1 = size(residuals(205:end,:),1);
sigma = residuals(205:end,:)'*residuals(205:end,:)/(T1 - 1 - p - N*p);
ivirf_scaled_sw = zeros(N, hor+1);

% Compute the unit variance structural shock
mp_uv_sw = delta'./sqrt(delta'*sigma*delta)*residuals(205:end,:)';
check_uv_sw = (mp_uv_sw*mp_uv_sw')/(T1 - 1 - p - N*p);

for h=1:hor+1
    ivirf_scaled_sw(:,h) = wold(:,:,h)*sigma*(delta./sqrt(delta'*sigma*delta));
end

for i = 1:N
    subplot(3,2,i)
    plot(ivirf_scaled_sw(i,:))
end


%% Replication of: Diego Känzig -- The macroeconomic effects of oil supply news: Evidence from OPEC announcements
% American Economic Review 2021
% Idea: much like the central bank controls interest rates, OPEC controls
% the oil price to a near discretionary extent. Känzig uses OPEC
% announcements and the variation in oil futures contracts in a short
% window around these announcements to construct a proxy series for the oil
% news shock. This is assumed to be orthogonal to other souces of
% fluctutations and can therefore be used in the Proxy-SVAR machinery we
% have seen.

% Replication files can be found here: https://github.com/dkaenzig/replicationOilSupplyNews
% Copyright is of Diego Känzig

clear; clc; close all;

% Load the data for the VAR, these are monthly data from 1960M1 to 2017M12
% these data are in levels
load OilDataM.mat
dates = datetime('1960-01-01'):calmonths(1):datetime('2017-12-01');

% Load the data for the instrument
load OilSurprisesMLog.mat
datesproxy = datetime('1974-01-01'):calmonths(1):datetime('2017-12-01');

prox = 15; % We use the first principal component of the 1-12 months futures surprises
iv_kanzig = oilProxiesWTIM(:,prox);

% Specify the variables in the VAR. This includes the deflated oil price,
% world oil production, world oil inventories, world industrial production,
% US industrial production, and the US CPI
y = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 ...
    log(WORLDIP)*100 log(IP)*100 log(CPI)*100];

% Cut the data to the samples corresponding to the paper
% VAR: M1:1974 to M12:2017
% Proxy: M4:1983 to M12:2017
startVAR = find(dates == datetime('1974-01-01'));
endVAR = find(dates == datetime('2017-12-01'));

startPROX = find(datesproxy == datetime('1983-04-01'));
endPROX = find(datesproxy == datetime('2017-12-01'));

% Adjust the instrument and VAR data
iv_kanzig_final = iv_kanzig(startPROX:endPROX);
finaldata = y(startVAR:endVAR,:);
[T, N] = size(finaldata);

% Estimate the VAR
c = 1;
p = 12;
[beta, residuals] = VAR(finaldata, p, c);

% Truncate the residuals to the sample of the proxy to compute sigma to get
% variance representative of that sub-period

% Compute Wold IRFs
hor = 48;
wold = woldirf(beta, c, p, hor);

% Stage 1: Adjust size of VAR residuals to size of instrument
u_p = residuals(:,1); % The instrument should be related to oil prices
u_p_final = u_p(end-length(iv_kanzig_final)+1:end); % Use this way since they end in the same period Dec. 2017
u_q_final = residuals(end-length(iv_kanzig_final)+1:end,2:end);

% Compute Cholesky factor
u = residuals(end-length(iv_kanzig_final)+1:end, :);
T1 = size(u,1);
sigma = (u' * u) ./ (T1 - 1 - p - N*p);
S = chol(sigma, 'lower'); 

% Stage 2: We regress the reduced form estimate u_p pn our instrument of the
% oil news shock 
[~, uhat, ~] = myols(u_p_final, iv_kanzig_final,0);

% Stage 3: We regress the other residuals from the VAR onto uhat to get sq/sp
sq_sp = (uhat'*uhat)\uhat'*u_q_final;

% Stage 4: We normalize the initial impulse of the policy variable to 1
s =[1 sq_sp(2) sq_sp(3) sq_sp(4) sq_sp(5) sq_sp(5)]';

% Stage 5: Now we can compute the structural IRFs as:
ivirf = zeros(N, hor+1);

for h=1:hor+1
    ivirf(:,h) = wold(:,:,h)*s;
end

%% We can bootstrap confidence bands for this procedure.
% For this we need to bootstrap the VAR (we will use a Wild bootstrap for
% this as recommended in Gertler and Karadi 2015) and we need to bootstrap
% the instrument as well
prc = 90;
[upper, lower, meanirf, medianirf] = ...
    bootstrapIV_corrected(finaldata,p,c,beta,residuals,iv_kanzig_final, 1000, 2000, [1,417], [100,516], 1, hor, prc);

%% Plot the IRFs
varnames = {'Oil Price', 'World Oil Prod.','World Oil Inven.',  'World IP', 'US IP', 'US CPI'};
shockname = "Oil News";
plotirf_partial(meanirf,upper,lower,varnames,shockname, prc)

%% Scale to unit variance to obtain FEVD
scaler = s./norm(S\s);

% Unit variance shock IRF
ivirf_scaled = zeros(N, hor+1);

for h=1:hor+1
    ivirf_scaled(:,h) = wold(:,:,h)*scaler;
end

% Unit variance shock 
oil_news_unit_var = scaler'*inv(sigma)*u';

% Check if this series has unit variance
check_uv = oil_news_unit_var*oil_news_unit_var'/(T1 - 1 - p - N*p);

%% Now we can compute the FEVD using this shock

% Get the total variation of the forecast errors
temp = 0;
denom = zeros(N,hor+1);
for h = 1:hor+1
    temp =  temp + wold(:,:,h) * sigma * wold(:,:,h)';
    denom(:,h) = diag(temp); % we only want the variances, not the covariances
end
    
% Get the variation in forecast errors due to a specific shock
irf_sq = ivirf_scaled.^2;
num = cumsum(irf_sq,2);

% Compute the share of the total
fevd_iv = num./denom;

% Plot it
plot_vardec(fevd_iv,varnames,shockname)

%% Alternative way: using Stock and Watson's (2012) method 

% Project: zt = delta'residuals + vt
[delta, zhat, vt] = myols(iv_kanzig_final,u,0);

% Compute the unit variance structural shock
sigma = (u'*u)/(T1 - 1 - p - N*p);
oil_news_uv_sw = delta'./sqrt(delta'*sigma*delta)*u';
check_uv_sw = (oil_news_uv_sw*oil_news_uv_sw')/(T1 - 1 - p - N*p);

% Compute the IRF
ivirf_scaled_sw = zeros(N, hor+1);

for h=1:hor+1
    ivirf_scaled_sw(:,h) = wold(:,:,h)*sigma*(delta./sqrt(delta'*sigma*delta));
end

figure;
for i = 1:N
    subplot(3,2,i)
    plot(ivirf_scaled_sw(i,:))
end

%% Replication: Plagborg-Moller and Wolf -- Local Projections and VARs Estimate the Same Impulse Responses
% Econometrica 2021
% Idea: The main contribution of the paper is a demonstration that local
% projections and SVARs can be used interchangibly as they measure the same
% population moments, the IRF being one of them. There is another result
% which is an illustration how the Proxy-SVAR methodology can be used
% similarly to the narrative method as something the authors refer to as
% "internal instrument" SVAR. This is simply a Cholesky SVAR where the
% proxy is ordered first. Of course this means that we need as many
% observations of the proxy as of the other variables in the VAR. When we
% looked at MAR and GK, we used a larger pre-sample to obtain the residuals
% and were not too bothered with the instrument being a relatively short
% series. The following is merely meant to illustrate the internal
% instrument SVAR procedure, which is just a Cholesky SVAR.

% Order the instrument first in a Cholesky VAR
startVAR = find(dates == datetime('1983-04-01'));
endVAR = find(dates == datetime('2017-12-01'));

startPROX = find(datesproxy == datetime('1983-04-01'));
endPROX = find(datesproxy == datetime('2017-12-01'));

% Adjust the instrument and VAR data
iv_kanzig_final = iv_kanzig(startPROX:endPROX);
finaldata = [iv_kanzig_final, y(startVAR:endVAR,:)];

p = 12;
[beta, residuals] = VAR(finaldata, p, c);

% Compute Wold IRFs
hor = 48;
wold = woldirf(beta, c, p, hor);

% Compute Cholesky IRFs
sigma = (residuals' * residuals) ./ (T - 1 - p - N*p);
S = chol(sigma, 'lower');
cholirf = choleskyIRF(wold, S);

% Runkle (1987) Bootstrap for the Cholesky IRF
nboot = 1000;
prc = 68;
[bootchol, ~, ~, boot_beta] = bootstrapChol(finaldata,p,c,beta,residuals,nboot,hor,prc,'residual');

% Compute the bias-adjusted bootstrap of Kilian (1998)
[upper_corrected, lower_corrected, corrections] = bootstrapChol_corrected(beta,boot_beta,p,c,finaldata,residuals,nboot,hor,prc,'residual');

varnames = {'IV','Oil Price', 'World Oil Prod.','World Oil Inven.',  'World IP', 'US IP', 'US CPI'};
shockname = "Oil News";
shock = 1;
plotirfchol(cholirf,upper_corrected,lower_corrected,shock,varnames,shockname, prc);


% The responses are broadly the same as in the proxy approach, however, IP
% World and US IP are positively affected by the increase oil news shock.
% Oil prices increase so countries start depleting reserves. Once the
% acceleration in prices becomes smaller production gets ramped up again.
% The response in industrial production would merit some more
% disentangling: it could be that higher oil prices require subsitution to
% other forms of energy which may be beneficial for some industrialized
% nations (e.g. coal in India). The US economy itself could benefit from
% the higher prices as an oil producing country making oil production more
% profitable. In another paper Li, Plagborg-Moller and Wolf (2021)
% (https://arxiv.org/pdf/2104.00655.pdf) investigate the bias variance
% tradeoff form different estimation and identification methods and find
% that between the Proxy-SVAR and the internal-instrument SVAR approach the
% former may suffer from invertibility issues which can bias the estimates
% of the IRF at shorter horizons, but Proxy-SVARs have lower variance
% (they do not talk directly about inference, but this is likely to lead to
% tighter confidence bands, as we see in the IV applications).

%% Given that the Cholesky way does not give the same results as the IV
% procedure, two explanations come to mind: either, the sample for the
% Cholesky way is significantly different from the IV approach or the shock
% is non-fundamental. A test for this is described in Forni, Gambetti,
% Ricco (2022): External Instrument SVAR Analysis for Noninvertible Shocks:

clear; clc; close all;

% Load the data for the VAR, these are monthly data from 1960M1 to 2017M12
% these data are in levels
load OilDataM.mat
dates = datetime('1960-01-01'):calmonths(1):datetime('2017-12-01');

% Load the data for the instrument
load OilSurprisesMLog.mat
datesproxy = datetime('1974-01-01'):calmonths(1):datetime('2017-12-01');

prox = 15; % We use the first principal component of the 1-12 months futures surprises
iv_kanzig = oilProxiesWTIM(:,prox);

% Specify the variables in the VAR. This includes the deflated oil price,
% world oil production, world oil inventories, world industrial production,
% US industrial production, and the US CPI
y = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 ...
    log(WORLDIP)*100 log(IP)*100 log(CPI)*100];

% Cut the data to the samples corresponding to the paper
% VAR: M1:1974 to M12:2017
% Proxy: M4:1983 to M12:2017
startVAR = find(dates == datetime('1974-01-01'));
endVAR = find(dates == datetime('2017-12-01'));

startPROX = find(datesproxy == datetime('1983-04-01'));
endPROX = find(datesproxy == datetime('2017-12-01'));

% Adjust the instrument and VAR data
iv_kanzig_final = iv_kanzig(startPROX:endPROX);
finaldata = y(startVAR:endVAR,:);
[T, N] = size(finaldata);

% Estimate the VAR
c = 1;
p = 12;
[beta, residuals] = VAR(finaldata, p, c);


% Recall the sample adjusted residuals and the instrument
zt = iv_kanzig_final;
ut = residuals(end-length(iv_kanzig_final)+1:end,:);

% Regress zt onto its own past and obtain the residual
p = 12;
zy = zt(p+1:end);
zx = lagmakerMatrix(zt,p);
[~, ~, zhat] = myols(zy, zx, 0);

% Adjust ut accordingly
ut2 = ut(p+1:end,:);

% Create r leads of the residuals
r = 12;
treg = length(zhat);
zt_reg = zhat(1:end-r);
ut_reg = [];
ut_curr = ut2(1:end-r,:);
for rr = 1:r
    ut_reg = [ut_reg, ut2(1+rr:treg-r+rr,:)];
end
rhs = [ut_curr, ut_reg];

% Regress the proxy onto the current and future values of the VAR residuals
c = 0;
[delta, ~, ~, r2_u, fitted_partial] = myols(zt_reg, rhs, c);


% 1. Fundamentalness test: OLS regress zt onto current and future values of
% the VAR residuals. Apply F-test for joint significance. If we reject the
% Null that all coefficients are zero, the shock is non-fundamental.

% Regress zt on the restricted model where the RHS only has the
% contemporaneous values. For fundamentalness, we want the leads to be
% insignificant, ie. we do not want to reject the Null
[delta_aux, ~, ~, r2_r] = myols(zt_reg, ut_curr, c);
[T,m] = size(ut_reg);
df = T - size(rhs,2) - c;
[Fstat, pval] = do_ftest(r2_u, r2_r, m, df); % We cannot reject the Null of fundamentalness

% 2. Recoverability test: OLS regress zt onto current and future values of
% the VAR residuals. Apply a Ljung-Box test to the fitted value from this
% regression. If we reject the Null of no autocorrelation, the shock is
% not-recoverable. Then we can either add information to the VAR or be
% happy with relative IRFs. Otherwise, absolute measures of the unit
% variance shock, IRFs, FEVD, HD can be computed.
% Perform Ljung-Box test of autocorrelation

stdRes = normalize(fitted_partial);
autocorr(stdRes)
[h,pValue] = lbqtest(stdRes); % We reject the Null of no autocorrelation, and conclude that the absolute shock is not recoverable

% We conclude that the issue is likely just the difference in samples
