function n_xt = n_scott_gamma_kernel(x,t,nu,B,gam,kswitch)
% Scott's 1968 analytical solutions to 
% collection equation for Gamma (Gaussian-like)
% initial distribution and Golovin kernel.

xt = 10;



% Golovin's Kernel
if kswitch == 0
tau = 1 - exp(-0.00153.*t); % normalized time

%gam_series_old = zeros(length(x),length(t));
gam_series_old = zeros(1,length(t));
%gam_series = 100;
gam_tol = 1;
gam_tol_max = 1e-5;
k = 0;

%syms x

%digits(500)

%x_new = vpa(x);

x_new = x;

%if x_new < xt

while gam_tol > gam_tol_max
%while k <100 

gam_series  = gam_series_old +...
tau.^(k) .*nu.^(nu.*(k+1))./...
(factorial(k+1).*gammaz(nu.*(k+1))) .* x_new.^(k.*(nu+1));
    
%gam_series  = vpa(gam_series_old +...
%tau.^(k) .* x.^(k.*(nu+1)) .*nu.^(nu.*(k+1))./...
%(gammaz(k+2).*gammaz(nu.*(k+1))));

gam_tol = ...
max(abs(gam_series_old-gam_series)./gam_series_old,[],'omitnan');

gam_series_old = gam_series;

k = k + 1;
end


n_xt = (1-tau) .* x_new.^(nu-1) .* exp(-(nu+tau).*x_new) .*gam_series_old;

%else
    

if x_new>=xt
    
    n_xt= (1-tau).*exp(-(tau+nu).*...
                x_new+(nu+1).*x_new.*tau.^(1./(nu+1)))./(x_new.^(1.5).*...
                tau.^((2.*(nu-1)+3)./(2.*(nu-1)+4)).*(2.*pi.*(nu+1)/nu).^0.5).*...
                (1-((2.*(nu-1)+3).*(nu+2))./(24.*x_new.*tau.^(1/(nu+1)).*nu.*(nu+1)));


    
    
end

%end

display(k)

elseif kswitch == 1

T = 0.0429.*t;

gam_series_old = zeros(1,length(t));

%gam_series_old = zeros(length(x),length(t));
%gam_series = 100;
gam_tol = 1;
gam_tol_max = 1e-5;
k = 0;

x_new = x;

if x_new<xt

while gam_tol > gam_tol_max
%while k <100 

gam_series  = gam_series_old +...
((x_new.*nu).^(nu.*(k+1))./gamma(nu.*(k+1))).*(T./(T+2)).^k;
    
%gam_series  = vpa(gam_series_old +...
%tau.^(k) .* x.^(k.*(nu+1)) .*nu.^(nu.*(k+1))./...
%(gammaz(k+2).*gammaz(nu.*(k+1))));

gam_tol = ...
max(abs(gam_series_old-gam_series)./gam_series_old,[],'omitnan');

gam_series_old = gam_series;

k = k + 1;
end

n_xt = ((4.*exp(-nu.*x_new))./(x_new.*(T+2).^2)) .*gam_series_old;


%xt_new = x_new(x_new>=xt);


    
elseif x_new>=xt  
    
    ys = zeros(1,length(T));
    
for tt = 1 : length(T)    
    
ys_eq = @(ys) (T(tt)./(T(tt)+2))-ys.^(nu-1).*(ys-1./x_new);

ys(tt) = fsolve(ys_eq,1); 

end

% Not sure where the 0.9221 factor comes from but it seems to be necessary.


n_xt = 0.9221.*4.*(nu).*exp((ys-1).*(nu).*x_new)./((T+2).^2.*(ys.^(nu-1)).*(2.*pi.*nu.*(nu-(nu-1)./(ys.*x_new))).^(0.5));




end

display(k)   
    

% Product kernel
elseif kswitch == 2 
   
%t = t.*(7./40); % correction factor?
    
%T = 0.00959.*t;

T = (7/40).*0.00959.*t;

gam_series_old = zeros(1,length(t));

%gam_series_old = zeros(length(x),length(t));
%gam_series = 100;
gam_tol = 1;
gam_tol_max = 1e-5;
k = 0;

x_new = x;

while gam_tol > gam_tol_max
%while k <100 

gam_series  = gam_series_old +...
((x_new).^((nu-1)+k.*(nu+2)).*nu.^((nu+1).*(k+1))./(factorial(k+1).*gammaz((nu+1).*(k+1)))).*(T).^k;
    
%gam_series  = vpa(gam_series_old +...
%tau.^(k) .* x.^(k.*(nu+1)) .*nu.^(nu.*(k+1))./...
%(gammaz(k+2).*gammaz(nu.*(k+1))));

gam_tol = ...
max(abs(gam_series_old-gam_series)./gam_series_old,[],'omitnan');

gam_series_old = gam_series;

k = k + 1;
end


n_xt = (exp(-x_new.*(T+nu))) .*gam_series_old;

%xt_new = x_new(x_new>=xt);

if (x_new>=xt)
n_xt = ((nu+1).*exp(-x_new.*(T+nu)+x_new.*(nu+2).*T.^(1./(nu+2)).*(nu./(nu+1)).^((nu+1)./(nu+2)))./...
    (x_new.^(5./2).*(T.*(nu+1)./nu).^((2.*(nu+2)-1)./(2.*(nu+2))).*(2.*pi.*(nu.*(nu+2))).^(1./2))).*...
    (1-((nu+3).*(2.*(nu+2)-1))./(24.*x_new.*T.^(1./(nu+2)).*nu.^((nu+1)./(nu+2)).*(nu+1).^(1./(nu+2)).*(nu+2)));

end


display(k)   

% SBE solution following Feingold et al. 1987
elseif kswitch == 3

n0_xt = (1./gamma(nu)).*(nu.^nu).*x.^(nu-1).*exp(-nu.*x); % initial dist.
n_xt = (n0_xt+gam.*(exp(B.*gam.*t)-1).*exp(-gam.*x))./(1+(1./gam).*(exp(B.*gam.*t)-1));

% Steady state SCE/SBE solution following Feingold et al. 1987
elseif kswitch == 4
    
    
    
    %E = 1-B;
    
    E = 1 - t;
    
    eta = 0.5.*E.*gam.*(2-E);
    
    n_xt = (1-E).*gam.^2.*(besseli(0,eta.*x)-besseli(1,eta.*x)).*exp(-(gam-eta).*x);
    
    
    % PUT IN APPROX. FOR LARGE ARGUMENT.
    
 %   if x > 100
        
%        I0_large = (1./(sqrt(2.*pi))).*eta.*x.*exp((eta.*x));
        
 %       I1_large = (1./(sqrt(2.*pi.*eta.*x).*(1./(eta.*x).^2).^(1./4))).*exp(-asinh(1./(eta.*x))+(eta.*x).*sqrt(1+(1./(eta.*x)).^2));
        
 %       n_xt = (1-E).*gam.^2.*(I0_large-I1_large).*exp(-(gam-eta).*x);

        
        
 %   end
    
    
    
end

end