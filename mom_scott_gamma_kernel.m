function M_rt = mom_scott_gamma_kernel(r,t,nu,B,gam,kswitch)
% Scott's 1968 analytical solution to 
% collection equation for Gamma (Gaussian-like)
% initial distribution and Golovin kernel.

syms k
M_rt = zeros(1,length(t));
if kswitch == 0
tau = 1 - exp(-0.00153.*t); % normalized time

for tt = 1:length(t)
    
gam_series = symsum(tau(tt).^(k) .*nu.^(nu.*(k+1))./...
((tau(tt)+nu).^(r+nu+k.*(nu+1)) .* factorial(k+1)) .*...
(gamma(nu+r+k.*(nu+1))./...
 gamma(nu.*(k+1))),0,inf);

M_rt(tt) = (1-tau(tt)) .*gam_series;

end


elseif kswitch == 1
    
T =  0.0429.*t;
    
for tt = 1:length(t)
gam_series = symsum(gamma(r+nu.*(k+1)).*((T(tt)./(T(tt)+2)).^k)./gamma(nu.*(k+1)),0,inf);
 
M_rt(tt) = (4./(T(tt)+2).^2).*nu^(-r).*gam_series;

end
    
    

elseif kswitch == 2
    
%T =  0.00959.*t;

T =  (7/40).*0.00959.*t; % NOTE: I have no idea where this 7/40 correction factor comes from but apparently it's needed in order to match Scott's results.

%r = r-1;
    
for tt = 1:length(t)
    
gam_series = symsum((T(tt).^k).*nu.^((k+1).*(nu+1)).*(T(tt)+nu).^(-(nu+r+k.*(nu+2))).*gamma(r+nu+k.*(nu+2))./(factorial(k+1).*gamma((k+1).*(nu+1))),0,inf);
 
M_rt(tt) = gam_series;

end
  


elseif kswitch == 3

 
 M0r = nu.^(-r) .*gamma(nu+r)./gamma(nu);
    
 
 for tt = 1 :length(t)
 den = (1+(1./gam).*(exp(B.*gam.*t(tt))-1));   
M_rt(tt) = (M0r+gam.*(exp(B.*gam.*t(tt))-1).*gam.^(-1-r).*gamma(1+r))./den;
 end
 
 
 
 
 elseif kswitch == 4
     
     
     % SBE/SCE state state solution (with the help of Mathematica)
     
     %E = 1-B;
     
     E = 1 - t;
     
     
     eta = 0.5.*E.*gam.*(2-E);
     
    M_rt = 0.5.*(1-E).*gam.^2.*gamma(r+1).*(gam-eta).^(-r-2).*...
        (2.*(gam-eta).*hypergeom([0.5.*(r+1),0.5.*(r+2)],1,(eta./(gam-eta)).^2)-...
         eta.*(r+1).*sign(gam-eta).*hypergeom([0.5.*(r+2) 0.5.*(r+3)],2,(eta./(gam-eta)).^2));
 
 
 
 
 
    
end




 


end