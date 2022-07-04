% Scott_test_PSD_fitting.m
% Estimate various PSD parameters that represent
% Scott 1968 and Graham Feingold's 1987 analytical solution to SCE
% (SCE/SBE)


tic

opts = optimoptions(@lsqnonlin,'Display','iter-detailed','MaxFunctionEvaluations',10000,'MaxIterations',10000,...
                                                              'SpecifyObjectiveGradient',false,'CheckGradients',false,...
                                                             'functiontolerance',1e-16,'steptolerance',1e-16,...
                                                             'optimalitytolerance',1e-16);

mom_switch = 0;
                                                         
                                                         
ncase = 'GIG2';

Kernel_switch = 2;
x = logspace(-4,3,500);

xlow = -2;
xhigh = 4;


if Kernel_switch == 4
    
xlow = -4;
xhigh = 3;

tsind = 1;



elseif Kernel_switch == 2
    
    xlow = -2;
    xhigh = 3;
    
   % tsind = 50;
    
    
    tsind = 1;
    
    
elseif Kernel_switch == 0
    
    xlow = -2;
    xhigh = 4;
    
    %tsind = 50;
    
    tsind = 1;
    
    x = logspace(xlow,xhigh,500);
    
    
elseif Kernel_switch == 1
    
    xlow = -2;
    xhigh = 3;
    
    tsind = 1;
    
end

nu = 4;
B = 0.9;
C = 0.1;

gam = 10;

%t = 0:10:2000;

t = [0 50 100 200 400];

%t = 0:5:400;

%t = logspace(1,log10(400),30);

%t = [0 400 800 1600];

%t = [0 500 1000 1500 2000];

%t = [0 10 50 100 800];

if Kernel_switch == 4
    
    nu = 1;
    
    %t = 0;
    
    %t = 0.01:0.01:0.95;
    
    t = [0.1 0.3 0.5 0.7 0.9];
    
    % t = b;
    
end

spacing = int32(((t/max(t))).*255+1);

spacing(t==0)=1;
         
cmap = flipud(cmap_blueteal);

%newcolors = cmap(spacing,:);       
         

newcolors = [0 0 0; 0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19; 0.49 0.18 0.56];



if Kernel_switch ==0
tau = 1 - exp(-0.00153.*t); % normalized time
elseif Kernel_switch == 1
    T = 0.00959.*t;
end


N = 8;

r = 0:1:N-1;

%r = 0:0.1:N-1;

mom3 = find(r==3);
mom4 =  find(r==4);

if mom_switch == 1

% Scott's moments
M_rt = NaN(N,length(t));
M_rt2 = NaN(N,length(t));


for rdum = 1 : N
    
M_rt(rdum,:) = mom_scott_gamma_kernel(r(rdum),t,nu,B,gam,Kernel_switch);

fprintf('r =%i\n',r(rdum));

end

end


%% Calculate synthetic PSD moment "errors" 

rel_error = 0.01;

M_rt_min = (1-rel_error).*M_rt;
M_rt_max = (1+rel_error).*M_rt;

cost_min = nansum(rootnd(M_rt_min,M_rt));
cost_max = nansum(rootnd(M_rt_max,M_rt));

cost_error = 0.5.*(cost_min+cost_max); % Total "error"



% Scott's solution

if mom_switch == 1
n_xt = NaN(length(x),length(t));

for xtemp = 1:length(x)

n_xt(xtemp,:) = n_scott_gamma_kernel(x(xtemp),t,nu,B,gam,Kernel_switch);

end

n_xt(n_xt<0)=0;

end

% mean volume diameter
Dmv = M_rt(5,:)./M_rt(4,:);




% Normalize by number
M_rt_norm = M_rt./M_rt(1,:);

% My solution using the H-function

H_xt     = NaN(length(x),length(t));
H_mxt = NaN(length(x),length(t));

Mon = M_rt_norm';

Mrn = NaN(length(r),length(t));


%% Optimize

for time = 1 : length(t)
    
    disp(time)
    
    
    [Mn,nt,lb,lu,x0,xn,pnames] = choose_case(Dmv(time),r,x,ncase);
    
    if time == 1
        
        params = NaN(length(x0),length(t));
        resnorm = NaN(1,length(t));
        residual = NaN(N,length(t));

        
    end
    
     
     cost = @(z) rootnd(Mon(time,:),Mn(z));
     
     
     [param_temp,resnorm(time),residual(:,time),~] = lsqnonlin(@(z) cost(z),x0,lb,lu,opts);
    
     
     
     
     %H_xt(:,time) = M_rt(1,time).*nt(param_temp);
     
     H_xt(:,time) = nt(param_temp);

     params(:,time) = param_temp;
     
     
     Mrn(:,time) = Mn(param_temp);

end


%[A,c] = ellip_space(Pspace,cost,cost_error);

%figure;
%subplot(1,2,1)
%plot(x,n_xt);
%hold on;
%subplot(1,2,2);
%plot(x,H_xt,'--');
%xlim([0 4.8])

%set(gca,'xscale','log')


%% Plotting
hfig = figure;
set(hfig,'units','normalized','position',[0.5076    0.0600    0.2681    0.8333])
subplot(3,1,1);
hold on;
for ii = 1 : tsind:length(t)

if Kernel_switch ~=4
plot(log10(x),x' .*n_xt(:,ii)./M_rt(1,ii),'color',newcolors(ii,:),'linewidth',1);

%plot(log10(x),x' .*H_xt(:,ii)./M_rt(1,ii),'--','color',newcolors(ii,:),'linewidth',1);

plot(log10(x),x' .*H_xt(:,ii)./Mrn(1,ii),'--','color',newcolors(ii,:),'linewidth',1);


else
    
    % For SCE/SBE steady state test, normalize by fragment count
plot(log10(x),x' .*n_xt(:,ii)./M_rt(1,ii),'color',newcolors(ii,:),'linewidth',1);



%plot(log10(x),x' .*H_xt(:,ii)./M_rt(1,ii),'--','color',newcolors(ii,:),'linewidth',1);



plot(log10(x),x' .*H_xt(:,ii)./Mrn(1,ii),'--','color',newcolors(ii,:),'linewidth',1);
    
    
end


end
set(gca,'xtick',[-3 -2 -1 0 1 2],'xticklabel',...
   {'0.001' '0.01' '0.1' '1.0' '10.0' '100.0'})
xlim([xlow xhigh]);
%ylim([0 1])


subplot(3,1,2);
hold on;
for ii = 1 : tsind:length(t)

if Kernel_switch~=4
plot(log10(x),x'.^2 .*n_xt(:,ii)./M_rt(2,ii),'color',newcolors(ii,:),'linewidth',1);



%plot(log10(x),x'.^2 .*H_xt(:,ii)./M_rt(2,ii),'--','color',newcolors(ii,:),'linewidth',1);

plot(log10(x),x'.^2 .*H_xt(:,ii)./Mrn(2,ii),'--','color',newcolors(ii,:),'linewidth',1);



else
    
plot(log10(x),x'.^2 .*n_xt(:,ii)./M_rt(2,ii).*gam,'color',newcolors(ii,:),'linewidth',1);



%plot(log10(x),x'.^2 .*H_xt(:,ii)./M_rt(2,ii).*gam,'--','color',newcolors(ii,:),'linewidth',1);



plot(log10(x),x'.^2 .*H_xt(:,ii)./Mrn(2,ii).*gam,'--','color',newcolors(ii,:),'linewidth',1);

    
end
end
set(gca,'xtick',[-3 -2 -1 0 1 2],'xticklabel',...
   {'0.001' '0.01' '0.1' '1.0' '10.0' '100.0'})
xlim([xlow xhigh]);
%ylim([0 0.7])



subplot(3,1,3);
hold on;
for ii = 1 : tsind:length(t)

if Kernel_switch~=4
plot(log10(x),x'.^3 .*n_xt(:,ii)./M_rt(3,ii),'color',newcolors(ii,:),'linewidth',1);



%plot(log10(x),x'.^3 .*H_xt(:,ii)./M_rt(3,ii),'--','color',newcolors(ii,:),'linewidth',1);



plot(log10(x),x'.^3 .*H_xt(:,ii)./Mrn(3,ii),'--','color',newcolors(ii,:),'linewidth',1);

else
    
plot(log10(x),x'.^3 .*n_xt(:,ii)./M_rt(3,ii).*gam,'color',newcolors(ii,:),'linewidth',1);


%plot(log10(x),x'.^3 .*H_xt(:,ii)./M_rt(3,ii).*gam,'--','color',newcolors(ii,:),'linewidth',1);


plot(log10(x),x'.^3 .*H_xt(:,ii)./Mrn(3,ii).*gam,'--','color',newcolors(ii,:),'linewidth',1);

    
end
end
set(gca,'xtick',[-3 -2 -1 0 1 2],'xticklabel',...
   {'0.001' '0.01' '0.1' '1.0' '10.0' '100.0'})
xlim([xlow xhigh]);
%ylim([0 0.7])

toc

function [Mn,nx,lb,lu,x0,xn,pnames] = choose_case(Dmv,n,x,ncase)

switch ncase
    
    case 'gamma'
        
            %plen = 3;
   
            
            pnames = {'nu' 'b'};
            
            xn = @(z) Dmv.*gamma(3./(z(2))+(z(1)))./gamma(4./(z(2))+(z(1)));
            
            lb = [0.1 0.1];
            lu = [10 10];
           x0 = [8 1];
           
           
           Mn = @(z) xn(z).^n.*gamma(n./(z(2))+z(1))./gamma(z(1));
           
           nx =  @(z) 1./(gamma(z(1))).*z(2).*(1./xn(z)).*(x./xn(z)).^(z(1).*z(2)-1).*exp(-(x./xn(z)).^(z(2)));
           
           
          
%             nu = linspace(lb(1),lu(1),50);
%             b = linspace(lb(2),lu(2),50);
%             
%             [P1,P2] = meshgrid(nu,b);
%             
%             
%             Pspace = [P1(:) P2(:)]';
           
           %logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2));
           
           
           %stepsize = 5;
           
           
       case 'gamma_inc'
        
            %plen = 3;
   
            
            pnames = {'nu' 'b' 'm0' 'mn'};
            
            %xn = @(z) Dmv.*gamma(3./(z(2))+(z(1)))./gamma(4./(z(2))+(z(1)));
            
            lb = [-10 0.1 1e-8 0.001];
            lu = [10  10 0.001 1000];
           x0 = [8 1 1e-5 10];
           
           
          % Mn = @(z) xn(z).^n.*gamma(n./(z(2))+z(1))./gamma(z(1));
           
           Mn = @(z) z(4).^n .* (gamma(z(1)+n./z(2)))./gamma(z(1)).*gammainc(z(1)+(n./z(2)),z(3)./z(4))./gammainc(z(1),z(3)./z(4));
           
           nx =  @(z) 1./(gamma(z(1)).*gammainc(z(1),z(3)./z(4))).*z(2).*(1./z(4)).*(x./z(4)).^(z(1).*z(2)-1).*exp(-(x./z(4)).^z(2));
           
           xn = @(z) z(4);
          
%             nu = linspace(lb(1),lu(1),50);
%             b = linspace(lb(2),lu(2),50);
%             
%             [P1,P2] = meshgrid(nu,b);
%             
%             
%             Pspace = [P1(:) P2(:)]';
           
           %logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2));
           
           
           %stepsize = 5;
                  
           
           
           
           
           
    case 'exp2'
        
        
            pnames = {'f' 'fn'};
        
        
        
            xn = @(z) (Dmv./4).*(z(1)+(1-z(1)).*z(2).^(-3))./(z(1)+(1-z(1)).*z(2).^(-4));

        
            lb = [0 0.1];
            lu = [1 1];
            x0 = [0.5 1];
     
        
          % logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2));
           
          Mn = @(z) xn(z).^n.*(z(1).*gamma(n+1)+(z(1)+(1-z(1)).*z(2).^(-n)));
          
           nx = @(z) (z(1)./xn(z)).*exp(-x./xn(z))+((1-(z(1)))./xn(z)).*z(2).*exp(-z(2).*x./xn(z));
        
           %stepsize = 5;
           
       case 'gamexp'
        
            pnames = {'f' 'fn' 'nu'};
        
       
            xn = @(z) Dmv.*((z(1).*gamma(z(3)+3)./gamma(z(3)))+((1-z(1)).*(z(3)./(z(2))).^3 .* 6))./(((z(1).*gamma(z(3)+4)./gamma(z(3))))+((1-z(1)).*(z(3)./(z(2))).^4 .* 24));

        
            lb = [0 0.1 0];
            lu = [1 1 10];
            x0 = [0.5 1 1];
     
        
           Mn = @(z) ((z(1).*gamma(z(3)+n)./gamma(z(3)))+((1-z(1)).*(z(3)./(z(2))).^n .* gamma(1+n)));
                   
           
           nx =  @(z) (1./gamma(z(3))).*(z(1)./xn(z)).*(x./xn(z)).^(z(3)-1).*exp(-x./xn(z))+...
                     ((1-(z(1)))./xn(z)).*(z(2)./z(3)).*exp(-((z(2))./z(3)).*x./xn(z));
           
      case 'expgam'
        
        
            pnames = {'f' 'fn' 'nu'};
        
       
            xn = @(z) Dmv.*((z(1).*6)+((1-z(1)).*(1./(z(3).*z(2))).^3 .* gamma(z(3)+3)./gamma(z(3))))./(((z(1).*24))+((1-z(1)).*(1./(z(3).*z(2))).^4 .* gamma(z(3)+4)./gamma(z(3))));

        
            lb = [0 0.1 0];
            lu = [1 100 10];
            x0 = [0.5 1 1];
     
        
            %logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2))&& (x(3)>=lb(3))&&(x(3)<=lu(3));    
           
           Mn = @(z) ((z(1).*gamma(1+n)+((1-z(1)).*(1./(z(3).*z(2))).^n .* gamma(z(3)+n)./gamma(z(3)))));
                   
           
           nx =  @(z) (z(1)./xn(z)).*exp(-x./xn(z))+...
                       (1./gamma(z(3))).*(x./xn(z)).^(z(3)-1).*((1-(z(1)))./xn(z)).*(z(2).*z(3)).^z(3).*exp(-((z(2).*z(3))).*x./xn(z));
           
                                 
                   
           %stepsize = 4;      

     case 'gam2'
        
        
            pnames = {'f' 'fn' 'nu1', 'nu2'};
        
       
            xn = @(z) Dmv.*((z(1).*gamma(z(3)+3)./gamma(z(3)))+((1-z(1)).*(z(3)./(z(4).*z(2))).^3 .* gamma(z(4)+3)./gamma(z(4))))./...
                                       (((z(1).*gamma(z(3)+4)./gamma(z(3))))+((1-z(1)).*(z(3)./(z(4).*z(2))).^4 .* gamma(z(4)+4)./gamma(z(4))));

        
            lb = [0 0 0 0];
            lu = [1 10 10 10];
            x0 = [0.5 0.4 2 3];
     
            
            
            %lb = [0 0 0 0];
            %lu = [1 10 30 30];
            %x0 = [0.5 0.4 2 3];
        
            %logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2))&& (x(3)>=lb(3))&&(x(3)<=lu(3))&&...
            %                    (x(4)>=lb(4))&&(x(4)<=lu(4));    
           
           %Mn = @(n) (1./gamma(z(3))).*(z(1)./xn(z)).*(D_temp./xn(z)).^(z(3)-1).*exp(-D_temp./xn(z))+...
           %            (1./gamma(z(4))).*((1-(z(1)))./xn(z)).*(z(2).*z(4)./z(3)).^z(4).*exp(-((z(2).*z(4))./z(3)).*D_temp./xn(z));
      
           Mn = @(z) xn(z).^n.*((z(1).*gamma(z(3)+n)./gamma(z(3)))+((1-z(1)).*(z(3)./(z(4).*z(2))).^n .* gamma(z(4)+n)./gamma(z(4))));
                   
           
           nx =  @(z) (1./gamma(z(3))).*(z(1)./xn(z)).*(x./xn(z)).^(z(3)-1).*exp(-x./xn(z))+...
                       (1./gamma(z(4))).*(x./xn(z)).^(z(4)-1).*((1-(z(1)))./xn(z)).*(z(2).*z(4)./z(3)).^z(4).*exp(-((z(2).*z(4))./z(3)).*x./xn(z));
           
           %stepsize   = 2;     
                   
    case 'GIG2'
        
        
           pnames = {'nu' 'fn'};
        
        
        
           % xn = @(z) (Dmv.*besselk(3+(z(1)),2.*z(2).^(0.5))./besselk(4+(z(1)),2.*z(2).^(0.5))).^2;

            xn = @(z) (Dmv.*besselk(z(1)+3,2.*z(2))./besselk(z(1)+4,2.*z(2)));
        
        
            lb = [-10 0];
            lu = [10 10];
            x0 = [3 3];  
     
        
          % logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2));                        
                     
           %Mn = @(z) xn(z).^(n).*besselk(z(1)+n,2.*sqrt(z(2)))./besselk(z(1),2.*sqrt(z(2)));
           
           Mn = @(z) xn(z).^(n).*besselk(z(1)+n,2.*z(2))./besselk(z(1),2.*z(2));
           
           
            
                   
          % stepsize = 5;
          
          
          %nx = @(z) 1./(2.*xn(z).^z(1).*besselk(z(1),2.*sqrt(z(2)))).*x.^(z(1)-1).*exp(-(sqrt(z(2))./xn(z)).*(x+xn(z).^2./x));

          nx = @(z) 1./(2.*besselk(z(1),2.*z(2))).*(1./xn(z)).*(x./xn(z)).^(z(1)-1).*exp(-(z(2)).*((x./xn(z))+(xn(z)./x)));

          
          
          
     case 'genGIG3'
           
           
           %mccount = 2000000;
        
           nutemp = @(z) sign(z(1)).*abs(z(1)).^(1./z(3));
           
           fntemp = @(z) z(2).^(1./z(3));
        
           pnames = {'nu' 'fn' 'b'};
           
            %lb = [-10 0. 0.1];
            %lu = [10 10 5];
           % x0 = [3 3 1];  
            
            
            lb = [-15 0. 0.1];
            lu = [15 10 5.0];
            x0 = [3 3 1];  
        
            
            xn = @(z) Dmv.*besselk(nutemp(z)+3./z(3),2.*fntemp(z))./besselk(nutemp(z)+4./z(3),2.*fntemp(z));
              
            Mn = @(z) xn(z).^(n).*besselk(nutemp(z)+n./z(3),2.*fntemp(z))./besselk(nutemp(z),2.*fntemp(z));
           
            nx = @(z) z(3)./(2.*besselk(nutemp(z),2.*fntemp(z))).*(1./xn(z)).*(x./xn(z)).^(nutemp(z).*z(3)-1).*exp(-(fntemp(z)).*((x./xn(z)).^z(3)+(xn(z)./x).^z(3)));

           
           
              
           
        
        
    case 'GIG3'
        
        
           pnames = {'nu' 'fn' 'b'};
        
           %% DOUBLE CHECK!
            % sqrt(sn)
            
            % OLD
            %xn = @(z) (Dmv.*besselk(3+(z(1)./z(3)),2.*z(2).^(0.5))./besselk(4+(z(1)./z(3)),2.*z(2).^(0.5)));

            % NEW?
            %xn = @(z) (Dmv.*besselk(z(1)+3./z(3),2.*z(2).^(0.5))./besselk(z(1)+4./z(3),2.*z(2).^(0.5)));
               
            xn = @(z) (Dmv.*besselk(z(1)+3./z(3),2.*z(2))./besselk(z(1)+4./z(3),2.*z(2)));
        
            lb = [-10 0.00001 0.1];
            lu = [10 10 5];
            x0 = [3 3 1];  
     
        
           %logprior = @(x) (x(1)>=lb(1))&&(x(1)<=lu(1))&& (x(2)>=lb(2))&&(x(2)<=lu(2))&& (x(3)>=lb(3))&&(x(3)<=lu(3));                        
                     
           %Px = @(z) D_temp.^(z(1).*z(3)-1).*exp(-((z(2)./xn(z)).^(0.5).*(D_temp.^z(3)+xn(z)./D_temp.^(z(3)))));

           
          % Normalized Mn: Mn/M0
           %Mn = @(z) xn(z).^(n).*besselk(z(1)+n./z(3),2.*sqrt(z(2)))./besselk(z(1),2.*sqrt(z(2)));
           
           Mn = @(z) xn(z).^(n).*besselk(z(1)+n./z(3),2.*z(2))./besselk(z(1),2.*z(2));
           
           
           %nx = @(z) z(3)./(2.*xn(z).^z(1).*besselk(z(1),2.*sqrt(z(2)))).*x.^(z(1).*z(3)-1).*exp(-(sqrt(z(2))./xn(z)).*(x.^z(3)+xn(z).^2./x.^z(3)));

           
           %nx = @(z) z(3)./(2.*xn(z).^z(1).*besselk(z(1),2.*z(2))).*x.^(z(1).*z(3)-1).*exp(-(z(2)./xn(z)).*(x.^z(3)+xn(z).^2./x.^z(3)));

            nx = @(z) z(3)./(2.*besselk(z(1),2.*z(2))).*(1./xn(z)).*(x./xn(z)).^(z(1).*z(3)-1).*exp(-(z(2)).*((x./xn(z)).^z(3)+(xn(z)./x).^z(3)));

           
           %stepsize = 3;
           
 
          
           
end




end




function cost = rootnd(Mon,Mn)


cost = (Mon-Mn).^2./max(Mon.*Mn,1e-60);

end


