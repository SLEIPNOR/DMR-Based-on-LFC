% Performace Setting
disp('********************************************************')
disp('*                                                      *')
disp('*                  Designed Result                     *')
disp('*                                                      *')
disp('********************************************************')
desired_precision = 50;

fprintf('\n')
disp('Performance Requirements:')

N_p = 1.04;  
AoV = pi/180*30;
Rgmax = 10;
Rgmin = 7;
P = 10^-1.4;
PR = P/N_p;
fprintf('AoV = %.10f\n',AoV)
fprintf('Rgmax = %.10f\n',Rgmax)
fprintf('Rgmin = %.10f\n',Rgmin)
fprintf('P = %.10f\n\n',P)

PS = 0.56e-6;


beta = vpa(1/Rgmax +1/Rgmin,desired_precision);

gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);

delta = vpa(Rgmax-PR,desired_precision);

PIP = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision)

P_MLA = vpa(PS*PIP,desired_precision)

F = vpa(((P*PIP)/(Rgmax^2*gamma)-1)*(2/beta),desired_precision)
D_M = vpa((2*PIP*P_MLA)/(F*gamma),desired_precision)
f_M = vpa(2*F/(beta*F+2),desired_precision)
D_S = 2*F*tan(AoV/2)
f_MLA = (f_M/D_M)*P_MLA
%%
desired_precision = 50;
PS = 0.56e-6;
F = 53e-3;
P_MLA = 13.89e-6;
D_M = 25e-3;
f_M = 50e-3;
D_S = 6.14e-3;
PIP = 10;
f_sst = 1/P_MLA;
f_suv = PIP/D_M;

Rgmax = vpa((f_M*f_sst*F)/(f_sst*F-f_M*(f_sst+f_suv)),desired_precision)
Rgmin = vpa((f_M*f_sst*F)/(f_sst*F-f_M*(f_sst-f_suv)),desired_precision)

AN = vpa(f_M/D_M,desired_precision);
Rmax = Rgmax^2/(Rgmax-f_M);
DoF = vpa((2*f_M^2*AN*P_MLA*Rmax^2)/(f_M^4),desired_precision)


imdist = vpa(F*(f_sst*D_M)/((f_suv+f_sst)*D_M-2),desired_precision)
O = vpa(f_M*imdist/(imdist-f_M),desired_precision);
PR = vpa(Rgmax-O,desired_precision)
P = max(DoF,PR)
AoV = 180*2*atan(D_S/(2*F))/pi
%% mu_p

fprintf('\n')
disp('Performance Requirements:')
desired_precision = 50;
AoV = pi/180*20;
Rgmax = 2;
Rgmin = 0.5;
P = 0.1;
PS = 0.56e-6;

i = 1;
min = 1.1;
max = 1.4;
num =100;
for N_p = min:(max - min)/(num-1) :max
    PR = P/N_p;

    beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
    
    gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
    
    delta = vpa(Rgmax-PR,desired_precision);
    
    PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
    
    P_MLA(i) = vpa(PS* PIP(i),desired_precision);
    
    F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
    D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
    f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
    D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
    f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);

    i = i+1
end


% %% max depth
% 
% fprintf('\n')
% disp('Performance Requirements:')
% desired_precision = 50;
% AoV = pi/180*30;
% P = 0.05;
% N_p = 1.06;
% PR = P/N_p;
% 
% 
% PS = 0.56e-6;
% 
% i = 1;
% min = 3;
% max = 20;
% num =200;
% for Rgmax = min:(max - min)/(num-1) :max
% 
%    
%     Rgmin = Rgmax-2;
%  
% 
%     beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
%     
%     gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
%     
%     delta = vpa(Rgmax-PR,desired_precision);
%     
%     PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
%     
%     P_MLA(i) = vpa(PS* PIP(i),desired_precision);
%     
%     F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
%     D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
%     f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
%     D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
%     f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);
% 
%     i = i+1
% end
% 
% %% Range
% 
% fprintf('\n')
% disp('Performance Requirements:')
% desired_precision = 50;
% AoV = pi/180*30;
% P = 0.05;
% N_p = 1.06;
% PR = P/N_p;
% Rgmax = 10;
% 
% PS = 0.56e-6;
% 
% i = 1;
% min = 0.1;
% max = 5;
% num =200;
% for Rg = min:(max - min)/(num-1) :max
% 
%    
%     Rgmin = Rgmax-Rg;
%  
% 
%     beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
%     
%     gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
%     
%     delta = vpa(Rgmax-PR,desired_precision);
%     
%     PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
%     
%     P_MLA(i) = vpa(PS* PIP(i),desired_precision);
%     
%     F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
%     D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
%     f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
%     D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
%     f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);
% 
%     i = i+1
% end
% 
% %% Pm
% 
% fprintf('\n')
% disp('Performance Requirements:')
% desired_precision = 50;
% AoV = pi/180*30;
% N_p = 1.06;
% Rgmax = 10;
% Rgmin = 8;
% PS = 0.56e-6;
% 
% 
% min = -1.5;
% max = -1;
% num =200;
% values = logspace(min, max, num);
% 
% for i = 1:num
%     P = values(i);
%     PR = P/N_p;
% 
%     beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
%     
%     gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
%     
%     delta = vpa(Rgmax-PR,desired_precision);
%     
%     PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
%     
%     P_MLA(i) = vpa(PS* PIP(i),desired_precision);
%     
%     F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
%     D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
%     f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
%     D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
%     f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);
% 
%     i
% end
% 
% %% AoV
% 
% fprintf('\n')
% disp('Performance Requirements:')
% desired_precision = 50;
% AoV = pi/180*30;
% Rgmax = 10;
% Rgmin = 8;
% P = 0.05;
% PS = 0.56e-6;
% N_p = 1.06;
% PR = P/N_p;
% 
% 
% i = 1;
% min = 1;
% max = 90;
% num =200;
% for Av = min:(max - min)/(num-1) :max
%    
%     AoV = pi/180*Av;
%     beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
%     
%     gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
%     
%     delta = vpa(Rgmax-PR,desired_precision);
%     
%     PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
%     
%     P_MLA(i) = vpa(PS* PIP(i),desired_precision);
%     
%     F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
%     D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
%     f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
%     D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
%     f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);
% 
%     i = i+1
% end
% 
% %% PS 
% 
% fprintf('\n')
% disp('Performance Requirements:')
% desired_precision = 50;
% AoV = pi/180*30;
% Rgmax = 10;
% Rgmin = 8;
% P = 0.05;
% 
% N_p = 1.06;
% PR = P/N_p;
% AoV = pi/180*30;
% 
% i = 1;
% min = 0.01e-6;
% max = 2e-6;
% num =200;
% for PS = min:(max - min)/(num-1) :max
%    
%     
%     beta = vpa(1/Rgmax +1/Rgmin,desired_precision);
%     
%     gamma = vpa(1/Rgmin-1/Rgmax,desired_precision);
%     
%     delta = vpa(Rgmax-PR,desired_precision);
%     
%     PIP(i) = vpa((2*delta*gamma)/(delta*gamma-delta*beta+2),desired_precision);
%     
%     P_MLA(i) = vpa(PS* PIP(i),desired_precision);
%     
%     F(i) = vpa(((P* PIP(i))/(Rgmax^2*gamma)-1)*(2/beta),desired_precision);
%     D_M(i) = vpa((2* PIP(i)*P_MLA(i))/(F(i)*gamma),desired_precision);
%     f_M(i) = vpa(2*F(i)/(beta*F(i)+2),desired_precision);
%     D_S(i) = vpa(2*F(i)*tan(AoV/2),desired_precision);
%     f_MLA(i) = vpa((f_M(i)/D_M(i))*P_MLA(i),desired_precision);
% 
%     i = i+1
% end
% 

%% plot 1

figure(1)

N_p_T = min:(max - min)/(num-1) :max;
% Rgmax_T = min:(max - min)/(num-1) :max;
% Rg_T = min:(max - min)/(num-1) :max;
% P_T = logspace(min, max, num);
% AoV_T = min:(max - min)/(num-1) :max;
% PS_T= min:(max - min)/(num-1) :max;
% PS_T = PS_T *10^6;
% semilogx(P_T,F*10^3,'LineWidth', 2);
plot(N_p_T ,F*10^3,'LineWidth', 2);
hold on 
% semilogx(P_T,f_M*10^3,'LineWidth', 2);
plot(N_p_T ,f_M*10^3,'LineWidth', 2);
hold on 
% semilogx(P_T,D_M*10^3,'LineWidth', 2);
plot(N_p_T ,D_M*10^3,'LineWidth', 2);
hold on 
% semilogx(P_T,D_S*10^3,'LineWidth', 2);
plot(N_p_T ,D_S*10^3,'LineWidth', 2);
grid on;

pbaspect([1 1 1]);

xlabel('\mu_p');
% xlabel('{\it Rg}_{max} (m)')
% xlabel('{\it A}_{V} (°)')
% xlabel('{\it S}_{px} (µm)')
% xlabel('{\it Rg}_{max}-{\it Rg}_{min} (m)')
% xlabel('{\it P}_{m} (m)')
ylabel('(mm)');
% xticks(Rgmax_T(1):2 :Rgmax_T(end));
% yticks(floor(F(1)):1:ceil(F(end)));

legend('\it F', '{\it f}_M', '{\it D}_M', '{\it D}_S','Location', 'Best','FontSize', 12);

% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; {\it S}_{px} = 0.56µm" )
% title("{\it Rg}_{max}-{\it Rg}_{min} = 2m; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("{\it Rg}_{max} = 10m; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm" )
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; \mu_p = 1.06; {\it A}_{V} = 30°")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [0.5m,2m]; \n {\it P}_{m} = 0.1m; {\it A}_{V} = 20°; {\it S}_{px} = 0.56µm" )
title("LFC 1" )
set(gca, 'FontSize', 20);


%% plot 2
yyaxis left;
% semilogx(P_T ,P_MLA*10^6,'LineWidth', 2);
plot(N_p_T ,P_MLA*10^6,'LineWidth', 2);
ylabel('{\it P}_{MLA}(µm)');

% 添加第二个 Y 轴
yyaxis right;
% semilogx(P_T ,f_MLA*10^6,'LineWidth', 2);
plot(N_p_T ,f_MLA*10^6,'LineWidth', 2);

ylabel('{\it f}_{MLA}(µm)');
grid on;

% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; {\it S}_{px} = 0.56µm" )
% title("{\it Rg}_{max}-{\it Rg}_{min} = 2m; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("{\it Rg}_{max} = 10m; {\it P}_{m} = 0.05m; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it A}_{V} = 30°; \mu_p = 1.06; {\it S}_{px} = 0.56µm" )
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; \mu_p = 1.06; {\it S}_{px} = 0.56µm")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [8m,10m]; {\it P}_{m} = 0.05m; \mu_p = 1.06; {\it A}_{V} = 30°")
% title("[{\it Rg}_{min},{\it Rg}_{max}] = [0.5m,2m]; {\it P}_{m} = 0.1m; {\it A}_{V} = 20°; {\it S}_{px} = 0.56µm" )
title("LFC 6" )
set(gca, 'FontSize', 20);
pbaspect([1 1 1]);
xlabel('\mu_p');
% 
% xlabel('{\it S}_{px} (µm)')
% xlabel('{\it A}_{V} (°)')
% xlabel('{\it Rg}_{max} (m)')
% xlabel('{\it Rg}_{max}-{\it Rg}_{min} (m)')
% xlabel('{\it P}_{m} (m)')

%%
dF = F*0.008;

h = (f_M*1.7*F)/(1.7*F-f_M)-(f_M*(1.7*F+dF))/(1.7*F+dF-f_M);
(h*(1.7*F)^2-2*h*f_M*1.7*F+h*f_M^2)/(F*(f_M^2-h*1.7*F+h*f_M))






