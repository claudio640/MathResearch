%% Defining values and stream function
%close all; clear;
% Testing values where N is the number of points, L is the length of the
% interval, r is the radius of modon and c is its speed
N = 64;
L = 3;
R = 1;
c = 1;

% h is the number of intervals. Build x and y points and get their sizes
    h = (2*L)/(N-1);
    x = -L:h:L;  y = x;
   [yy, xx] = meshgrid(y,x);
    
    xSize = size(x,2);
    ySize = size(y,2);
   
% Construction of stream function
    Bes = @(x)besselj(1,x);
    z1 = fzero(Bes,3);
    Sqrt_a = z1/R;
    Jprime = (besselj(0,z1) - besselj(2,z1))/2;
% Get b and K as defined in the lectures
    b = 2*c/(Sqrt_a*Jprime);
    K = -c*R^2;
% Build the PSI in out components
    P_tilde_in  = (b*yy.*besselj(1,Sqrt_a.*sqrt(xx.^2 + yy.^2)))./(sqrt(xx.^2 +yy.^2)) -c*yy;
    P_tilde_out = (K.*yy)./(xx.^2 + yy.^2);
% only consider the values where they are defined
    P_tilde_inR = P_tilde_in.*(xx.^2 +yy.^2 <= R^2);
    P_tilde_inR(isnan(P_tilde_inR)) = 0;
    
    
    P_tilde_outR = P_tilde_out.*(xx.^2 +yy.^2 > R^2);
    P_tilde_outR(isnan(P_tilde_outR)) = 0;
    
% join PSI together and reshape into a vector to compute u, v and w   
    P_tilde = P_tilde_inR + P_tilde_outR;
    P_tilde_fullR = P_tilde(2:end-1, 2:end-1);
    
    CA = load('CA.mat');
    CA = struct2cell(CA);
    CA = cat(2,CA{:});
    
    PSI = c*yy(2:end-1,2:end-1) + P_tilde_fullR +8*CA/norm(CA);
    
figure(99);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),PSI, [-3:0.1:3])
axis equal
title('cy + $\tilde{\Psi}$ + 8*CA/norm(CA)', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

figure(999);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),PSI-8*CA/norm(CA), [-3:0.1:3])
axis equal
title('cy + $\tilde{\Psi}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar
    
    P_tilde_full = reshape(P_tilde_fullR, [(xSize-2)*(ySize-2) 1]);
    
    % Get the d/dx , d/dy and laplacian matrices    
    Mx = dxMatrix(N,h,xSize); My = dyMatrix(N,h,xSize); ML = LMatrix(N,h);
% Compute u, v and w and shape them into a matrix to use as coefficients    
    u_tilde = My*P_tilde_full; v_tilde = -1*(Mx*P_tilde_full); w_tilde = -1*(ML*P_tilde_full);
    
%% Velocity u

 PSI_in_dy = (Sqrt_a*b*yy.^2 .*(besselj(0,Sqrt_a.*sqrt(xx.^2 + yy.^2)) - besselj(2,Sqrt_a.*sqrt(xx.^2 + yy.^2))))./(2*(xx.^2 + yy.^2)) - (b*yy.^2 .*besselj(1,Sqrt_a.*sqrt(xx.^2 + yy.^2)))./(xx.^2 +yy.^2).^(3/2)  + b*besselj(1,Sqrt_a.*sqrt(xx.^2 + yy.^2))./sqrt(xx.^2 +yy.^2)  -c; 
 PSI_out_dy = (K*(xx.^2 -yy.^2))./(xx.^2 + yy.^2).^2;
 
 PSI_in_dy = PSI_in_dy.*(xx.^2 + yy.^2 <= R^2);
 PSI_out_dy = PSI_out_dy.*(xx.^2 + yy.^2 > R^2);
 
 PSI_full_dy = PSI_in_dy + PSI_out_dy;

u_tildeR = reshape(u_tilde, [xSize-2 ySize-2]);

u_tildeR(:,1) = u_tildeR(:,1) -P_tilde(2:end-1,1)/(2*h);
u_tildeR(:,end) = u_tildeR(:,end) +P_tilde(2:end-1,end)/(2*h);

figure(203);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),PSI_full_dy(2:end-1,2:end-1),21)
axis equal
title(' Calculus $\tilde{u}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

figure(202);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),u_tildeR,21)
axis equal
title('Discrete $\tilde{u}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar   

figure(205);  clf
surf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),abs(PSI_full_dy(2:end-1,2:end-1)-u_tildeR),'edgecolor', 'none')
title('Discrete $\tilde{u}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
%% Velocity v

PSI_in_dx = (Sqrt_a*b.*xx.*yy.*(besselj(0,Sqrt_a.*sqrt(xx.^2 + yy.^2)) - besselj(2,Sqrt_a.*sqrt(xx.^2 + yy.^2))))./(2*(xx.^2 + yy.^2)) - (b*xx.*yy.*besselj(1,Sqrt_a.*sqrt(xx.^2 + yy.^2)))./(xx.^2 +yy.^2).^(3/2);
PSI_in_dx = -1*PSI_in_dx;
PSI_out_dx = (2*K*xx.*yy)./(xx.^2 + yy.^2).^2;
 
 PSI_in_dx = PSI_in_dx.*(xx.^2 + yy.^2 <= R^2);
 PSI_out_dx = PSI_out_dx.*(xx.^2 + yy.^2 > R^2);
 
 PSI_full_dx = PSI_in_dx + PSI_out_dx;

v_tildeR = reshape(v_tilde, [xSize-2 ySize-2]);

v_tildeR(1,:) = v_tildeR(1,:) +P_tilde(1,2:end-1)/(2*h);
v_tildeR(end,:) = v_tildeR(end,:) -P_tilde(end,2:end-1)/(2*h);

figure(200);  clf
contourf(xx(3:end-2,3:end-2),yy(3:end-2,3:end-2),PSI_full_dx(3:end-2,3:end-2),21)
axis equal
title(' Calculus $\tilde{v}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

figure(201);  clf
contourf(xx(3:end-2,3:end-2),yy(3:end-2,3:end-2),v_tildeR(2:end-1,2:end-1),21)
axis equal
title('Discrete $\tilde{v}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

figure(207);  clf
surf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),abs(PSI_full_dx(2:end-1,2:end-1)-v_tildeR),'edgecolor', 'none')
title('Discrete $\tilde{u}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
%% Vorticity w

w_tilde_full = Sqrt_a^2 *(c.*yy + P_tilde_inR).*(xx.^2 +yy.^2 <= R^2);

figure(200);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),w_tilde_full(2:end-1,2:end-1),21)
axis equal
title(' Calculus $\tilde{w}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

w_tildeR = reshape(w_tilde, [xSize-2 ySize-2]);
w_tildeR(:,1) = w_tildeR(:,1) -P_tilde(2:end-1,1)/(h^2);
w_tildeR(:,end) = w_tildeR(:,end) -P_tilde(2:end-1,end)/(h^2);
w_tildeR(1,:) = w_tildeR(1,:) -P_tilde(1,2:end-1)/(h^2);
w_tildeR(end,:) = w_tildeR(end,:) -P_tilde(end,2:end-1)/(h^2);

figure(201);  clf
contourf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),w_tildeR,21)
axis equal
title('Discrete $\tilde{w}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')
colorbar

figure(205);  clf
surf(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1),abs(w_tilde_full(2:end-1,2:end-1)-w_tildeR),'edgecolor', 'none')
title('Discrete $\tilde{u}$', 'interpreter', 'latex')
xlabel('x'); ylabel('y')