%% Project 2: Usuki Transfer Matrix Technique
clear all
% Constants (everything is in SI units)
h = 6.626e-34;
hbar = h/(2*pi);
e = 1.602e-19;
a = 2e-9;
m0 = 9.11e-31;
mstar_gaas = 0.067*m0;
numpoints = 100;
Nx = 2*numpoints;
Ny = 2*numpoints;
Channel_Width = floor(Nx/10);
x = linspace(0,a*Nx,Nx);
y = linspace(0,a*Ny,Ny);
numpoints = 2*numpoints;
t = (hbar^2)/(2*mstar_gaas*a^2);
numfermis = 1;
numwidths = 10;
widths = linspace(10e-9,700e-9,numwidths);
fermis = 1e-3*e*linspace(5,5,numfermis);

%potential = 0;
for z = 1:numfermis
ef = fermis(z);
width = 1e9*widths(2);
% Define potential energy to create loop
Potential = 1e4*t*ones(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if (i<(Nx-floor(Nx/4)+Channel_Width) && i>floor(Nx/4)-Channel_Width && j<(Ny-floor(Ny/4)+Channel_Width) && j>floor(Ny/4)-Channel_Width)
            Potential(i,j) = 0;
        end
    end
end

for i = 1:Nx
    for j = 1:Ny
        if (i<(Nx-floor(Nx/4)) && i>floor(Nx/4) && j<(Ny-floor(Ny/4)) && j>floor(Ny/4))
            Potential(i,j) = 1e4*t;
        end
    end
end

for i = 1:Nx
    for j = 1:Ny
        if (i>=(Nx-floor(Nx/4)) || i<=floor(Nx/4))
            if j>(floor(Ny/2)-Channel_Width/2) && j < floor(Ny/2)+Channel_Width/2
            Potential(i,j) = 0;
            end
        end
    end
end
Potential(:,1) = 1e4*t;
Potential(:,Ny) = 1e4*t;
Potential = rot90(Potential);
Nb = 60;
Bvec = linspace(0,0.0001,Nb);
for n = 1:Nb
B = Bvec(n); %tesla
Hamiltonian_l_minus1 = t*eye(Ny);
% Hamiltonian_l_plus1 = t*eye(Ny);
for i = 1:numpoints
    Hamiltonian_l_minus1 (i,i) = -t*exp(2*1j*pi*B*i);
end
Hamiltonian_l_plus1 = conj(Hamiltonian_l_minus1);

for l = 1:Nx
    Hamiltonian_l = zeros(Ny,Ny);
    for i = 1:Ny
        Hamiltonian_l(i,i) = 4*t+Potential(i,l);
        if i == 1
            Hamiltonian_l(i,i+1) = -t;
        end
        if i == Ny
            Hamiltonian_l(i,i-1) = -t;
        end
        if i > 1 && i < Ny
            Hamiltonian_l(i,i+1) = -t;
            Hamiltonian_l(i,i-1) = -t;
        end
    end
%     if l == 100
%         check = (Hamiltonian_l);
%     end
%     Hamiltonian_l = rot90(Hamiltonian_l);
    if l == 1
        T_0 = [zeros(Ny),eye(Ny);-inv(Hamiltonian_l_plus1)*Hamiltonian_l_minus1,...
            inv(Hamiltonian_l_plus1)*(Hamiltonian_l-ef*eye(Ny))];
        [somevectors,somevalues] = eig(T_0);
        % The sort_eig function takes the unsorted eigenvectors/values from
        % eig and sorts them into forward propagating, forward decaying,
        % back propagating and back evanescent.
        [sortedvalues,sortedvectors,forwardmodes] = sort_eig2(somevectors,somevalues);
        T_0new = (sortedvectors);
        C2l = zeros(Ny);
        C1l = eye(Ny);
        P2l = inv([T_0new(numpoints+1:2*numpoints,1:numpoints)*C2l+...
            (T_0new(numpoints+1:2*numpoints,numpoints+1:2*numpoints))]);
        P1l = -P2l*T_0new(numpoints+1:2*numpoints,1:numpoints)*C1l;
        Tl = T_0new;
    else
    
   % propagate through slices
       lnew = Tl*[C1l,C2l;zeros(Ny),eye(Ny)]*[eye(Ny),zeros(Ny);P1l,P2l];
       C1l = lnew(1:numpoints,1:numpoints);
       C2l = lnew(1:numpoints,numpoints+1:2*numpoints);
       if l ~= numpoints
        Tl = [zeros(Ny),eye(Ny);-inv(Hamiltonian_l_plus1)*Hamiltonian_l_minus1,...
                (Hamiltonian_l_plus1)\(Hamiltonian_l-ef*eye(Ny))];
       else
           Tl = inv(T_0new);
       end
       P2l = inv([Tl(numpoints+1:2*numpoints,1:numpoints)*C2l+...
              (Tl(numpoints+1:2*numpoints,numpoints+1:2*numpoints))]);
       P1l = -P2l*Tl(numpoints+1:2*numpoints,1:numpoints)*C1l;
       p2density(l,:,:) = P2l;
       p1density(l,:,:) = P1l;
    end
end
% T0_end = [zeros(Ny),inv(T_0new(Ny+1:2*Ny,1:Ny));eye(Ny),-T_0new(1:Ny,1:Ny)*...
%    inv(T_0new(Ny+1:2*Ny,1:Ny))];
T0_end = inv(T_0);
lnew =Tl*[C1l,C2l;zeros(Ny),eye(Ny)]*[eye(Ny),zeros(Ny);P1l,P2l];
propagating_speeds_in = hbar*log(sortedvalues(1:forwardmodes))/(1j*a*mstar_gaas);
N_prop_forward = forwardmodes;
C1final = lnew(1:forwardmodes,1:forwardmodes);
C1final2 = lnew(Ny+1:Ny+forwardmodes,Ny+1:Ny+forwardmodes);
%Calculate T matrix out
[somevectorsout,somevaluesout] = eig(inv(T0_end));
% The sort_eig function takes the unsorted eigenvectors/values from
% eig and sorts them into forward propagating, forward decaying,
% back propagating and back evanescent.
[sortedvaluesout,sortedvectorsout,ignorethis] = sort_eig2(somevectorsout,somevaluesout);
propagating_speeds_out = hbar*log(sortedvaluesout(1:forwardmodes))/(1j*a*mstar_gaas);
% propagating_speeds_out = propagating_speeds_in;
T = zeros(N_prop_forward,N_prop_forward);
for i = 1:N_prop_forward
    for j = 1:N_prop_forward
        T(i,j) = (abs(propagating_speeds_in(i)/(propagating_speeds_out(j))))*abs(C1final(i,j))^2;
    end
end
Gfinal(n) = (sum(sum(T)));
disp("progress is: "+string(n/Nb)); 
end

end

plot(Bvec ,Gfinal)
%% Finding Density

x = linspace(0,a*Nx,Nx);
y = linspace(0,a*Ny,Ny);
% Initialize phi and set the last slice of phi equal to the last slice p
phi = zeros(Nx,numpoints,numpoints);
phi(Nx,:,:) = p1density(Nx,:,:);
% back-propagate phi
for i = Nx-1:-1:1
    p1densitytemp = reshape(p1density(i,:,:),numpoints,numpoints);
    p2densitytemp = reshape(p2density(i,:,:),numpoints,numpoints);
    phitemp  = reshape(phi(i+1,:,:),numpoints,numpoints);
    phi(i,:,:) = p1densitytemp+p2densitytemp*phitemp;
end
% using phi, calculate electron density for each mode
for i = 1:N_prop_forward
for j = 1:numpoints
for k = 1:numpoints
    nmode(k,j,i) = abs(phi(j,k,i))^2;
end
end
end
% Sum all modes together for electron density
density = zeros(numpoints,numpoints);
for i = 1:N_prop_forward
    density = density + nmode(:,:,i);
end
figure
imagesc(x,y,density(:,1:numpoints))
title("Electron Density")
xlabel("x-position (m)")
ylabel("y-position (m)")