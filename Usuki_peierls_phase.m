%% Project 2: Usuki Transfer Matrix Technique
clear all
close all
% Constants (everything is in SI units)
h = 6.626e-34;
hbar = h/(2*pi);
e = 1.602e-19;
q = -e;
a = 1e-9;
m0 = 9.11e-31;
mstar_gaas = 0.067*m0;
Disorder_num = 0;
numpoints = 100;
numwidths = 10;
widths = linspace(10e-9,700e-9,numwidths);
width = 100;
Nx = 2*numpoints;
Ny = 2*numpoints;
x = linspace(0,a*Nx,Nx);
y = linspace(0,a*Ny,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
numpoints = Nx;
t = (hbar^2)/(2*mstar_gaas*a^2);
numfermis = 1;
num_potentials = 1;
A_potential = zeros(numpoints,numpoints,3);
for k = 1:Nx
A_potential(:,k,1) = 0e-11*(k-numpoints/2);
end
% for i = 1:Nx
%     for j = 1:Ny
%         if i == 1
%             A_potential(i,j,:) = 0;
%         end
%     end
% end
% A_potential(:,:,1) = rot90(A_potential(:,:,1));
% A_potential(:,:,2) = rot90(A_potential(:,:,2));
phi_potential = zeros(numpoints,numpoints);
Tmatrix = zeros(30,30,numfermis);
% fermi_level = 1e-3*e*25;

B =0; %tesla
for i = 1:Ny
Hamiltonian_l_minus1(i,i) = -t*exp(-1j*(q/hbar)*A_potential(10,i,1)*dx);
Hamiltonian_l_plus1(i,i) = conj(Hamiltonian_l_minus1(i,i));
end
fermis = 1e-3*e*linspace(2,2,numfermis);

% potentials = linspace(1.5,1.5,num_potentials);
for k = 1:num_potentials


Potential = 0*ones(Nx,Ny);
% rectangle channel potential
for i = 1:Nx
    for j = 1:Ny
        if (floor(Nx/2)-floor(width/2)+25 < i && floor(Nx/2)+floor(width/2)-25 > i && (floor(Ny/2)-floor(width/4) > j+20 || floor(Ny/2)+floor(width/4) < j+20))
            Potential(i,j) = 1*1e4*t;
        end
    end
end

for i = 1:Nx
    for j = 1:Ny
        if (floor(Nx/2)-floor(width/2)+25 < i && floor(Nx/2)+floor(width/2)-25 > i && (floor(Ny/2)-floor(width/4) > j-20 || floor(Ny/2)+floor(width/4) < j-20))
            Potential(i,j) = 1*1e4*t;
        end
    end
end
% for i = 1:Nx/2
%     for j = 1:Ny
%         if j > 1.00*i && j < Ny-1.00*i
%             Potential(i,j) = fermi;
%             
%         end
%     end
% end
% 
% for i = Nx/2:Nx
%     for j = 1:Ny
%         if j > 1.00*(Nx-i) && j < (Ny - 1.00*(Nx-i))
%             Potential(i,j) = fermi-potentials(k)*e;
%         end
%     end
% end

Potential(:,1) = 1e4*t;
Potential(:,Ny) = 1e4*t;
Potential = rot90(Potential);
% for i = 1:Disorder_num
%     idx(i)=randperm(length(Potential),1);
%     idx2(i)=randperm(length(Potential),1);
%     if floor(Nx/2)-floor(width/2)+25 < idx2(i) && idx2(i) < floor(Nx/2)+floor(width/2)-25
%         Potential(idx(i),idx2(i)) = 1e4*t;
%     end
% end
for z = 1:numfermis
ef = fermis(z);

% Define potential energy

Hamiltonian_l = zeros(Ny,Ny);    
for l = 1:Nx
%     for i = 1:numpoints
%         if l == 1
%             Hamiltonian_l_plus1(i,i) = -t*exp(1j*(q/hbar)*A_potential(l+1,i,1)*dx);
%             Hamiltonian_l_minus1(i,i) = conj(Hamiltonian_l_plus1(i,i));
%         elseif l == Nx
%             Hamiltonian_l_minus1(i,i) = -t*exp(-1j*(q/hbar)*A_potential(l-1,i,1)*dx);
%             Hamiltonian_l_plus1(i,i) = conj(Hamiltonian_l_minus1(i,i));
%         else
%             Hamiltonian_l_minus1(i,i) = -t*exp(-1j*(q/hbar)*A_potential(l-1,i,1)*dx);
%             Hamiltonian_l_plus1(i,i) = -t*exp(1j*(q/hbar)*A_potential(l+1,i,1)*dx);
%         end
%     end
        
    
    for i = 1:Ny
        Hamiltonian_l(i,i) = 4*t+Potential(l,i)+phi_potential(l,i);
        if i == 1
            Hamiltonian_l(i,i+1) = -t*exp(1j*(q/hbar)*A_potential(l,i+1,2)*dy);
        end
        if i == Ny
            Hamiltonian_l(i,i-1) = -t*exp(-1j*(q/hbar)*A_potential(l,i-1,2)*dy);
        end
        if i > 1 && i < Ny
            Hamiltonian_l(i,i+1) = -t*exp(1j*(q/hbar)*A_potential(l,i+1,2)*dy);
            Hamiltonian_l(i,i-1) = -t*exp(-1j*(q/hbar)*A_potential(l,i-1,2)*dy);
        end
    end
%     Hamiltonian_l(1,2)-Hamiltonian_l(2,1)
     if l == 1
      check = (Hamiltonian_l);
      check_minus1 = Hamiltonian_l_minus1;
      check_plus1 = Hamiltonian_l_plus1;
     end
%     Hamiltonian_l = rot90(Hamiltonian_l);
    if l == 1
        T_0 = [zeros(Ny),eye(Ny);-inv(Hamiltonian_l_plus1)*Hamiltonian_l_minus1,...
            (Hamiltonian_l_plus1)\(Hamiltonian_l-ef*eye(Ny))];
        [somevectors,somevalues] = eig(T_0);
        % The sort_eig function takes the unsorted eigenvectors/values from
        % eig and sorts them into forward propagating, forward decaying,
        % back propagating and back evanescent.
        [sortedvalues,sortedvectors,forwardmodes] = sort_eig2(somevectors,somevalues);
        T_0new = (sortedvectors);
        test = abs(sortedvalues);
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
       Tl = [zeros(Ny),eye(Ny);-eye(Ny),...
                inv(Hamiltonian_l_plus1)*(Hamiltonian_l-ef*eye(Ny))];
       else
           Tl = inv(T_0new);
       end
%        test_val = sum(sum(abs(Tl)))
       P2l = inv([Tl(numpoints+1:2*numpoints,1:numpoints)*C2l+...
              (Tl(numpoints+1:2*numpoints,numpoints+1:2*numpoints))]);
       P1l = -P2l*Tl(numpoints+1:2*numpoints,1:numpoints)*C1l;
       p2density(l,:,:) = P2l;
       p1density(l,:,:) = P1l;
    end
end
T0_end = ([zeros(Ny),eye(Ny);-inv(Hamiltonian_l_plus1)*Hamiltonian_l_minus1,...
    (Hamiltonian_l_plus1)\(Hamiltonian_l-ef*eye(Ny))]);
[somevectorsout,somevaluesout] = eig((T0_end));
[sortedvaluesout,sortedvectorsout,ignorethis] = sort_eig2(somevectorsout,somevaluesout);
Tl_end = inv(sortedvectorsout);
% T0_end = inv(T_0);    
lfinal = Tl_end*[C1l,C2l;zeros(Ny),eye(Ny)]*[eye(Ny),zeros(Ny);P1l,P2l];
for i = 1:forwardmodes
    for j = 1:Nx
        propagating_speeds_in(i) = (2*t*dx/(hbar))*sum(abs(sortedvectors(j,:)*transpose(sortedvectors(j,:))).^2)*...
            sin((log(sortedvalues(i))/(1j)-(e/hbar)*A_potential(1,j,1))*dx);
    end
end
N_prop_forward = forwardmodes;
C1final = lfinal(1:forwardmodes,1:forwardmodes);
C1final2 = lfinal(Ny+1:Ny+forwardmodes,Ny+1:Ny+forwardmodes);
%Calculate T matrix out
% The sort_eig function takes the unsorted eigenvectors/values from
% eig and sorts them into forward propagating, forward decaying,
% back propagating and back evanescent.

for i = 1:forwardmodes
    for j = 1:Nx
        propagating_speeds_out(i) = (2*t*dx/(hbar))*sum(abs(sortedvectorsout(j,:)*transpose(sortedvectorsout(j,:))).^2)*...
            sin((log(sortedvaluesout(i))/(1j)-(e/hbar)*A_potential(Nx,j,1))*dx);
    end
end
%propagating_speeds_out = propagating_speeds_in;
T = zeros(N_prop_forward,N_prop_forward);
for i = 1:N_prop_forward
    for j = 1:N_prop_forward
        T(i,j) = (abs(propagating_speeds_out(i)/(propagating_speeds_in(j))))*abs(C1final(i,j))^2;
    end
end
Tmatrix(1:forwardmodes,1:forwardmodes,z) = T;
Gfinal(z) = (sum(sum(T)));
end
plot(fermis,real((Gfinal)));

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
for l = 1:numpoints
    nmode(l,j,i) = abs(phi(j,l,i))^2;
end
end
end
% Sum all modes together for electron density
density = zeros(numpoints,numpoints);
for i = 1:N_prop_forward
    density = density + nmode(:,:,i);
end
figure
imagesc(x,y,(density(:,1:numpoints)))
title("Electron Density")
xlabel("x-position (m)")
ylabel("y-position (m)")

% current(k) = calculate_ballistic_current(e*potentials(k),Tmatrix);
% percent_complete = 100*k/num_potentials;
% disp("Percent complete: " +percent_complete)
end
% figure
% plot(potentials,current)