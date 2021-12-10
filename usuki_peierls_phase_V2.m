function [usuki_currentx,usuki_currenty] = usuki_peierls_phase(Ax,Ay,Az,phi,numpoints)

% Constants (everything is in SI units)
h = 6.626e-34;
hbar = h/(2*pi);
e = 1.602e-19;
q = e;
a = 2e-9;               % size in meter
m0 = 9.11e-31;          % [kg]
mstar_GaAs = 0.067*m0;  % [kg]
% Disorder_num = 0;

% numpoints = numpoints; %size(Ax2Usuki,2);
Nx =  numpoints;
Ny =  numpoints;
usuki_currentx=zeros(numpoints,numpoints);
usuki_currenty=zeros(numpoints,numpoints);
x = linspace(0,a*Nx,Nx);
y = linspace(0,a*Ny,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
T = (hbar^2)/(2*mstar_GaAs*a^2); % kinetic energy operator  % Hopping Energy * 10,000 [J]
%% Fermi level configuration  Boundary condition quantum wires are modeled as square well waveguides
% Nano-Electronic Devices: Semiclassical and Quantum Transport Modeling
% edited by Dragica Vasileska, Stephen M. Goodnick
numfermis = 1;


% fermis = 1e-3*e*linspace(2,2,numfermis);
fermis = 1e-3*e*linspace(5,5,numfermis); %just use one to get current density
%% A vector potential
A_pot_valu = 10;
A_potential = zeros(numpoints,numpoints,3);

% B = 0.00023; %T


A_potential(:,:,1) = Ax(:,:,numpoints/2); % Assuming middle of z-direction
A_potential(:,:,2) = Ay(:,:,numpoints/2); % 
A_potential(:,:,3) = Az(:,:,numpoints/2); %

% figure
% A2d = A_potential(:,:,1);
% plot(A2d)
% for i = 1:Nx
%     for j = 1:Ny
%         if i == 1
%             A_potential(i,j,:) = 0;
%         end
%     end
% end
% A_potential(:,:,1) = rot90(A_potential(:,:,1));
% A_potential(:,:,2) = rot90(A_potential(:,:,2));
phi_potential = phi;
Tmatrix = zeros(30,30,numfermis);
% fermi_level = 1e-3*e*25;

% imagesc(real(Hamiltonian_l_minus1));colorbar
%% initialize potential based on kinetic energy

% potentials = linspace(1.5,1.5,num_potentials);
% for k = 1:num_potentials
% potential_valu = 0;
% Potential_op = 0*potential_valu*ones(Nx,Ny);
% rectangle channel potential
% % for i = 1:Nx
% %     for j = 1:Ny
% %         if (floor(Nx/2)-floor(width/2)+25 < i && floor(Nx/2)+floor(width/2)-25 > i ...
% %                 && (floor(Ny/2)-floor(width/4) > j+20 || floor(Ny/2)+floor(width/4) < j+20))
% %             Potential_op(i,j) = 1*1e4*T;
% %         end
% %     end
% % end
% %
% % for i = 1:Nx
% %     for j = 1:Ny
% %         if (floor(Nx/2)-floor(width/2)+25 < i && floor(Nx/2)+floor(width/2)-25 > i ...
% %                 && (floor(Ny/2)-floor(width/4) > j-20 || floor(Ny/2)+floor(width/4) < j-20))
% %             Potential_op(i,j) = 1*1e4*T;
% %         end
% %     end
% % end
Lnth = floor(Nx/4);
width = floor(Nx/10);
[Domain_air,Domain_FET] = Domain_Rectangle_Mosfet(Nx,Ny,2, Nx/2);

Potential_op = 0*T .* Domain_air + 1e4*T.* Domain_FET ;
Potential_op = rot90(Potential_op);
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
% hardwall potential on top and bottom of Nanowire
% Potential_op(:,1) = 1e4*T;
% Potential_op(:,Ny) = 1e4*T;
% Potential_op = rot90(Potential_op);
% %  Create potential operator for a random disordered system
% for i = 1:Disorder_num
%     idx =randperm(length(Potential_op),2);  % choose random locations to define different potentials at those points
%     %     idx2(i)=randperm(length(Potential_op),1); % choose random locations to define different potentials at those points
%     Potential_op(idx(1),idx(2)) = 1e4*T;
% % end
% figure
% imagesc((Potential_op));colorbar
% title('rectangle channel Potential op')
% print(gcf,'rectangle channel_Potential_op.jpg','-djpeg');
%% Hamiltonian
Gfinal = zeros(numfermis,1);
for z = 1:numfermis
    ef = fermis(z);
    
    % Define potential energy
    
    
    % Define potential energy in Hamiltonian
    fprintf('Getting Eigs and Sorting Modes...');
    disp(['Fermi level number ', num2str(z)])
    for l = 1:Nx
        Hamiltonian_l = zeros(Ny,Ny);
        %% Hamiltonian
        Hamiltonian_l_minus1 = T*eye(Ny);
        for i = 1:Ny
                Hamiltonian_l_minus1(i,i) = -T*exp(1j*(q/hbar)*A_potential(i,l,1)*dx); %sam
        %     Hamiltonian_l_minus1(i,i) = -T*exp(1j*2*pi*B*i);

        end
        Hamiltonian_l_plus1  = conj(Hamiltonian_l_minus1 );
        for i = 1:Ny
            Hamiltonian_l(i,i) = 4*T+Potential_op(i,l)+phi_potential(i,l);
            if i == 1
                Hamiltonian_l(i,i+1) = -T*exp(1j*(q/hbar)*A_potential(i+1,l,2)*dy);
            end
            if i == Ny
                Hamiltonian_l(i,i-1) = -T*exp(-1j*(q/hbar)*A_potential(i-1,l,2)*dy);
            end
            if i > 1 && i < Ny
                Hamiltonian_l(i,i+1) = -T*exp(1j*(q/hbar)*A_potential(i+1,l,2)*dy);
                Hamiltonian_l(i,i-1) = -T*exp(-1j*(q/hbar)*A_potential(i-1,l,2)*dy);
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
        if forwardmodes == 0
            disp("No propagating modes!!")
            usuki_currentx(:,:) = 0;
            usuki_currenty(:,:) = 0;
            return
        end
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
    T0_end = ([zeros(Ny),eye(Ny);-inv(Hamiltonian_l_plus1)*Hamiltonian_l_minus1,...
        (Hamiltonian_l_plus1)\(Hamiltonian_l-ef*eye(Ny))]);
%     T0_end = inv(T_0);
    [somevectorsout,somevaluesout] = eig((T0_end));
    [sortedvaluesout,sortedvectorsout,ignorethis] = sort_eig2(somevectorsout,somevaluesout);
    Tl_end = inv(sortedvectorsout);
    lfinal = Tl_end*[C1l,C2l;zeros(Ny),eye(Ny)]*[eye(Ny),zeros(Ny);P1l,P2l];
    for i = 1:forwardmodes
        for j = 1:Nx
            propagating_speeds_in(i) = (2*T*dx/(hbar))*sum(abs(sortedvectors(j,:)*transpose(sortedvectors(j,:))).^2)*...
                sin((log(sortedvaluesout(i))/(1j)-(e/hbar)*A_potential(j,1,1))*dx);
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
            propagating_speeds_out(i) = (2*T*dx/(hbar))*sum(abs(sortedvectorsout(j,:)*transpose(sortedvectorsout(j,:))).^2)*...
                sin((log(sortedvalues(i))/(1j)-(e/hbar)*A_potential(j,Nx,1))*dx);
        end
    end
    %propagating_speeds_out = propagating_speeds_in;
    T_matrix = zeros(N_prop_forward,N_prop_forward);
    for i = 1:N_prop_forward
        for j = 1:N_prop_forward
            T_matrix(i,j) = (abs(propagating_speeds_out(i)/(propagating_speeds_in(j))))*abs(C1final(i,j))^2;
        end
    end
    Tmatrix(1:forwardmodes,1:forwardmodes,z) = T_matrix;
    Gfinal(z) = (sum(sum(T_matrix)));
end
% figure
% plot(fermis,real((Gfinal)));
% ylabel('Conductance')
% xlabel('Gate Voltage')
% print(gcf,'Rectangles_Conductance.jpg','-djpeg');
%% Finding Density

x = linspace(0,a*Nx,Nx);
y = linspace(0,a*Ny,Ny);
% Initialize phi and set the last slice of phi equal to the last slice p
phi_current = zeros(Nx,numpoints,numpoints);
phi_current(Nx,:,:) = p1density(Nx,:,:);
% back-propagate phi
for i = Nx-1:-1:1
    p1densitytemp = reshape(p1density(i,:,:),numpoints,numpoints);
    p2densitytemp = reshape(p2density(i,:,:),numpoints,numpoints);
    phitemp  = reshape(phi_current(i+1,:,:),numpoints,numpoints);
    phi_current(i,:,:) = p1densitytemp+p2densitytemp*phitemp;
end
% using phi, calculate electron density for each mode
for i = 1:N_prop_forward
    for j = 1:numpoints
        for l = 1:numpoints
            nmode(l,j,i) = abs(phi_current(j,l,i))^2;
        end
    end
end
% Sum all modes together for electron density
density = zeros(numpoints,numpoints);
for i = 1:N_prop_forward
    density = density + nmode(:,:,i);
end
current_densityx = zeros(numpoints,numpoints);
current_densityy = zeros(numpoints,numpoints);
for i = 1:N_prop_forward
    [gradphix,gradphiy] = gradient(phi_current(:,:,i));
    [gradconjx,gradconjy] = gradient(conj(phi_current(:,:,i)));
    current_densityx = current_densityx+(hbar/(2*m0*1j))*(conj(phi_current(:,:,i)).*...
        gradphix-phi_current(:,:,i).*gradconjx);
    current_densityy = current_densityy+(hbar/(2*m0*1j))*(conj(phi_current(:,:,i)).*...
        gradphiy-phi_current(:,:,i).*gradconjy);
end
electron_density = figure;
imagesc(x,y,(density(:,2:numpoints)))
title("Electron Density")
xlabel("x-position (m)")
ylabel("y-position (m)")
pause(0.5)
close(electron_density)
% print(gcf,'rectangle_density.jpg','-djpeg');
%writematrix(density,'rectangleDensity.txt')
% current = -q*calculate_ballistic_current(fermis(numfermis),Tmatrix);
% percent_complete = 100*k/num_potentials;
% disp("Percent complete: " +percent_complete)
% end
% figure()
% imagesc(x,y,rot90(real(current_densityx(:,2:numpoints))));
usuki_currentx = 2*rot90(current_densityx); %indices were switched during calculation of phi so switching back with rot90
usuki_currenty = 2*rot90(current_densityy);