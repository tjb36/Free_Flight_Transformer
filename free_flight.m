function psi_out = free_flight( psi_in , DX1o, DX2o, DX3o, mult1, mult2, mult3, t_flight )
 
M1 = size(psi_in,2); % Number of points in x direction
M2 = size(psi_in,1); % Number of points in y direction
M3 = size(psi_in,3); % Number of points in z direction
M = M1 * M2 * M3;    % Total number of grid points

num_mu = mult1 * mult2 * mult3;   % Total number of terms in sum

% Original lattice lengths
L1o = DX1o * M1;
L2o = DX2o * M2;
L3o = DX3o * M3;

% Expanded lattice real space grid spacings
DX1 = DX1o * mult1;
DX2 = DX2o * mult2;
DX3 = DX3o * mult3;

% Expanded lattice box lengths
L1 = DX1 * M1;
L2 = DX2 * M2;
L3 = DX3 * M3;

% Expanded lattice k space grid spacings
DK1 = 2*pi / L1;
DK2 = 2*pi / L2;
DK3 = 2*pi / L3;

% Original lattice k space grid spacings
DK1o = DK1 * mult1;
DK2o = DK2 * mult2;
DK3o = DK3 * mult3;

% Original lattice offsets
a1o = -L1o / 2;
a2o = -L2o / 2;
a3o = -L3o / 2;

% Expanded lattice offsets
a1 = - L1 / 2;
a2 = - L2 / 2;
a3 = - L3 / 2;

% Original lattice real space grid vectors
x1o = ( 0:DX1o:(M1-1)*DX1o ) + a1o;
x2o = ( 0:DX2o:(M2-1)*DX2o ) + a2o;
x3o = ( 0:DX3o:(M3-1)*DX3o ) + a3o;
x2o = x2o.';                        % Change dimension orientation to allow implicit expansion in Matlab
x3o = permute( x3o ,[1,3,2]);       % Change dimension orientation to allow implicit expansion in Matlab

% Expanded lattice real space grid vectors
x1vals = x1o * mult1;
x2vals = x2o * mult2;
x3vals = x3o * mult3;

% Expanded lattice k space grid vectors
k1vals = ifftshift( (-pi/DX1):DK1:(pi/DX1 - DK1) );
k2vals = ifftshift( (-pi/DX2):DK2:(pi/DX2 - DK2) );
k3vals = ifftshift( (-pi/DX3):DK3:(pi/DX3 - DK3) );

% Original lattice k space grid vectors
k1o = k1vals * mult1;
k2o = k2vals * mult2;
k3o = k3vals * mult3;
k2o = k2o.';                  % Change dimension orientation to allow implicit expansion in Matlab
k3o = permute( k3o ,[1,3,2]); % Change dimension orientation to allow implicit expansion in Matlab

%%% PREFACTORS WHICH DO NOT DEPEND ON LOOP VARIABLES %%%          
Ekfac1 = exp( -i/2 * t_flight * k1o.^2 );  % Prefactors used in exponential in Eq. (42c) 
Ekfac2 = exp( -i/2 * t_flight * k2o.^2 );
Ekfac3 = exp( -i/2 * t_flight * k3o.^2 );
kdafac1 = exp( i * k1o * (a1-a1o) );       % Prefactors used in exponential in Eq. (42c) 
kdafac2 = exp( i * k2o * (a2-a2o) );
kdafac3 = exp( i * k3o * (a3-a3o) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi_out = zeros(M,1); % Preallocate

for mu1 = 0:(mult1-1)
    for mu2 = 0:(mult2-1)
        for mu3 = 0:(mult3-1)
            
            
            %%% PREFACTORS WHICH DEPEND ON LOOP VARIABLES %%%          
            xofac1 = exp( -i*x1o*mu1*DK1o/mult1 ); % Prefactors used in exponential in Eq. (42d) 
            xofac2 = exp( -i*x2o*mu2*DK2o/mult2 );
            xofac3 = exp( -i*x3o*mu3*DK3o/mult3 );
            
            kofac1 = exp( -i*k1o*t_flight*mu1*DK1o/mult1 ); % Prefactors used in exponential in Eq. (42c) 
            kofac2 = exp( -i*k2o*t_flight*mu2*DK2o/mult2 );
            kofac3 = exp( -i*k3o*t_flight*mu3*DK3o/mult3 );
            
            xfac1  = exp( i*x1vals*mu1*DK1o/mult1 ); % Prefactors used in exponential in Eq. (42b) 
            xfac2  = exp( i*x2vals*mu2*DK2o/mult2 );
            xfac3  = exp( i*x3vals*mu3*DK3o/mult3 );
            
            % Prefactor used in exponential in Eq. (42b) 
            dtfac  = exp( -i/2 * t_flight *( (mu1*DK1o/mult1)^2 + (mu2*DK2o/mult2)^2 + (mu3*DK3o/mult3)^2 ) );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            fieldmu = psi_in .* xofac1 .* xofac2 .* xofac3;        % produces kernel for A in Eq. (42d)
            fieldmu = fftn(fieldmu);	                           % produces A in Eq. (42d)
            fieldmu = fieldmu .* (Ekfac1.*Ekfac2.*Ekfac3) .* (kdafac1.*kdafac2.*kdafac3) .* (kofac1.*kofac2.*kofac3); % kernel for B in Eq. (42c)
            fieldmu = ifftn(fieldmu);	                           % produces B in Eq. (42c)
            fieldmu = reshape( permute(fieldmu,[3,1,2]) , [M 1] ); % Convert 3D form back to 1D form, ready for last loop

            % Final loop calculates f in Eq. 42b and accumulates to sum in Eq. 42a
            % (THIS PART NOT VECTORIZED YET %
            for cur_x_i=1:M
                j1 = floor(floor((cur_x_i-1)/M3)/M2);   % index on new lattice
                j2 = mod( floor((cur_x_i-1)/M3) , M2 ); % -1 and +1 used to allow Matlab's [1] start indexing, but still retain Deuar's formulas
                j3 = mod( (cur_x_i-1) , M3 );
                jbar1 = mult1*j1;		                % index jbar on the "big" virtual lattice
                jbar2 = mult2*j2;	            
                jbar3 = mult3*j3;               
                n1prime = mod(jbar1,M1);                % index n' on the transfromed mu lattice part
                n2prime = mod(jbar2,M2);        
                n3prime = mod(jbar3,M3);        
                primed_x_i = mod( (n3prime+M3*(n2prime+M2*n1prime)) , M );       % index in Eq. (42e)
                f = xfac1(j1+1) * xfac2(j2+1) * xfac3(j3+1) * dtfac;             % produces f in Eq. (42b)
                psi_out(cur_x_i) = psi_out(cur_x_i) + f * fieldmu(primed_x_i+1); % produces psi in Eq. (42a)
            end
            
        end
    end
end

end
