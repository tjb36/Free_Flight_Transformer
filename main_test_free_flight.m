% Test the function free_flight.m

clear;clc;close all;

Nx = 8;           
Ny = 6; 
Nz = 4; 
M = Nx*Ny*Nz;

mult1 = 2;
mult2 = 5;
mult3 = 3;

t_flight = 12.5;

% Set up real space grid
dx = 0.25;
dy = 0.3;
dz = 2/3;
x  = -(Nx*dx/2):dx:(Nx*dx/2-dx);
y  = -(Ny*dy/2):dy:(Ny*dy/2-dy);
z  = -(Nz*dz/2):dz:(Nz*dz/2-dz);
[X, Y, Z]  = meshgrid(x, y, z);

% Generate a random input wavefunction psi
rng('default');
rng(1990);
psi_3D_real = 2*(rand(Ny,Nx,Nz)-0.3);
rng(2008);
psi_3D_imag = i*(rand(Ny,Nx,Nz)-0.5);
psi_3D = psi_3D_real + psi_3D_imag;

tic
psi_1D_out = free_flight(psi_3D, dx, dy, dz, mult1, mult2, mult3, t_flight)
toc
