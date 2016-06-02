clear

// Calculates the positions (x) for given particles 
// using the langevin equation and a second order integrator

// CONSTANTS

k_B = 1

function f = force(x)
// Returns the value for the force f = -k*x due to the potential V = (1/2)*k*x^2

// Input:   x: vector (coordinates)
// Output:  f: vector
// CONSTANTS: k
    k = 2
    f = -k*x
endfunction

function x = velocity_verlet(dt,n,T,g,m,x0,v0,forcefield)
// Returns a vector containing all positions obtained using 
// a generalization of the velocity verlet algorithm in order to
// integrate the langevin equation in one dimension

// Input:   dt: number (timestep)
//          n: number (amount of steps)
//          T: number (temperature)
//          g: number (gamma factor: friction coefficient)
//          m: number (mass of particle)
//          x0: rowvector (initial position)
//          v0: rowvector (initial velocity)
//          forcefield: function (a function which returns the values for the force)
// Output:  x: nxd matrix (coordinates where d is the dimension)
//          v: nxd matrix (velocities where d is the dimension)
    sigma = sqrt(2*k_B*T*g*m^(-1))
    d = length(x0)
    x(1,1:d) = x0
    for i = 1:1:n-1
        xi_n = grand(1,d,"nor",0,1)
        eta_n= grand(1,d,"nor",0,1)
        v_half = v0 + (1/2)*dt*(1/m)*forcefield(x(i,1:d)) - (1/2)*dt*g*v0 + (1/2)*sqrt(dt)*sigma*xi_n - (1/8)*dt^2*g*((1/m)*forcefield(x(i,1:d)) - g*v0) - (1/4)*dt^(3/2)*g*sigma*((1/2)*xi_n + (1/sqrt(3))*eta_n)
        x(i+1,1:d) = x(i,1:d) + dt*v_half + dt^(3/2)*sigma*(1/sqrt(12))*eta_n
        v0 = v_half + (1/2)*dt*(1/m)*forcefield(x(i+1,1:d)) - (1/2)*dt*g*v_half + (1/2)*sqrt(dt)*sigma*xi_n - (1/8)*dt^2*g*((1/m)*forcefield(x(i+1,1:d)) - g*v_half) - (1/4)*dt^(3/2)*g*sigma*((1/2)*xi_n + (1/sqrt(3))*eta_n)
    end
endfunction

// PARAMETERS

dt = 5*10^-3
n = 10000
T = 0.5
g = 0.3
m = 1
x0 = 0
v0 = 0
k = 2

x = velocity_verlet(dt,n,T,g,m,x0,v0,force)
