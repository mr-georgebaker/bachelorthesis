# Bachelor thesis

This repository contains the code for my bachelor thesis in physics at the University of Vienna about molecular dynamics based on the langevin equation. The detailed progress can be found here: [progress.pdf](progress.pdf)

<h3>Milestones</h3>

<ul>
<li>Implement and test Langevin Equation using velocity verlet integration</li>
<li>Langevin Dynamics for Lennard-Jones-Fluid</li>
<li>Implement Replica Exchange Algorithm </li>
<li>Multiprocessing implementation using OpenMP</li>
</ul>

<h4>Implement and test Langevin Equation using velocity verlet integration</h4>

<h5>Implementation</h5>
For a implementation where the force doesn't depend on the position of the other particles see [vv_xv.sce](velocity_verlet/vv_xv.sce) (positions and velocities returned) and [vv_x.sce](velocity_verlet/vv_x.sce) (faster version, only positions returned).

<h5>Test</h5>

If ![alt tag](http://latex.codecogs.com/png.latex? V(x) = \\frac{1}{2} k x^2 \\text{ and } \\gamma = 0) the equation reduces to a simple harmonic oscillator where one knows the probability density of the position and that the energy is conserved. Results are available in [vv_xv_nonboltzmann.sce](velocity_verlet/vv_xv_nonboltzmann.sce)

For the case where ![alt tag](http://latex.codecogs.com/png.latex? V(x) = \\frac{1}{2} k x^2 \\text{ and } \\gamma  \\neq 0) the probability distribution is a Boltzmann distribution. The implementation is available in [vv_x_boltzmann.sce](velocity_verlet/vv_x_boltzmann.sce)
