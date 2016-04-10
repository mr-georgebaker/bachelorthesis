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
See [velocity_verlet/velocity_verlet.sce](velocity_verlet/velocity_verlet.sce) and [velocity_verlet/velocity_verlet_fast.sce](velocity_verlet/velocity_verlet_fast.sce)

<h5>Test</h5>

If ![alt tag](http://latex.codecogs.com/png.latex? V(x) = \\frac{1}{2} k x^2 \\text{ and } \\gamma = 0) the equation reduces to a simple harmonic oscillator where one knows the probability density of the position and that the energy is conserved. Results are available in [velocity_verlet/velocity_verlet_energy_conservation.sce](velocity_verlet/velocity_verlet_energy_conservation.sce)
