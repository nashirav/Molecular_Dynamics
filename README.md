Molecular_Dynamics
==================
Implementation of the simulated annealing and replica exchange Monte Carlo algorithms for protein folding in the HP model
in Python and using NumPy library.

Hydrophobic-polar protein folding (HP) model is used in the study of the general principles of protein folding.
The idea of the HP model is based on the observation that a key role in the process of folding
has the hydrophobic effect - tendency of hydrophobic amino acids to aggregate and 'hide' from the water molecules.
Amino acids are over the alphabet {H,P}, where H is hydrophobic and P polar amino acid
and there are on the square lattice.

Markov chain Monte Carlo (MCMC) sampling using a Metropolis-Hastings allows to 
sampling the set of possible configurations, according to the Boltzmann distribution propability (in case of Metropolis-Hastings).

For more details see [1].

<h2> 1. Replica exchange Monte Carlo </h2>

Usage:
<pre><code>$ python replica_exchange.py</pre></code>
Output is the trajectory of replica at the lowest temperature in the pdb format.


<h2> 2. Simulated annealing </h2>

Usage:
<pre><code>$ python simulated_annealing.py</pre></code>
Output is the trajectory in the pdb format.

<h2> 3. Visualization </h2>

You can visualizate trajectory in the pdb format using PyMOL:
<pre><code>$ pymol trajectory.pdb</pre></code>
and then in PyMOL console:

<pre><code>PyMOL> run show</pre></code>

<hr>

[1] Dill, KA, Bromberg, S, Yue, K, Fiebig, KM, Yee, DP, Thomas, PD, Chan, HS (1995). Principles of protein folding--a perspective from simple exact models. Protein Sci., 4, 4:561-602.

