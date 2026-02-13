This project provides a minimal reproducible example for an issue I encountered when using the composite linear solver in AMReX. I found that when solving the equation

$$ \nabla \cdot (\sigma \nabla \phi) = 0 , $$

for a certain (challenging) choice of $\sigma$ and boundary conditions, the solver performs very differently depending on whether the simulation is run with or without AMR - even if the AMR is applied across the entire domain and the effective resolution is the same as in the unigrid case.

To see this in action, set the `AMREX_HOME` and `HYPRE_DIR` variables in `GNUmakefile` and compile using `make`. Then, run using `./main2d.gnu.MPI.ex settings`. The code will first solve the problem on a $128 \times 128$ mesh without AMR, and then on a $32 \times 32$ mesh with two levels of AMR which are applied across the entire domain, bringing the effective resolution up to $128 \times 128$ cells.

You should see that the first simulation converges within only three iterations, while the second simulation takes 10877 iterations to complete. The output data can be found in `test_unigrid_plot` and `test_AMR_plot`. The final results are very similar, although the AMR results exhibits some banding for $y < 0$ that becomes visible when plotting `phi` on a logarithmic scale.

Now, my question is - why does this happen? Why is the performance so different between these two runs?