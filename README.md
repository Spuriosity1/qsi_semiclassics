Quantum spin ice classical dynamics simulation, originally written by Atilla Szabó 2018
modified by Alaric Sanders in 2022.

A detailed description of the calculations used can be found [on arxiv](https://arxiv.org/abs/1902.08641v2)
Recommended reading: [Semiclassical approach to quantum spin ice](http://arxiv.org/abs/1609.03079)

# Compiling 

This project is built with [meson](https://mesonbuild.com/).
Build the project with the standard litany,
```bash
meson setup build
ninja -C build
```
To move all the executables where `driver` expects them to be (a directory called `bin` in the project root), run
```bash
ninja -C build install
```

# Running a Benchmark
On a simple machine, simply run
`parallel < static_corr.plan`

For SLURM, simply run
`sbatch runme.slurm`


# Indexing convention

The 3D arrangement of pyrochlore spins and dual pyrochlore sites (i.e. plaquettes) are by necessity serialised in all output files. The hierarchy to have in mind is
```
    Cubic unit cell -> FCC site[4] -> pyrochlore sublattice[4]
```
i.e. every unit cell has 16 sites' worth of information. The convention is to index the pyrochlores from the 'up' sublattice of the diamond lattice. There are four FCC sites in a cubic unit cell.
The cubic cell at cubic-lattice Bravais site `[x y z]` has a serialised index (for a given system size `L`) `[x y z] -> L*(L*x + y) + z`. The same convention applies for dual lattice sites, see [the visualisation here](https://spuriosity1.github.io/2022-03-22-diamondrender/)


