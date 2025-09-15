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
Via SLURM, the toolchain is
```bash
cd driver
./plan_static_corr.sh
# Importing temperature points from file 'temps'
# wrote static_corr.plan
./generate_slurmfile.sh static_corr.plan
# Wrote job script tmp/1757944459.slurm
# Run with
# sbatch --array=1-304 --cpus-per-task=1 tmp/1757944459.slurm
sbatch --array=1-304 --cpus-per-task=1 tmp/1757944459.slurm
```

To scale this up/down, modify the system_size and n_sample in input/SSF/test.toml:

```python
system_size=    4   # <--- linear dimension execution time scales like system_size^3
n_burnin=       4
n_sample=       24  # <--- bigger is better
```

# Indexing convention

The 3D arrangement of pyrochlore spins and dual pyrochlore sites (i.e. plaquettes) are by necessity serialised in all output files. The hierarchy to have in mind is
```
    Cubic unit cell -> FCC site[4] -> pyrochlore sublattice[4]
```
i.e. every unit cell has 16 sites' worth of information. The convention is to index the pyrochlores from the 'up' sublattice of the diamond lattice. There are four FCC sites in a cubic unit cell.
The cubic cell at cubic-lattice Bravais site `[x y z]` has a serialised index (for a given system size `L`) `[x y z] -> L*(L*x + y) + z`. The same convention applies for dual lattice sites, see [the visualisation here](https://spuriosity1.github.io/2022-03-22-diamondrender/)


