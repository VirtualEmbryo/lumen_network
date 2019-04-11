# Configuration file for a chain network, N=100, with PBC.
[files]
outdir = ../outputs/
tpl = example.ini
subdir = config
nsim = 1

[network]
topology = chain
nlayers = 100
pbc = True
nbicellular = 0
seed = 0

[tensions]
gamma_border = 1.00
gamma_c_border = 1.00

[noisy]
noisy = False			# Not implemented for chain network
lumen_pos_avg = 0.
lumen_pos_std = 0.1

[volume]
vol_avg = 1
vol_std = 0.1

[chimera]
chimeras = False
gamma_chimera = 0.5

[display]
show = False