# Configuration file for the typical embryo-like network: hexagonal with 2 layers.
[files]
outdir = ../outputs/
tpl = example.ini
subdir = config
nsim = 1

[network]
topology = hexagonal
nlayers = 2
pbc = False
nbicellular = 0
seed = None

[tensions]
gamma_border = 0.90
gamma_c_border = 1.00

[noisy]
noisy = False
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