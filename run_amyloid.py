"""Simulate the amyloid/inflammation model over 60 years and plot results."""
import os

# Make Perl and BioNetGen visible regardless of how the IDE sets up PATH
PERL_BIN = r'C:\Users\corbi\strawberry-perl\perl\bin'
BNG_DIR = r'C:\Users\corbi\BioNetGen-2.9.3'
os.environ['PATH'] = PERL_BIN + os.pathsep + os.environ.get('PATH', '')
os.environ['BNGPATH'] = BNG_DIR

import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from amyloid_model import model

years = np.linspace(0, 60, 601)
sim = ScipyOdeSimulator(model, tspan=years, compiler='python')
res = sim.run()
df = res.dataframe

fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

ax = axes[0]
ax.plot(years, df['Ab_monomer'], label='Abeta monomer')
ax.plot(years, df['Ab_oligomer'], label='Abeta oligomer (toxic)')
ax.plot(years, df['Plaque'], label='Plaque')
ax.set_ylabel('Amount (a.u.)')
ax.set_title('Amyloid-beta')
ax.legend()

ax = axes[1]
ax.plot(years, df['Mg_resting'], label='Microglia: resting')
ax.plot(years, df['Mg_active'], label='Microglia: activated')
ax.plot(years, df['Mg_senescent'], label='Microglia: senescent')
ax.plot(years, df['Cytokine'], '--', label='Cytokines')
ax.set_ylabel('Fraction / amount')
ax.set_title('Neuroinflammation')
ax.legend()

ax = axes[2]
ax.plot(years, df['N_total'], label='Living neurons')
ax.plot(years, df['N_stressed'], label='Stressed neurons')
ax.plot(years, df['N_dead'], label='Cumulative deaths')
ax.set_ylabel('Fraction of initial')
ax.set_xlabel('Years from simulation start')
ax.set_title('Neurons')
ax.legend()

fig.tight_layout()
fig.savefig('amyloid_results.png', dpi=120)

for t in (0, 10, 20, 30, 40, 50, 60):
    row = df.iloc[int(t * 10)]
    print(f"yr {t:2d}: plaque={row['Plaque']:.3f} olig={row['Ab_oligomer']:.3f} "
          f"Mg_act={row['Mg_active']:.2f} Mg_sen={row['Mg_senescent']:.2f} "
          f"cyt={row['Cytokine']:.3f} neurons={row['N_total']:.3f}")

plt.show()