"""
Beta-amyloid accumulation, neuroinflammation, and neurodegeneration (PySB).

A phenomenological model: all amounts are normalized (arbitrary units),
time is in YEARS. Rates are illustrative, not fitted to data.

Modules
  1. Abeta production, aggregation (monomer -> oligomer -> plaque)
  2. Clearance by microglia (impaired when activated or senescent)
  3. Inflammation: microglial activation, cytokine release, positive feedback
  4. Neurons: healthy -> stressed -> dead; dead-cell debris fuels inflammation
"""
from pysb import Model, Monomer, Parameter, Initial, Rule, Observable

Model()

# ---------------------------------------------------------------- monomers
Monomer('Abm')                                   # soluble Abeta monomer
Monomer('Abo')                                   # soluble oligomer (toxic)
Monomer('Plq')                                   # fibrillar plaque mass
Monomer('Mg', ['s'], {'s': ['rest', 'act', 'sen']})  # microglia
Monomer('Cyt')                                   # pro-inflammatory cytokines (IL-1b, TNF-a)
Monomer('N', ['s'], {'s': ['h', 'st']})          # neurons: healthy / stressed
Monomer('Deb')                                   # debris from dead neurons (DAMPs)
Monomer('Dead')                                  # cumulative neuron deaths (counter)

# --------------------------------------------------------------- parameters
# Abeta production & aggregation
Parameter('k_prod', 1.0)      # APP -> Abeta (beta/gamma-secretase), per year
Parameter('k_deg', 0.5)       # enzymatic degradation (neprilysin, IDE)
Parameter('k_olig', 1.0)      # 2 Abm -> Abo
Parameter('k_dis', 0.5)       # Abo -> 2 Abm
Parameter('k_nuc', 0.02)      # Abo -> plaque nucleation
Parameter('k_grow', 2.0)      # Abo + Plq -> 2 Plq (templated growth)

# Microglial clearance
Parameter('k_clr_m', 5.0)     # resting microglia clear monomer
Parameter('k_clr_o', 2.0)     # resting microglia clear oligomer
Parameter('f_act', 0.3)       # activated microglia clear at 30% efficiency
Parameter('k_clr_p', 0.05)    # activated microglia clear plaque (slow)

# Inflammation
Parameter('k_act_o', 2.0)     # oligomer activates microglia
Parameter('k_act_p', 0.5)     # plaque activates microglia
Parameter('k_act_c', 1.0)     # cytokine activates microglia (feedback)
Parameter('k_act_d', 1.0)     # neuron debris activates microglia
Parameter('k_res', 1.5)       # activated -> resting (resolution)
Parameter('k_sen', 0.03)      # resting -> senescent (aging)
Parameter('k_cyt_a', 1.0)     # cytokine release by activated microglia
Parameter('k_cyt_s', 0.3)     # cytokine release by senescent microglia (SASP)
Parameter('k_cyt_d', 5.0)     # cytokine decay

# Neurons
Parameter('k_str_o', 0.5)     # oligomer-induced neuronal stress
Parameter('k_str_c', 0.5)     # cytokine-induced neuronal stress
Parameter('k_rec', 1.0)       # stressed -> healthy recovery
Parameter('k_die', 0.1)       # stressed neuron death
Parameter('k_dcl', 2.0)       # debris cleared by resting microglia
Parameter('k_ddeg', 0.5)      # debris passive decay

# Initial conditions
Parameter('Mg_0', 1.0)
Parameter('N_0', 1.0)

Initial(Mg(s='rest'), Mg_0)
Initial(N(s='h'), N_0)

# -------------------------------------------------------------------- rules
# 1. Production, degradation, aggregation
Rule('Ab_production', None >> Abm(), k_prod)
Rule('Ab_degradation', Abm() >> None, k_deg)
Rule('oligomerization', Abm() + Abm() | Abo(), k_olig, k_dis)
Rule('plaque_nucleation', Abo() >> Plq(), k_nuc)
Rule('plaque_growth', Abo() + Plq() >> Plq() + Plq(), k_grow)

# 2. Microglial clearance (resting = efficient, activated = impaired)
from pysb import Expression
Expression('k_clr_m_act', k_clr_m * f_act)
Expression('k_clr_o_act', k_clr_o * f_act)
Rule('clear_Abm_rest', Mg(s='rest') + Abm() >> Mg(s='rest'), k_clr_m)
Rule('clear_Abo_rest', Mg(s='rest') + Abo() >> Mg(s='rest'), k_clr_o)
Rule('clear_Abm_act', Mg(s='act') + Abm() >> Mg(s='act'), k_clr_m_act)
Rule('clear_Abo_act', Mg(s='act') + Abo() >> Mg(s='act'), k_clr_o_act)
Rule('clear_Plq_act', Mg(s='act') + Plq() >> Mg(s='act'), k_clr_p)

# 3. Inflammation
Rule('activate_by_Abo', Mg(s='rest') + Abo() >> Mg(s='act') + Abo(), k_act_o)
Rule('activate_by_Plq', Mg(s='rest') + Plq() >> Mg(s='act') + Plq(), k_act_p)
Rule('activate_by_Cyt', Mg(s='rest') + Cyt() >> Mg(s='act') + Cyt(), k_act_c)
Rule('activate_by_Deb', Mg(s='rest') + Deb() >> Mg(s='act') + Deb(), k_act_d)
Rule('resolution', Mg(s='act') >> Mg(s='rest'), k_res)
Rule('senescence', Mg(s='rest') >> Mg(s='sen'), k_sen)
Rule('cytokine_act', Mg(s='act') >> Mg(s='act') + Cyt(), k_cyt_a)
Rule('cytokine_sen', Mg(s='sen') >> Mg(s='sen') + Cyt(), k_cyt_s)
Rule('cytokine_decay', Cyt() >> None, k_cyt_d)

# 4. Neurons
Rule('stress_by_Abo', N(s='h') + Abo() >> N(s='st') + Abo(), k_str_o)
Rule('stress_by_Cyt', N(s='h') + Cyt() >> N(s='st') + Cyt(), k_str_c)
Rule('recovery', N(s='st') >> N(s='h'), k_rec)
Rule('neuron_death', N(s='st') >> Deb() + Dead(), k_die)
Rule('debris_clear', Mg(s='rest') + Deb() >> Mg(s='rest'), k_dcl)
Rule('debris_decay', Deb() >> None, k_ddeg)

# -------------------------------------------------------------- observables
Observable('Ab_monomer', Abm())
Observable('Ab_oligomer', Abo())
Observable('Plaque', Plq())
Observable('Mg_resting', Mg(s='rest'))
Observable('Mg_active', Mg(s='act'))
Observable('Mg_senescent', Mg(s='sen'))
Observable('Cytokine', Cyt())
Observable('N_healthy', N(s='h'))
Observable('N_stressed', N(s='st'))
Observable('N_total', N())
Observable('N_dead', Dead())
