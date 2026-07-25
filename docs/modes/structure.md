# Model structure

The model structure defines the process representations used by a FUSE
simulation. Rather than embedding these choices within a fixed model,
FUSE represents the model structure as a series of modelling decisions,
each corresponding to one aspect of the hydrologic response.

The model structure is specified through the decision file, which is
identified by the `decisions_file` entry in the control file. The
decision file specifies one selected process representation for each
modelling decision, thereby defining a specific model instantiation.

The modelling decisions implemented in FUSE are described by
[Clark et al. (2008)](https://doi.org/10.1029/2007WR006735), with the
snow modelling decision (Decision 9) described by
[Henn et al. (2015)](https://doi.org/10.1002/2014WR016736).

The Bow River test case uses the decision file

```text
test/CAN_05BB001/settings/fuse_v2/fuse_zDecisions_2.txt
```

The decision file begins by specifying one selected process
representation for each modelling decision. These entries are followed
by a list of all available process representations, allowing users to
construct alternative model instantiations by selecting different
process representations.

The current version of FUSE defines nine modelling decisions.

| Decision | Description |
|----------|-------------|
| 1 | Rainfall error model |
| 2 | Upper-layer storage architecture |
| 3 | Lower-layer storage architecture and baseflow |
| 4 | Surface runoff formulation |
| 5 | Percolation |
| 6 | Evaporation |
| 7 | Interflow |
| 8 | Runoff routing |
| 9 | Snow model |

For example, the following entries define a model with

- multiplicative rainfall error;
- a single upper-layer tension store;
- two parallel lower-layer reservoirs;
- ARNO/VIC surface runoff;
- field-capacity-based percolation;
- sequential evaporation;
- no interflow;
- gamma-distribution routing; and
- a temperature-index snow model.

```text
multiplc_e  RFERR
tension1_1  ARCH1
tens2pll_2  ARCH2
arno_x_vic  QSURF
perc_f2sat  QPERC
sequential  ESOIL
intflwnone  QINTF
rout_gamma  Q_TDH
temp_index  SNOWM
```

The remainder of the decision file lists all available process
representations for each modelling decision. Alternative model
instantiations can be created by selecting different process
representations from these lists.
