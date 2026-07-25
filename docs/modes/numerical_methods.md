# Numerical solution methods

The numerical solution methods define how the governing hydrologic model
equations are integrated through time. By separating the governing
equations from their numerical solution, FUSE allows the same model
structure and parameter set to be evaluated using different
time-stepping methods.

The numerical solution methods are specified through the numerics file,
which is identified by the `numerics_file` entry in the control file.

The Bow River test case uses the numerics file

```text
test/CAN_05BB001/settings/fuse_v2/fuse_zNumerix.txt
```

The numerical methods implemented in FUSE were developed and evaluated
by

> Clark, M. P., & Kavetski, D. (2010). Ancient numerical daemons of
> conceptual hydrological modeling: 1. Fidelity and efficiency of time
> stepping schemes. *Water Resources Research*, 46, W10510.
> https://doi.org/10.1029/2009WR008894

This work demonstrated that numerical approximation errors can be
comparable to, or even exceed, model structural errors when
inappropriate time-stepping schemes are used. It also identified
practical numerical methods that provide a good balance between
computational efficiency and numerical reliability.

## Numerics file entries

The numerics file specifies

- the numerical time-stepping method;
- whether fixed or adaptive time stepping is used;
- the absolute and relative error tolerances for adaptive time
  stepping; and
- additional solver settings.

The available time-stepping methods are

| Option | Description |
|--------|-------------|
| Explicit Euler | First-order explicit Euler method. |
| Explicit Heun | Second-order explicit Heun method. |
| Implicit Euler | First-order implicit Euler method. |
| Implicit Heun | Second-order implicit Heun method. |
| Semi-implicit Euler | First-order semi-implicit Euler method. |

## Physics branches

FUSE currently provides two physics branches.

- **fuse_v1** is the legacy implementation and supports all numerical
  methods defined in the numerics file, including both fixed-step and
  adaptive time-stepping schemes.

- **fuse_v2** is the differentiable implementation. To simplify the
  temporal propagation of parameter sensitivities and reduce software
  complexity, this branch currently uses a fixed-step Implicit Euler
  method and does not expose alternative numerical solution methods.

## Fixed and adaptive time stepping

In **fuse_v1**, each numerical method may be used with either fixed or
adaptive time stepping.

With **fixed time stepping**, the governing equations are integrated
using a prescribed step size throughout the simulation. Numerical
accuracy therefore depends on whether the selected time step is
appropriate for the model structure, parameter values, forcing data,
and evolving model states.

With **adaptive time stepping**, FUSE automatically subdivides the model
time step to satisfy specified absolute and relative error tolerances.
Substeps that exceed the prescribed error tolerances are repeated using
a smaller step size, while the step size is increased when the
estimated numerical error is sufficiently small.

Adaptive time stepping generally provides more consistent numerical
accuracy than selecting a uniformly small fixed time step.

## Solution constraints

Hydrologic state variables are subject to physical constraints. For
example, water storage cannot be negative and may be bounded by a
maximum storage capacity. FUSE checks the feasibility of the numerical
solution and applies constraint handling where required.

## Guidance

Clark and Kavetski (2010) showed that fixed-step explicit methods can
produce large and uncontrolled numerical errors over substantial regions
of the FUSE parameter space. Among the methods evaluated, the
fixed-step **Implicit Euler** method and the **adaptive Explicit Heun**
method provided the best overall balance between numerical fidelity and
computational efficiency for most hydrologic applications.

The numerical method and related solver settings are specified through
the numerics file. In **fuse_v1**, the control file selects the
numerical method and, when adaptive time stepping is used, specifies the
absolute and relative error tolerances. In **fuse_v2**, the numerical
solution is fixed to the Implicit Euler method and these options are not
currently exposed.
