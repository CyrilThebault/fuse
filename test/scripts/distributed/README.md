# Prepare a mizuRoute domain

The recommended workflow is to use `prepare_mizuroute_domain.sh`, which

1. creates a complete hydrologic fabric from MERIT-Basins shapefiles;
2. prepares a compact mizuRoute hydrofabric; and
3. creates a runoff-remapping file for lumped hydrology with distributed
   routing.

The third step supports the test configuration in which the hydrologic model
produces a single basin-average runoff time series while mizuRoute routes water
through a distributed river network. The mapping file assigns the single runoff
cell to every routing HRU.

The complete preprocessing workflow is

```text
MERIT-Basins shapefiles
           │
           ▼
create_hydrofabric.R
           │
           ▼
hydrofabric_merit.nc
           │
           ▼
prepare_mizuroute_hydrofabric.sh
           │
           ▼
hydrofabric_mizuRoute.nc
           │
           ▼
create_lumped_to_hru_mapping.sh
           │
           ▼
lumped_to_hru.nc
```

## Example (Bow River distributed-routing test case)

```bash
catchment_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_basin.shp
river_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_river.shp
metadata=test/metadata/merit_basins_shapefile_metadata.csv
hydrofabric_merit=test/CAN_05BB001/distributed/input/hydrofabric_merit.nc
hydrofabric_mizuRoute=test/CAN_05BB001/distributed/input/hydrofabric_mizuRoute.nc
mapping_file=test/CAN_05BB001/distributed/input/lumped_to_hru.nc
```

Run the complete workflow with

```bash
test/scripts/distributed/prepare_mizuroute_domain.sh \
    "${river_shp}" \
    "${catchment_shp}" \
    "${metadata}" \
    "${hydrofabric_merit}" \
    "${hydrofabric_mizuRoute}" \
    "${mapping_file}"
```

This workflow creates

```text
test/CAN_05BB001/distributed/input/hydrofabric_merit.nc
test/CAN_05BB001/distributed/input/hydrofabric_mizuRoute.nc
test/CAN_05BB001/distributed/input/lumped_to_hru.nc
```

Temporary files are written to

```text
test/CAN_05BB001/distributed/input/work/
```

during processing and are removed automatically when the workflow completes
successfully.

---

# Lumped hydrology with distributed routing

The runoff input contains a single spatial cell, whereas the routing network
contains multiple river-network HRUs. The runoff-remapping file therefore
assigns every routing HRU to runoff cell `(1,1)`:

```text
polyid      = hruId
nOverlaps   = 1
weight      = 1.0
i_index     = 1
j_index     = 1
```

Internally, mizuRoute stores the remapping information as a ragged array.

The variables

- `polyid`
- `nOverlaps`

contain one entry for each routing HRU.

The variables

- `weight`
- `i_index`
- `j_index`

contain one entry for every runoff-cell/HRU overlap. Consequently,

```text
size(data) = sum(nOverlaps)
```

For the lumped-routing configuration, every routing HRU overlaps the single
runoff cell exactly once. Therefore,

```text
nOverlaps(:) = 1
size(data) = number of routing HRUs
```

This configuration is appropriate when

- the runoff NetCDF file contains exactly one spatial cell;
- the runoff variable represents a basin-average runoff depth rather than a
  basin-total discharge;
- the routing HRUs collectively represent the same basin; and
- the routing HRU areas sum to the basin area represented by the lumped
  hydrologic model.

mizuRoute applies the same runoff depth to each routing HRU and converts that
depth to a local runoff volume using the HRU area. Consequently, the basin-total
runoff volume is preserved while allowing runoff to enter the distributed river
network at every routing HRU.

---

# Running the individual scripts

## Step 1. Create the complete hydrologic fabric

```bash
Rscript create_hydrofabric.R \
    "${catchment_shp}" \
    "${river_shp}" \
    "${metadata}" \
    "${hydrofabric_merit}"
```

## Step 2. Prepare the compact mizuRoute hydrofabric

```bash
prepare_mizuroute_hydrofabric.sh \
    "${hydrofabric_merit}" \
    "${hydrofabric_mizuRoute}"
```

## Step 3. Create the lumped-to-HRU runoff mapping

```bash
create_lumped_to_hru_mapping.sh \
    "${hydrofabric_mizuRoute}" \
    "${mapping_file}"
```

---

# Script summary

| Script | Purpose |
|---------|---------|
| `prepare_mizuroute_domain.sh` | Runs the complete preprocessing workflow (recommended). |
| `create_hydrofabric.R` | Creates the complete hydrologic fabric from the MERIT-Basins river and catchment shapefiles. |
| `prepare_mizuroute_hydrofabric.sh` | Creates the compact mizuRoute hydrofabric by extracting the variables required by mizuRoute and converting them to the required units. |
| `create_lumped_to_hru_mapping.sh` | Creates the mizuRoute runoff-remapping file by assigning every routing HRU to the single lumped runoff cell. |

---

# Requirements

The complete workflow requires

- R with the `sf` and `ncdf4` packages; and
- the NCO utilities
  - `ncap2`
  - `ncks`
  - `ncrename`

The shell scripts must be executable:

```bash
chmod +x \
    prepare_mizuroute_domain.sh \
    prepare_mizuroute_hydrofabric.sh \
    create_lumped_to_hru_mapping.sh
```
