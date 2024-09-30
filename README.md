# PRAS_HydroDispatchModifier

## Installation

This package is not on any julia registry so the way to add it is in the package manager using 
```
]add https://github.com/NREL/PRAS_HydroDispatchModifier.git
```

## Usage

Basic usage to get a modified PRAS SystemModel is

```
new_sys = modify_hydro(base_sys,hydro_data,hydro_flex)
```

`hydro_data` is a DataFrame of a long form table containing hydropower data including weekly or monthly generation, and minimum & maximum dispatch capacity. An [example input file](test/rts_hdgens.csv) can be found in the [test folder](test). It is up to the user to ensure that there is no overlap between weekly generators and monthly generators, as the code will throw an error otherwise.

`hydro_flex` can be one of either:
  - :high (represent reservoir units as GeneratorStorages)
  - :mid (same as above, but enforce min gen via capacity, energy,
          and regional net load adjustment)
  - :low (dispatch hydro at a fixed level consistent with energy and power constraints)

