# PIN1 auxin CUC model

Model consisting on a cell grid to study the interactions between the plant hormone auxin, the PIN1 auxin efflux carrier, and the CUC transciption factors.

The cell grid represents the margin epidermis of a developing leaf (middle domain flanked by abaxial and adaxial domains) during the satages when patterning of auxin sites occur.

## How to use it

Set parameters in file `params.py`. See comment inside the file for explanation of each parameter.

3 possible ways to run the model:

* Single
A single simulations that reads the parameter values from `params.py`

* Series
A series of simuations that share the same parameters except for one that will vary within a linear space defined in variable `params.series_param_a`.
Varying the value one parameter...

* Batch
Batch of simulations that share the same parameters except those defined in the dictionary `params.batch_params`. Example to compare the effect of no local auxin creation, auxin creation in the CUC domain, and auxi creation in the middle domain:
```
batch_params = {
    1: {
        'k_auxin_synth': .8,
        'k_cuc_auxin_synth': 0,
        'k_md_auxin_synth': 0  
    },
    2: {
        'k_auxin_synth': .8,
        'k_cuc_auxin_synth': 0.25,
        'k_md_auxin_synth': 0  
    },
    3: {
        'k_auxin_synth': .8,
        'k_cuc_auxin_synth': 0,
        'k_md_auxin_synth': 0.01
    }
}
```
To simulate, run `run_model.py`.

To view simulation results, browse the folders 
* images/
* graphs/

To reveiew simulation parameters / rerun a simulation:
* sim_logs/

