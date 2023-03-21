# PIN1 auxin CUC model

Model consisting on a cell grid to study the interactions between the plant hormone auxin, the PIN1 auxin efflux carrier, and the CUC transciption factors.

The cell grid represents the margin epidermis of a developing leaf (middle domain flanked by abaxial and adaxial domains) during the satages when patterning of auxin sites occur.

## How to use it

Tested with Python 3.8.2, numpy X.X.X, pandas X.X.X, seaborn X.X.X, matplotlib X.X.X, scypi X.X.X, shutil X.X.X

Set parameters in `params_template.py` or an identically formated .py file copy. See comment inside the file for explanation of each parameter.

Run model:
```
python3 run_model.py params.py
```

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

