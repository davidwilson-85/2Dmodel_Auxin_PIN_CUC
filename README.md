# PIN1 auxin CUC model

Model consisting on a cell grid to study the interactions between the plant hormone auxin, the PIN1 auxin efflux carrier, and the CUC transciption factors.

This model was used in [HERE REFERENCE TO PAPER] to study the patterning of auxin foci in the margin of plant leaves.

Tested with Python 3.8.2, scypi 1.3.3, numpy 1.23.3, pandas 1.2.3, matplotlib 3.4.0, seaborn 0.11.1, pillow 7.0.0.

## How to run models

Set parameters in `params.py` or an identically formated `.py` file copy. See comments inside the file for explanation of each parameter.

Run model:
```
python3 run_model.py params_file.py
```

There are 3 ways to run the model:

### Single

A single simulation that reads the parameter values from `params.py`

### Series

A series of simuations that share the same parameters except for one that will vary within a linear space. In the parameters file, set `is_series` to `True` and define the variable and the min and max values in the dictionary `series_param_a`. Example:

```
series_param_a = {
    'name': 'k_auxin_diffusion',
	'min': .1,
	'max': 1,
	'num_points': 11
```

A series can be also used to run replicates, for example when Auxin or CUC noise is used. Example to run 10 replicates:

```
series_param_a = {
    'name': 'dummy',
	'min': 0,
	'max': 0,
	'num_points': 10
```

### Batch

A batch of simulations that share the same parameters except those defined as follows: In the parameters file, set `is_batch` to `True` and define parameters and values in the dictionary `batch_params`. The informatoin in `batch_params` overrides the values specified normally in `params.py`.

Example to compare the effect of auxin creation in the CUC domain versus in the middle domain (dfined in `template_middle_domain`):

```
batch_params = {
    'A': {
        'k_cuc_auxin_synth': 0.25,
        'k_md_auxin_synth': 0  
    },
    'B': {
        'k_cuc_auxin_synth': 0,
        'k_md_auxin_synth': 0.01
    }
}
```

Batch and Series can be combined. In such case, a batch is hierarchiucally above the series.

## View output of simulations

To view simulation results, browse the folders 
* `images/`
* `graphs/`

If `is_batch` is set to `True`, the results appear in 

* `out_batch/`

To review simulation parameters / re-run a simulation:

* `sim_logs/`

