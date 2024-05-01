# NMR_scripts
Scripts for processing NMR data

# Error analysis for CPMG data

CPMG experiments often have duplicates performed at a few `ncyc` values. We find the best practice for estimating error to be to use

max(error from duplicates, error from spectrum noise)

Workflow:

1. Do processing in PINT using -noiseUncertainty flag to return error estimates from PINT that represent uncertainty 

-NB: The current PINT duplicate error processing uses pooled standard deviation across all sets of duplicates at the level of peak volumes / (intensities). However, using pooled standard deviation assumes that the variance across each set of duplicates is the same, which is not necessarily the case.

2. Use `process_CPMG_err.py` or `process_CPMG_err.ipynb` to calculate uncertainties for R_2,eff using duplicates, determine whether uncertainty from noise or duplicates are greater, and plot.

Simplest usage:

`python process_CPMG_err.py /path/to/PINT/out/ -T 0.04`

Would process all *.out files in PINT `out` dir, using a cpmg delay of 40 ms.

Outputs:

- Plots using largest error in `/output_plots/*pdf` using largest error

- Plots comparing error at `/output_plots/both_errs/*pdf`

- JSON file with raw data for further use in python: `output_plots/raw_data.json.zip`

NB: A useful functionality of PINT is to fit R_2,eff to the Carver-Richards equations to estimate p_B, kex, etc. These scripts do not currently do this.
