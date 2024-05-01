import numpy as np
import pandas as pd
import sys, os
from glob import glob

from argparse import ArgumentParser

# for plotting
import matplotlib.pyplot as plt

def read_PINT(file_dir, T):
    ''' Reads in files from PINT.
    Inputs:
    file_dir: path to directory containing *.out files from PINT.
    T: CPMG delay (in seconds)
    
    Outputs: 
    dataframe: dataframe with one row per peak dataset (assi)
    '''
    
    fils = sorted(glob(file_dir+'/*out'))

    output = []

    for fil in fils:
        d = np.loadtxt(fil)
        # sort by first column to get ncyc smallest to largest
        d = d[d[:, 0].argsort()]

        ncyc, intensity, intensity_err, volume, volume_err = d[:,0], d[:,1], d[:,2], d[:,3], d[:,4]
        field = ncyc/T
        
        assi = os.path.basename(fil).split('.')[0]
        
        output.append({'assi': assi, 'fields_inc_ref':field, 'intensity': intensity, 'intensity_err': intensity_err,
                       'volume': volume, 'volume_err': volume_err})
        
    dataframe = pd.DataFrame.from_records(output)      
    return dataframe

def check_has_duplicates(row):
    '''
    Check dataset has duplicates
    '''
    tmp = pd.DataFrame({k: row[k] for k in ['intensity','fields_inc_ref']})

    #get subset with duplicates
    tmp['n_duplicates'] = tmp.groupby('fields_inc_ref')['intensity'].transform('nunique')
    dup_subset = tmp.loc[tmp.n_duplicates>1]
    
    if len(dup_subset) > 1:
        return True
    else:
        return False
    
def get_largest_error(row,xaxis='ncyc', meas='volume',scale=False):
    '''
    Calculates error using max. Mean Absolute Deviation over duplicate sets
    Then selects which error to use based on which is larger, error from noise or error from duplicates
    
    Inputs:
    row: row of dataframe
    meas: input to use: volume, intensity, R2eff
    
    Outputs:
    Intensity_err_from_dup: vector of calculated error values on intensities
    Vol_err_from_dup: vector of calculated error values on volumes
    '''
    
    #convert row to its own dataframe
    tmp = pd.DataFrame({k: row[k] for k in [meas,meas+'_err',xaxis]})

    #get subset with duplicates
    tmp['n_duplicates'] = tmp.groupby(xaxis)[meas].transform('nunique')
    dup_subset = tmp.loc[tmp.n_duplicates>1]
    
    if len(dup_subset)==0:
        raise RuntimeError('this data does not appear to have duplicates.')

    dup_subset[meas+'_mean_at_arr'] = dup_subset.groupby(xaxis)[meas].transform('mean')
    
    if scale:
        dup_subset[meas+'_abs_dev_scaled'] = dup_subset.apply(lambda row: np.abs(row[meas] - row[meas+'_mean_at_arr'])/row[meas], axis=1)
        tmp[meas+'_MAD'] = tmp.apply(lambda row: row[meas] * np.max(dup_subset[meas+'_abs_dev_scaled']), axis=1)
    else:
        dup_subset[meas+'_abs_dev'] = dup_subset.apply(lambda row: np.abs(row[meas] - row[meas+'_mean_at_arr']), axis=1)
        tmp[meas+'_MAD'] = np.max(dup_subset[meas+'_abs_dev'])

    # select which error to use
    bigger_error = ['duplicates','noise'][np.argmax([tmp[meas+'_MAD'].mean(), tmp[meas+'_err'].mean()])]

    return tmp[meas+'_MAD'].values, bigger_error    
        
def calculate_R2eff(row, T, value='volume'):
    '''
    Calculates R2_eff either based on PINT-calculated peak intensity or volume.
    
    Inputs:
    row: row of dataframe
    T: delay (in s)
    value: column to use to calculate (`intensity` or `volume`)
    
    Outputs:
    R2eff: vector of calculated R2eff values
    R2eff_err: propagated R2eff error
    ncyc: values of ncyc to plot (not including references of 0)
    '''
    
    ref_inds = np.where(row['fields_inc_ref']==0)
    arr_inds = np.where(row['fields_inc_ref']!=0)
        
    I0 = np.mean(row[value][ref_inds])
    I0_err = np.sqrt(np.mean(np.square(row[value+'_err'][ref_inds])))
    
    I= row[value][arr_inds]
    I_err = row[value+'_err'][arr_inds]
    
    R2eff = np.log(I0/I)/T
    R2eff_err = np.sqrt((I0_err/I0)**2 + (I_err/I)**2) * R2eff
    
    return R2eff, R2eff_err, row['fields_inc_ref'][arr_inds]

def estimate_Rex(row,n_endpoints=2):
    R2eff_init = row['R2eff'][0]
    R2eff_final = np.mean(row['R2eff'][-n_endpoints:])
    Rex = R2eff_init - R2eff_final

    if row['larger_err'] == 'duplicates':
        R2eff_init_err = row['R2eff_dup_err'][0]
        R2eff_final_err = np.sqrt(np.sum(np.square(row['R2eff_dup_err'][-n_endpoints:])))

    else:
        R2eff_init_err = row['R2eff_err'][0]
        R2eff_final_err = np.sqrt(np.sum(np.square(row['R2eff_err'][-n_endpoints:])))

    Rex_err = np.sqrt(R2eff_init_err**2+R2eff_final_err**2)
    return Rex, Rex_err, R2eff_final, R2eff_final_err    

def plot_peak(row):
    
    os.makedirs('output_plots',exist_ok=True)
    plt.figure(figsize=(4,3))
    if row['larger_err']=='duplicates':
        plt.errorbar(row['ncyc'], row['R2eff'],yerr=row['R2eff_dup_err'],capsize=2,fmt='.')
    else:
        plt.errorbar(row['ncyc'], row['R2eff'],yerr=row['R2eff_err'],capsize=2,fmt='.')
    plt.xlabel(r'$\nu_{CPMG}$ (Hz)')
    plt.ylabel(r'$R_{2,eff}$')

    plt.axhline(row['R2eff_final'],linestyle=':',zorder=0,color='grey')
    ax = plt.gca()
    x0, x1 = ax.get_xlim()
    plt.fill_between([x0,x1],row['R2eff_final']-row['R2eff_final_err'],\
        row['R2eff_final']+row['R2eff_final_err'],alpha=0.1, zorder=0, color='grey',linewidth=0)

    ax.set_xlim([x0,x1])

    plt.title('%s, Rex = %.2f ± %.2f /s\nerr method: %s' % (row['assi'], row['Rex'], row['Rex_err'], row['larger_err']))
    plt.savefig('output_plots/%s.pdf' % row['assi'],bbox_inches='tight')
    plt.close()
    
def plot_both_errs(row):
    os.makedirs('output_plots/both_errs',exist_ok=True)
    plt.figure(figsize=(4,3))
    plt.errorbar(row['ncyc']-10, row['R2eff'],yerr=row['R2eff_err'],capsize=2,fmt='.',label='noise err')
    plt.errorbar(row['ncyc']+10, row['R2eff'],yerr=row['R2eff_dup_err'],capsize=2,fmt='.',label='duplicate err')
    plt.legend()
    plt.xlabel(r'$\nu_{CPMG}$ (Hz)')
    plt.ylabel(r'$R_{2,eff}$')
    plt.title('%s\nerr method: %s' % (row['assi'], row['larger_err']))
    plt.savefig('output_plots/both_errs/%s.pdf' % row['assi'],bbox_inches='tight')
    plt.close()
    
def write_for_chemex(df, filename='test'):
    '''Write R2eff from dataframe format to input format for fitting in ChemEx
    
    Inputs:
    df: dataframe
    filename: name of file to write
    '''
    
    os.makedirs('chemex_input_files',exist_ok=True)

    dat[['volume_dup_err','larger_vol_err']] = dat.apply(lambda row: get_largest_error(row,meas='volume',xaxis='fields_inc_ref'), axis=1,result_type='expand')

    for _, row in df.iterrows():
        with open('chemex_input_files/%s.out' % row['assi'],'w') as f:
            f.write("#CPMG delay\tVolume\tVolume_err_%s\n"%row['larger_vol_err'])
            for i, val in enumerate(row['fields_inc_ref']):
                if row['larger_vol_err'] == 'duplicates':
                    f.write('%.4E\t%.4E\t%.4E\n' % (val, row['volume'][i], row['volume_dup_err'][i]))
                else:
                    f.write('%.4E\t%.4E\t%.4E\n' % (val, row['volume'][i], row['volume_err'][i]))

if __name__=='__main__':

    p = ArgumentParser()
    p.add_argument("dir_path",
                   help="path to output files, i.e. 'out'")
    p.add_argument("-T", type=float,
                   help="cpmg delay in seconds (i.e., 0.04 for 40 ms)")
    p.add_argument("--write_chemex_input", action='store_true')

    args = p.parse_args()

    dat = read_PINT(args.dir_path,T=args.T)

    print('Found %d peak datasets, calculating R2eff...' % len(dat))

    #first calculate R2eff and propagate R2eff_err from noise.
    dat[['R2eff','R2eff_err','ncyc']] = dat.apply(lambda row: calculate_R2eff(row,T=args.T, value='volume'), axis=1,result_type='expand')

    has_duplicate = check_has_duplicates(dat.iloc[0])
    if has_duplicate:
        print('Dataset has duplicates. Identifying largest error estimate...')
        #then determine error from R2eff duplicates.
        dat[['R2eff_dup_err','larger_err']] = dat.apply(lambda row: get_largest_error(row,meas='R2eff',xaxis='ncyc'), axis=1,result_type='expand')
    else:
        print('no duplicates found.')

    # estimate Rex based on difference in beginning and end points.
    dat[['Rex','Rex_err','R2eff_final', 'R2eff_final_err']] = dat.apply(lambda row: estimate_Rex(row), axis=1,result_type='expand')

    # save Rex values to csv
    dat[['assi','larger_err','Rex','Rex_err']].to_csv('rex_vals.csv',index=False)
    print('Wrote R_ex values and uncertainties to `rex_vals.csv`')


    #create plots
    print('Creating plots in output_plots/*pdf ...')
    dat.apply(lambda row: plot_peak(row), axis=1)
    print('Creating plots with both error types to compare in output_plots/both_errs/*pdf ...')
    dat.apply(lambda row: plot_both_errs(row), axis=1)

    dat.to_json('output_plots/raw_data.json.zip')
    print('')
    print('Raw data is in output_plots/raw_data.json.zip')
    print("Reload this in python with `df = pd.read_json('raw_data.json.zip')`")

    if args.write_chemex_input:
        write_for_chemex(dat)
        print('Wrote ChemEx input files at `chemex_input_files/*out`')
    print('')
    print('How many peaks used duplicates vs. noise?')
    print(dat.groupby('larger_err').size().to_string())




