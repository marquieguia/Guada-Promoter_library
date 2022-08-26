import pandas as pd
from glob import glob
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import string
import time
start = time.process_time()


# Auxiliary functions
def p(s):
    print(s)


def extractor(f):

    OD = pd.read_excel(f ,header=14 ,usecols='B:M' ,nrows=8).to_numpy()
    RFU = pd.read_excel(f ,header=25 ,usecols='B:M' ,nrows=8).to_numpy()

    flat_RFU = RFU.flatten()
    flat_OD = OD.flatten()
    basename = os.path.basename(f)
    conc_ind = basename[basename.find('.') + 1]

    return flat_RFU ,flat_OD ,conc_ind


def did_grow(df ,OD_thres):
    bool_df = df.ge(OD_thres).all(axis=0)
    grown_df = df.drop(labels=bool_df[bool_df == False].index ,axis=1)
    dead_df = df.drop(labels=bool_df[bool_df == True].index ,axis=1)
    return grown_df ,dead_df


def raise_error(error_dic, label):
    for stock, percent in error_dic.items():
        p('')
        p('Alarm for high ' + label + ' ratio!!!')
        if label == 'dead':
            p('For ' + stock + ' a total of ' + str(percent) + '% have died.')
            p('')
        elif label == 'down':
            p('For ' + stock + ' a total of ' + str(percent) + '% have downward slope.')
            p('')


# Parameters for thresholds
OD_thres = float(input('OD Threshold for dead wells'))  # OD value that determines that didn't grow
dead_thres = 20  # if more than 20% didn't grow raise alarm!
downward_thres = 40  # if more than 40% have downward lineplot trend, raise alarm!

# Reading stocks files
stocks = []
files = glob(os.path.join('..', 'Sample_data' ,'*.xlsx'))
for f in files:
    base = os.path.basename(f)
    ind = base.find('.')
    stock = base[:ind]
    stocks.append(stock)

stocks = sorted(set(stocks))
stocks = [str(stock) for stock in stocks]
plate_dimensions = (8, 12)
rows = [char for char in string.ascii_uppercase[:plate_dimensions[0]]]
cols = [x+1 for x in range(plate_dimensions[1])]

# Auxiliary dictionaries to store across stocks
mother_dic = {}
dead_dic = {}
downward_dic = {}
dead_errors = {}
down_errors = {}
error_summary = np.ones((1, 2))
for stock in stocks:
    '''stock refers to the first two digits of the plate id,
       while conc refers to the digit following the plate name
       e.g. for 42.1.xlsx the stock == 42 and conc == 1'''

    columns_names = [stock + row + str(col) for row in rows for col in cols]

    ##### One stock at the time #####

    OD_dic = {}  # dic stores ODs of stock across the concentrations (ie plates)
    RFU_dic = {}  # dic stores RFUs of stock across the concentrations (ie plates)

    concs = sorted([f for f in files if os.path.basename(f).startswith(stock)])

    for conc in concs:
        ##### One concentration (ie plate) at the time ######
        flat_RFU ,flat_OD ,conc_ind = extractor(conc)
        OD_dic[conc_ind] = flat_OD
        RFU_dic[conc_ind] = flat_RFU

    ##### One stock at the time #####
    stock_OD = np.vstack((OD_dic['1'] ,OD_dic['2'] ,OD_dic['3'] ,OD_dic['4']))
    stock_RFU = np.vstack((RFU_dic['1'] ,RFU_dic['2'] ,RFU_dic['3'] ,RFU_dic['4']))
    stock_RFU_div_OD = stock_RFU / stock_OD
    stock_RFU_and_OD = np.vstack((stock_OD ,stock_RFU ,stock_RFU_div_OD))
    stock_df = pd.DataFrame(stock_RFU_and_OD ,columns=columns_names ,dtype=np.float64 ,
                            index=['OD.1' ,'OD.2' ,'OD.3' ,'OD.4' ,
                                   'RFU.1' ,'RFU.2' ,'RFU.3' ,'RFU.4' ,
                                   'RFU.1/OD.1' ,'RFU.2/OD.2' ,'RFU.3/OD.3' ,'RFU.4/OD.4'])
    # check how many grew
    grown_stock_df ,dead_stock_df = did_grow(stock_df ,OD_thres)
    dead_ratio = round((dead_stock_df.shape[1]/stock_df.shape[1])*100, 2)


    if dead_ratio >= dead_thres:
        dead_errors[stock] = round(dead_ratio)

    else:
        pass

    # check how many have downward slope (ie conc[x] < conc[x-1])
    breathing_room = 200
    RFU_div_OD_diff = grown_stock_df.loc[['RFU.1/OD.1' ,'RFU.2/OD.2' ,'RFU.3/OD.3' ,'RFU.4/OD.4']].diff()
    RFU_div_OD_diff += breathing_room

    bool_df = RFU_div_OD_diff.le(0).any(axis=0)
    downward_df = grown_stock_df.drop(labels=bool_df[bool_df == False].index ,axis=1)
    down_ratio = round((downward_df.shape[1]/grown_stock_df.shape[1])*100, 2)
    if down_ratio > downward_thres:
        down_errors[stock] = down_ratio

    else:
        pass

    error_summary = np.vstack((error_summary, (dead_ratio, down_ratio)))

    mother_dic[stock] = grown_stock_df
    dead_dic[stock] = dead_stock_df
    downward_dic[stock] = downward_df


# raise potential errors and error summary
raise_error(dead_errors, 'dead')
raise_error(down_errors, 'down')
error_summary = error_summary[1:, :]
error_df = pd.DataFrame(data=error_summary, dtype=np.float64, index=stocks, columns=['Dead (%)', 'Down (%)'])
p('ERROR SUMMARY:')
p(error_df)
error_summary_path = os.path.join('..', 'Sample_output', 'hit_data_error_summary.xlsx')
error_df.to_excel(error_summary_path)


# merge all stocks together and save
mother = pd.concat(mother_dic.values() ,axis=1)
dead = pd.concat(dead_dic.values() ,axis=1)
downward = pd.concat(downward_dic.values() ,axis=1)

# calculate fold change
mother.loc['Fold'] = mother.loc['RFU.4/OD.4'] / mother.loc['RFU.1/OD.1']

filter = (mother.loc['Fold'] > 3.0)
top_df = mother.loc[: , filter]

# export to excel
top_df.to_excel(os.path.join('..', 'Sample_output' ,'top_hits.xlsx'))
mother.to_excel(os.path.join('..', 'Sample_output','hit_data_mother.xlsx'))
dead.to_excel(os.path.join('..', 'Sample_output','dead.xlsx'))
downward.to_excel(os.path.join('..', 'Sample_output','downward.xlsx'))
TPA_conc = list(np.log10([0.0001 ,1 ,25 ,1000]))


# quick plotting
fold = top_df.loc['Fold']
top = round(fold.nlargest(50, keep='all'), 3)
fig2, (ax1, ax2) = plt.subplots(2,1,figsize=(8,12))
fig2.subplots_adjust(hspace=0.5)
sns.lineplot(data=top_df, x=[legend=False, ci=None, ax=ax1)
ax1.set_ylabel('RFU/OD')
ax1.set_xlabel('log10([TPA])')
ax2.hist(fold, bins=20)
ax2.set_xlim((0, top.max()+0.5))
ax2.set_ylabel('Frequency')
ax2.set_xlabel('Fold change')
plt.show()
plt.savefig(os.path.join('..', 'Sample_output', 'Figs', 'Fold_ge_3.png'))
plt.close()

p('')
print('Time taken:')
print(time.process_time() - start)



