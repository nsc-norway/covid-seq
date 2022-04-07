#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme()

try:
    os.mkdir("plots")
except OSError as e:
    if e.errno == 17: pass
    else:             raise e

# make the plots
samples = pd.read_csv("detailed_variant_table.csv")
variants = pd.read_csv("variants.csv")

def plot_sample(sample_data):
    sample_name = sample_data['Sample'].iloc[0]

    # First make a data frame twice as long (long format)
    sample_data_total = sample_data.copy()
    sample_data_total['Depth'] = sample_data_total.total_depth - sample_data_total.alt_depth
    sample_data_total['Type'] = 'Other'
    sample_data.rename(columns={'alt_depth': 'Depth'}, inplace=True)
    sample_data['Type'] = 'Variant (Alt)'
    df = pd.concat((sample_data, sample_data_total))
    fig = plt.figure()
    sns.barplot(data=sample_data_total, x='Variant', y='Depth', color='grey', label='All reads')
    sns.barplot(data=sample_data, x='Variant', y='Depth', color='red', label='Reads with variant')
    #df.plot(x='Variant', y='Depth', kind='bar', stacked=True)
    plt.xticks(rotation=90)
    plt.title("Sample {}".format(sample_name))
    plt.legend()

    plt.tight_layout()
    plt.savefig('plots/{}.pdf'.format(sample_name))
    plt.close()


samples.groupby(by="Sample").apply(plot_sample)
