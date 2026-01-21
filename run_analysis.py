from src import AlphaGenomePredictor, AlphaGenomePlotter
import config
from alphagenome.data import genome
from alphagenome.models import dna_client
import matplotlib.pyplot as plt
import numpy as np
import pysam
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr


model = dna_client.create(config.API_KEY)
vcf = pysam.VariantFile(config.VCF_PATH)

predictor = AlphaGenomePredictor(
    model, 
    dna_client, 
    genome, 
    window_size=config.WINDOW_SIZE
)
results_df = predictor.analyze_vcf(vcf, desc="Analysing")

#agg
summary_df = results_df.groupby('gene').agg({
    'symbol': 'first',
    'ref_val': 'sum',         
    'alt_val': 'sum', 
    'alt_val_max': 'max',
    'consequence': lambda x: '; '.join(filter(None, x.unique())) 
}).reset_index()

#merge
real_data = pd.read_csv(config.RNA_SEQ_PATH)
real_data['gene_id'] = real_data['gene_id'].str.split('.').str[0]

comparison_df = summary_df.merge(real_data, left_on='gene', right_on='gene_id')

plotter = AlphaGenomePlotter(comparison_df, window_size=WINDOW)
plotter.plot_correlation(pred_col='alt_val', title_suffix="Mean Signal")
plotter.plot_correlation(pred_col='alt_val_max', title_suffix="Max (Peak) Signal")
