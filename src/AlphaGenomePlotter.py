class AlphaGenomePlotter:
    def __init__(self, comparison_df, window_size):
        self.df = comparison_df
        self.window = window_size #just for plot name

    def _prepare_data(self, pred_col):
        """log10 transformation"""
        plot_df = self.df[(self.df[pred_col] > 0) & (self.df['tpm_unstranded'] > 0)].copy()
        plot_df['log_pred'] = np.log10(plot_df[pred_col])
        plot_df['log_tpm'] = np.log10(plot_df['tpm_unstranded'])
        return plot_df

    def plot_correlation(self, pred_col='alt_val', title_suffix="Mean Signal"):
        """
        Plot the correlation between a prediction column and real TPM.
        """
        plot_df = self._prepare_data(pred_col)
        

        #plot
        plt.figure(figsize=(10, 8))
 
        sns.regplot(
            data=plot_df, 
            x='log_pred', 
            y='log_tpm', 
            scatter_kws={'alpha': 0.4, 'edgecolor': 'w'},
            line_kws={'color': 'red', 'lw': 2})

        # Annotate top 10 genes from observed data
        top_genes = plot_df.nlargest(10, 'tpm_unstranded')
        for i, row in top_genes.iterrows():
            plt.text(
                x=row['log_pred'], 
                y=row['log_tpm'], 
                s=row['gene_name'], 
                fontsize=9, fontweight='bold', ha='right', va='bottom',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=0.5)
            )

        #top 10 genes from predicted data
        top_genes = plot_df.nlargest(10, pred_col)
        for i, row in top_genes.iterrows():
            plt.text(
                x=row['log_pred'], 
                y=row['log_tpm'], 
                s=row['gene_name'], 
                color='red',      
                fontsize=9, 
                fontweight='bold', 
                ha='right', 
                va='bottom',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=0.5)
    )

        #corr
        corr_val, p_val = spearmanr(plot_df[pred_col], plot_df['tpm_unstranded'])
        plt.gca().text(
            0.05, 0.95, f'Spearman ρ: {corr_val:.3f}\np-value: {p_val:.2e}',
            transform=plt.gca().transAxes,
            fontsize=12, verticalalignment='top', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )

        plt.xlabel(f'AlphaGenome {pred_col} (log10)')
        plt.ylabel('Real RNA-Seq (log10)')
        plt.title(f'Predictions vs. Real RNAseq ({title_suffix})\nWindow: {self.window}')
        plt.grid(True, which="both", ls="-", alpha=0.2)

        plt.show()
