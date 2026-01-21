class AlphaGenomePredictor:
    def __init__(self, model, dna_client, genome, window_size=131072, uberon = 'UBERON:0002048'):
        self.model = model
        self.dna_client = dna_client
        self.genome = genome
        self.half_window = int(window_size / 2)
        self.csq_columns = []
        self.uberon = uberon

    def _parse_csq_header(self, vcf_file):
        """Extracts CSQ string"""
        if 'CSQ' in vcf_file.header.info:
            desc = vcf_file.header.info['CSQ'].description
            if "Format: " in desc:
                format_str = desc.split("Format: ")[1].strip('">')
                self.csq_columns = format_str.split('|')

    def get_prediction(self, chrom, pos, ref, alt):
        start = pos - self.half_window
        end = pos + self.half_window
        interval = self.genome.Interval(chromosome=chrom, start=start, end=end)
        variant = self.genome.Variant(chromosome=chrom, position=pos, 
                                     reference_bases=ref, alternate_bases=alt)
        
        return self.model.predict_variant(
            interval=interval, variant=variant,
            ontology_terms=[self.uberon], 
            requested_outputs=[self.dna_client.OutputType.RNA_SEQ],
        )

    def parse_annotations(self, record):
        """Extracts genes and function"""
        annotations = []
        if 'CSQ' in record.info:
            for transcript in record.info['CSQ']:
                parts = transcript.split('|')
                data = dict(zip(self.csq_columns, parts)) if self.csq_columns else {}
                gene_id = data.get('Gene') or (parts[4] if len(parts) > 4 else None)
                symbol = data.get('SYMBOL') or (parts[3] if len(parts) > 3 else None)
                effect = data.get('Consequence') or (parts[1] if len(parts) > 1 else None)
                impact = data.get('IMPACT') or (parts[2] if len(parts) > 2 else None)

                if gene_id and gene_id.startswith('ENSG'):
                    annotations.append({
                        'gene': gene_id,
                        'symbol': symbol,
                        'consequence': effect,
                        'impact': impact
                    })
        return annotations

    def analyze_vcf(self, vcf_file, desc="Analysing expression"):
        self._parse_csq_header(vcf_file)
        rows_list = []
        
        #fpr process bar
        total_records = sum(1 for _ in vcf_file)
        vcf_file.reset()

        for record in tqdm(vcf_file, total=total_records, desc=desc):
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]
            
            #get variant annotations
            annots = self.parse_annotations(record)
            if not annots: continue 
            
            outputs = self.get_prediction(chrom, pos, ref, alt)

            #signal extraction
            ref_sig, alt_sig = outputs.reference.rna_seq.values, outputs.alternate.rna_seq.values
            ref_mean, alt_mean = np.mean(ref_sig), np.mean(alt_sig)
            ref_val_max, alt_val_max = np.max(ref_sig), np.max(alt_sig)
            
            for a in annots:
                row = a.copy()
                row.update({
                    'chrom': chrom,
                    'pos': pos,
                    'ref_val': ref_mean, 
                    'alt_val': alt_mean,
                    'ref_val_max':ref_val_max,
                    'alt_val_max':alt_val_max,
                    'diff': alt_mean - ref_mean
                })
                rows_list.append(row)

        df = pd.DataFrame(rows_list)
        return df
