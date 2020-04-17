from pybiomart import Server
import os
import pandas as pd
from typing import Tuple

from sciutil import SciUtil, SciException


class Biomart:

    def __init__(self, sciutil=None):
        self.u = SciUtil() if not sciutil else sciutil
        self.gene_annot_dict = None
        self.ens_to_name = None

    def build_add_metadata(self, df: pd.DataFrame, organism_lbl: str, ensembl_id_lbl: str, output_dir: str, gene_ids=None):
        """
        Wrapper around the other methods so you can make the calls in one line.
        Builds the metadata info and saves the information to a file before generating the data dictionary.
        Parameters
        ----------
        df
        organism_lbl
        ensembl_id_lbl
        output_dir
        gene_ids

        Returns
        -------

        """
        if gene_ids is None:
            gene_ids = []
            for g in df[ensembl_id_lbl].values:
                gene_ids.append(g)

        gene_info_file, meta_df = self.build_gene_info_file(organism_lbl, gene_ids, output_dir)
        self.build_gene_annot_dict(gene_info_file)
        return self.add_gene_metadata_to_df(df)

    def build_gene_info_file(self, organism_lbl: str, gene_ids: None, output_dir: '.', print_info=False) -> Tuple[str, pd.DataFrame]:
        """
        https://jrderuiter.github.io/pybiomart/

        organism_lbl must be one of the datasets available in ensembl i.e. hsapiens_gene_ensembl,
        ToDo: allow the user to have choice in their query
        Parameters
        ----------
        organism_lbl

        Returns
        -------

        """
        cut_gene_ids = []
        # we ned to remove the version if it exists
        if gene_ids is not None and len(gene_ids) > 0:
            for g in gene_ids:
                cut_gene_ids.append(g.split('.')[0])
        server = Server(host='http://www.ensembl.org')

        dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets[organism_lbl])
        # Prints the dataset information available
        if print_info:
            print(dataset.list_filters())
            print('\n'.join(list(dataset.list_attributes().values[:, 0])))

        if len(cut_gene_ids) > 0:
            gene_info_df = dataset.query(attributes=["ensembl_gene_id", "external_gene_name",  "percentage_gene_gc_content",
                                      "chromosome_name", "start_position","end_position", "strand", "go_id",
                                      "entrezgene_id"],
                          filters={'link_ensembl_gene_id': cut_gene_ids})
        else:
            # This will run for longer as we're getting all genes
            gene_info_df = dataset.query(
                attributes=["ensembl_gene_id", "external_gene_name", "percentage_gene_gc_content",
                            "chromosome_name", "start_position", "end_position", "strand", "go_id",
                            "entrezgene_id"])
        output_path = self.u.generate_label([output_dir, "gene_info"], ".tsv")
        # Save this to a tsv file
        self.u.save_df(gene_info_df, output_path, sep='\t')
        return output_path, gene_info_df

    def build_gene_annot_dict(self, gene_info_file: str) -> Tuple[dict, dict]:
        if not os.path.exists(gene_info_file):
            msg = self.u.msg.msg_file_not_found("build_gene_annot_dict", gene_info_file)
            self.u.err_p([msg])
            raise SciException(msg)
        # Load csv just created, we want a dictionary on gene id, with values for start, end and go terms (as a list)
        gene_dict = {}
        ens_to_name = {}
        ensembl_id, id_idx, gc_idx, chr_idx, start_idx, end_idx, strand_idx, go_idx, ncbi_idx = 0, 1, 2, 3, 4, 5, 6, 7, 8
        end_err, chr_err, start_err, strand_err = 0, 0, 0, 0

        with open(gene_info_file, 'r+') as fp:
            hdr = True
            for line in fp:
                line = line.split('\t')
                try:
                    if not hdr:
                        g_id = line[id_idx].strip()
                        ens_id = line[ensembl_id].strip()
                        # Here if we have multiple go terms they wil be separated by a comma which we change to a pipe
                        # so we can save it as a tab
                        go_term = line[go_idx].strip().replace(',', '|')
                        chrom = 'chr' + line[chr_idx].strip()
                        # i.e. we're only keeping normal chrs.
                        if chrom:
                            gc_content = float(line[gc_idx].strip())
                            start = int(line[start_idx].strip())
                            end = int(line[end_idx].strip())
                            strand = int(line[strand_idx].strip())
                            ncbi = line[ncbi_idx].strip()
                            # Check if we've already added it
                            gene = gene_dict.get(g_id)
                            if gene:
                                # Check if we have the same chrom
                                if gene['chr'] == chrom:
                                    if gene['start'] == start:
                                        if gene['end'] == end:
                                            if gene['direction'] == strand:
                                                if len(go_term) > 1 and len(go_term.split(':')) > 1 and 'GO' in go_term:
                                                    if go_term not in gene['go_terms']:
                                                        gene['go_terms'].append(
                                                            int(go_term.split(':')[1]))  # go_term.split('/')[-1])
                                            else:
                                                strand_err += 1
                                        else:
                                            end_err += 1
                                    else:
                                        start_err += 1
                                else:
                                    chr_err += 1

                            else:
                                ens_to_name[ens_id] = g_id
                                gene_dict[g_id] = {
                                    'id': g_id,
                                    'chr': chrom,
                                    'gc': gc_content,
                                    'start': start,
                                    'end': end,
                                    'direction': strand,
                                    'go_terms': [],
                                    'ncbi': ncbi # Used for mapping to kegg pathways.
                                }
                                if len(go_term.split(':')) > 1 and 'GO' in go_term:
                                    gene_dict[g_id]['go_terms'].append(
                                        int(go_term.split(':')[1]))
                    hdr = False
                except Exception as e:
                    print(line, e)
        # Print out our summary
        self.u.warn_p(
            ["Sorted our gene information, gene dictionary length: ", len(gene_dict), "\n Error log: \t chr: ", chr_err,
             "\t end: ", end_err, "\t start: ", start_err])

        self.gene_annot_dict = gene_dict
        self.ens_to_name = ens_to_name

        return gene_dict, ens_to_name

    def build_roi(self):
        """
        Here we make regions of interest for our peaks to map to.

        We sort this in the same way that our peak files are sorted, otherwise we would not be going
        through it efficiently.
        """
        chr_dict = {}
        start_err = 0

        for gene_id, values in self.gene_annot_dict.items():
            chrom = values['chr']
            start = values['start']
            # Lets make this have a "fake" start based on the TSS
            if values['direction'] < 0:
                start = values['end']

            # Check if we already have the chrom.
            if chr_dict.get(chrom):
                # Check if we have that start already
                if chr_dict[chrom].get(start):
                    chr_dict[chrom][start].append(gene_id)
                else:
                    chr_dict[chrom][start] = [gene_id]
            else:
                chr_dict[chrom] = {}
                chr_dict[chrom][start] = [gene_id]
        # Again, lets use the ordering of the index keys
        chrs_sorted = list(chr_dict.keys())
        chrs_sorted.sort()
        # Now lets make a list of the gene regions of interest
        gene_rois = []
        for c in chrs_sorted:
            starts_labels = list(chr_dict[c].keys())
            starts_labels.sort()
            for i in starts_labels:
                genes = chr_dict[c][i]
                for g_id in genes:
                    self.gene_annot_dict[g_id]['chr'] = self.gene_annot_dict[g_id]['chr']
                    gene_rois.append(self.gene_annot_dict[g_id])

        return gene_rois

    def add_gene_metadata_to_df(self, df: pd.DataFrame, drop_empty_rows=False) -> pd.DataFrame:
        names = []
        terms = []
        gc_content = []
        ncbi = []
        for g_id in df['id'].values:
            g = g_id.split('.')[0]
            if self.ens_to_name.get(g) is None:
                names.append(None)
                terms.append(None)
                gc_content.append(None)
                ncbi.append(None)
            else:
                gene_info = self.gene_annot_dict[self.ens_to_name.get(g)]
                names.append(self.ens_to_name.get(g))
                terms.append(gene_info['go_terms'])
                gc_content.append(gene_info['gc'])
                ncbi_saved = gene_info['ncbi']
                if ncbi_saved != 'NA':
                    ncbi.append(ncbi_saved)
                else:
                    ncbi.append(None)

        df['gene_id'] = names
        df['gc'] = gc_content
        df['go_terms'] = terms
        df['ncbi'] = ncbi
        if drop_empty_rows:
            df = df.dropna()
        return df