###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
import pandas as pd
import os
import unittest
import tempfile
import shutil

from sciutil import Biomart
from sciutil import SciUtil
from sciutil import SciException


class TestSciUtil(unittest.TestCase):

    def setUp(self):
        self.sciutil = SciUtil()
        self.biomart = Biomart()
        # Setup temp dir
        self.local = False
        if not self.local:
            self.tmp_dir = tempfile.mkdtemp(prefix='scidatannotate_tmp_')
        else:
            self.tmp_dir = '../tests/data/tmp/'

    def tearDown(self):
        if not self.local:
            # Delete temp dr
            shutil.rmtree(self.tmp_dir)

    def test_build_gene_info_file(self):
        gene_ids = ['ENSG00000156575', 'ENSG00000116745']
        file_path = self.sciutil.generate_label([self.tmp_dir, "gene_info"], ".tsv")
        self.biomart.build_gene_info_file('hsapiens_gene_ensembl', gene_ids, self.tmp_dir)
        self.assertEqual(os.path.exists(file_path), True)
        self.assertEqual(len(pd.read_csv(file_path)), 53)

    def test_build_gene_annot_dict(self):
        gene_ids = ['ENSG00000156575', 'ENSG00000116745']
        file_path = self.sciutil.generate_label([self.tmp_dir, "gene_info"], ".tsv")

        # Check exception is raised
        with self.assertRaises(SciException):
            self.biomart.build_gene_annot_dict(file_path)

        self.biomart.build_gene_info_file('hsapiens_gene_ensembl', gene_ids, self.tmp_dir)
        # Now we want to build a dictionary using this information (condensing it for each gene)
        gene_dict, ens_to_name = self.biomart.build_gene_annot_dict(file_path)
        self.assertEqual(ens_to_name.get('ENSG00000116745'), 'RPE65')
        self.assertEqual(gene_dict['PRG3']['gc'], 49.61)

    def test_add_gene_metadata_to_df(self):
        gene_ids = ['ENSG00000156575', 'ENSG00000116745']
        file_path = self.sciutil.generate_label([self.tmp_dir, "gene_info"], ".tsv")
        self.biomart.build_gene_info_file('hsapiens_gene_ensembl', gene_ids, self.tmp_dir)
        # Now we want to build a dictionary using this information (condensing it for each gene)
        gene_dict, ens_to_name = self.biomart.build_gene_annot_dict(file_path)
        # Now lets add some info to a dataframe with only one gene id
        df = pd.DataFrame()
        df['id'] = ['ENSG00000156575']
        df['sample_1'] = [2.88]
        df = self.biomart.add_gene_metadata_to_df(df)
        vals = df.values
        self.assertEqual(vals[0][2], 'PRG3')
        self.assertEqual(vals[0][-1], '10394')
        self.assertEqual(len(vals), 1)

    def test_build_add_metadata(self):
        gene_ids = ['ENSG00000156575', 'ENSG00000116745']
        file_path = self.sciutil.generate_label([self.tmp_dir, "gene_info"], ".tsv")
        df = pd.DataFrame()
        df['id'] = ['ENSG00000156575']
        df['sample_1'] = [2.88]
        annot_df = self.biomart.build_add_metadata(df, 'hsapiens_gene_ensembl', 'id', self.tmp_dir, gene_ids)
        self.assertEqual(annot_df.values[0][2], 'PRG3')
        self.assertEqual(annot_df.values[0][1], 2.88)
