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
import argparse
from sciutil import __version__
import sys
from sciutil.util import bcolors


def print_help():
    lines = [f'sciutil {__version__}',
             'usage: [--help] [--version]',
             'sciutil provides utility functions used across the sci* packages.',
             f'\t {bcolors.WARNING} warn_p(msg:str, sep:str): prints warning message in yellow {bcolors.ENDC}',
             f'\t {bcolors.FAIL} err_p(msg:str, sep:str): prints error message in red {bcolors.ENDC}',
             f'\t {bcolors.OKBLUE} dp(msg:str, sep:str): prints display message in blue {bcolors.ENDC}',
             '\t  get_date_str(): returns the current date string as yyyymmdd',
             '\t  save_df_json(df:pd.DataFrame, output_file:str, keep_index_b=False:Boolean): '
             'saves a dataframe (default no index)',
             '\t  save_plt(fig:matplotlib.plt, output_path:str, dpi=100:int): saves a figure as a PNG',
             '\t  generate_label(label_lst:list, postfix="":str, sep="_":str): generates a label with the date string '
             'postfix is the ending e.g. .csv, separator is the separator for elements in the list ',
             '\t  save_svg(fig:matplotlib.plt, label_lst:list, dpi=100:int): saves a figure as a PNG'
             ]
    print('\n'.join(lines))


def main(args=None):
    parser = argparse.ArgumentParser(description='Utility functions for sci* packages')
    parser.add_argument('help', type=str, help='Prints help descriptions')
    parser.add_argument('version', type=str, help='Prints the version of sciutil')

    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'sciutil v{__version__}')
        sys.exit(0)
    else:
        print_help()


if __name__ == "__main__":
    main()