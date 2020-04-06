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

from datetime import date


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class SciUtil:

    def __init__(self, fig_dir=None, data_dir=None, debug_on=True, plot_on=True, warn_color=bcolors.WARNING,
                 err_color=bcolors.FAIL, msg_color=bcolors.OKBLUE, sep="\t", user_date=None, save_fig=True):
        self.debug_on = debug_on
        self.save_fig = save_fig
        self.plot_on = plot_on
        self.fig_dir = fig_dir
        self.data_dir = data_dir
        self.warn_color = warn_color
        self.err_color = err_color
        self.msg_color = msg_color
        self.sep = sep
        self.date = user_date

    @staticmethod
    def print_msg(msg_lst, sep, color):
        msg = ""
        for i in msg_lst:
            msg += str(i) + sep

        print(color + "-" * 80 + bcolors.ENDC)
        print(color + msg.center(80, ' ') + bcolors.ENDC)
        print(color + '-' * 80 + bcolors.ENDC)

    def warn_p(self, msg_lst, sep=None, color=None):
        """
        Prints an error message. ToDo: Extend this to print to a log file as well.

        Parameters
        ----------
        self
        msg_lst
        sep
        color

        Returns
        -------

        """
        if color is None:
            color = self.warn_color
        if sep is None:
            sep = self.sep
        self.print_msg(msg_lst, sep, color)

    def err_p(self, msg_lst, sep=None, color=None):
        """
        Prints an error message. ToDo: Extend this to print to a log file as well.

        Parameters
        ----------
        self
        msg_lst
        sep
        color

        Returns
        -------

        """
        if color is None:
            color = self.err_color
        if sep is None:
            sep = self.sep

        self.print_msg(msg_lst, sep, color)

    def dp(self, msg_lst, sep=None, color=None):
        """
        Prints the message in a common debug format.
        Has a flag to stop printing too.
        Parameters
        ----------
        sep
        msg_lst
        color

        Returns
        -------

        """
        if self.debug_on:
            if sep is None:
                sep = self.sep
            if color is None:
                color = self.msg_color

            self.print_msg(msg_lst, sep, color)

    def get_date_str(self):
        """
        Helper funtion that returns the date - used typically when saving files.
        Returns
        -------

        """
        if not self.date:
            self.date = date.today().strftime(("%Y%m%d"))
        return self.date

    @staticmethod
    def save_df_json(data_df, outfile_path_str):
        """

        Parameters
        ----------
        data_df
        outfile_path_str

        Returns
        -------

        """

        data_df.to_json(outfile_path_str, orient='index')

    @staticmethod
    def save_df(data_df, outfile_path_str, keep_index_b=False):
        """
        By default don't keep the index, it is just annoying.

        Parameters
        ----------
        data_df
        outfile_path_str
        keep_index_b

        Returns
        -------

        """

        data_df.to_csv(outfile_path_str, index=keep_index_b)

    @staticmethod
    def save_plt(fig, name, dpi=100):
        """
        Save a figure.
        Parameters
        ----------
        fig
        name
        dpi

        Returns
        -------

        """
        fig.savefig(name, format = "png", bbox_inches = 'tight', pad_inches = 0, dpi=dpi)

    def generate_label(self, label_lst, postfix='', sep='_'):
        """

        Parameters
        ----------
        label_lst
        postfix
        sep

        Returns
        -------

        """
        date_time = self.get_date_str()
        label = ''
        for l in label_lst:
            if isinstance(l, str) and l[-1] != '/':
                label += str(l) + sep
            else:
                label += str(l)
        label += sep + date_time
        return label + postfix

    def save_svg(self, plt, label_lst):
        """
        Saves a figure as SVG.

        Parameters
        ----------
        plt
        label_lst

        Returns
        -------

        """
        if self.save_fig:
            label = self.generate_label(label_lst)
            self.dp(["Saving plot", label])
            plt.savefig(label + ".svg")