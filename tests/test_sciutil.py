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
import unittest
from sciutil import SciUtil
from datetime import date


class TestSciUtil(unittest.TestCase):

    def setUp(self):
        self.sciutil = SciUtil()

    def test_dp(self):
        self.sciutil.dp(["Test printing default message: test", 1])

    def test_warn_p(self):
        self.sciutil.warn_p(["Test printing warning message: test", 2])

    def test_err_p(self):
        self.sciutil.err_p(["Test printing error message: test", 3], sep='|')

    def test_get_date_str(self):
        self.assertEqual(self.sciutil.get_date_str(), date.today().strftime(("%Y%m%d")))

    def test_generate_label(self):
        date_str = date.today().strftime(("%Y%m%d"))
        self.assertEqual(self.sciutil.generate_label(['l1', 2], '.csv'),
                         f'l1_2_{date_str}.csv')

    def test_check_dir_format(self):
        
        dir_str_good = 'adir' + self.sciutil.dir_sep
        self.assertEqual(dir_str_good, self.sciutil.check_dir_format(dir_str_good))

        dir_str_bad = 'adir'
        self.assertEqual(dir_str_good, self.sciutil.check_dir_format(dir_str_bad))