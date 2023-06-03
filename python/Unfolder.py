import ROOT as rt
from collections import namedtuple

Bins = namedtuple('Bins', ['folded', 'unfolded'])


class Unfolder:
    def __init__(self, response_matrix, reg_mode, constraint_area, density_mode,
                 unfolded_bin, folded_bin):

        # 2D unfold as default
        self.bins = Bins(folded_bin, unfolded_bin)
        self.axis1_label = self.bins.unfolded.GetDistributionAxisLabel(0)
        self.axis2_label = self.bins.unfolded.GetDistributionAxisLabel(1)
        self.tau = 0
        # TODO option for 1D or 2D unfold
        # TODO option for TUnfoldIterativeEM
        self.unfold = rt.TUnfoldDensity(response_matrix, rt.TUnfold.kHistMapOutputHoriz,
                                        reg_mode,
                                        constraint_area,
                                        density_mode,
                                        self.bins.unfolded, self.bins.folded)

    def set_background_hist(self, th1_hist_list, name_list):
        pass

    def get_output_hist(self):
        # axis1 name axis2 name
        axis_steering = "dipt[O];dimass[UOC" + str(index) + "]"
