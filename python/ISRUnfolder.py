import Unfolder
import ROOTFiles


class ISRUnfolder:
    def __init__(self):
        self.data_handle = ROOTFiles()
        self.unfolders = {"dipt": None, "dimass": None}
