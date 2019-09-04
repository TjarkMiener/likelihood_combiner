class gloryduckInfo:
    def __init__(self):
        ######################################################
        #              Modify information here:              #
        ######################################################

        # List of the collaborations:
        self.collaborations = ['FermiLAT','HAWC','HESS','MAGIC','VERITAS']

        # List of sources:
        self.sources = ['BootesI', 'CanesVenaticiI', 'CanesVenaticiII', 'Carina', 'ComaBerenices', 'Draco', 'Fornax', 'Hercules', 'LeoI', 'LeoII', 'LeoIV', 'LeoT', 'LeoV', 'Sculptor', 'Segue1', 'Segue2', 'Sextans', 'UrsaMajorI', 'UrsaMajorII', 'UrsaMinor']
        
        # List of annihilation channels:
        self.channels = ['bb', 'tautau', 'mumu', 'WW', 'gammagamma', 'ZZ','ee']

        # Corresponding LaTex notation
        self.channels_LaTex = {'bb':'b\\bar{b}', 'tautau':'\\tau^{+}\\tau^{-}', 'mumu':'\mu^{+}\mu^{-}', 'WW':'W^{+}W^{-}', 'gammagamma':'\gamma\gamma', 'ZZ':'ZZ', 'ee':'e^{+}e^{-}'}
    
        #######################################################
        #                     Thank you!                      #
        #######################################################
