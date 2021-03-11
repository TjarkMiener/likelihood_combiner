"""
jfactor.py
==========
Collection of classes to handle the J-Fator.
"""

import numpy as np
from scipy.interpolate import interp1d
import tables
import os

__all__ = [
    'JFactor',
    'Bonnivard',
    'Custom',
    'GeringerSameth'
]

class JFactor:
    """
    Abstract class for handling the J-Factor.
    """

    def __init__(self,
                channel=None,
                sources=None,
                collaborations=None,
                resource=None,
                combination_data=None,
                precision=2,
                logJ=None,
                DlogJ=None,
                jnuisance=True):
        """
        Parameters
        ----------
        channel: str
            name of the channel.
        sources: `list of str`
            list of the sources.
        collaborations: dict
            settings of the collaborations, following the skeleton: {'name of collaboration' : angular cut `float`}
        resource : path
            path to any file, which might be needed.
        combination_data : path
            path to the dataset, which will be used in the combination.
        precision: int
            presicion of the J-Factor and it's uncertainty.
        logJ: dict
            custom dictionary holding the hardcoded log J-Factor to by-pass the build-in J-Factor sets, following the skeleton:
            {'name of the source': {'name of collaboration' : logJ value `float`}}
        DlogJ: dict
            custom dictionary holding the hardcoded log J-Factor uncertainties to by-pass the build-in J-Factor sets, following the skeleton:
            {'name of the source': {'name of collaboration' : DlogJ value `float`}} 
        jnuisance: bool
            boolean to enable the J-Factor uncertainty as nuisance parameter in the analysis.
        """
    
        self.channel = channel
        self.sources = sources
        self.collaborations = collaborations
        self.resource = resource
        self.combination_data = combination_data
        self.precision = precision
        self.logJ = logJ
        self.DlogJ = DlogJ
        self.jnuisance = jnuisance
        self.DlogJ_comb = None
        self.angular_separations = None
        self.logJ_profile = None
        self.DlogJ_profile = None
        self.combination_info = None
        
    def get_channel(self):
        return self.channel

    def get_sources(self):
        return self.sources

    def get_collaborations(self):
        return self.collaborations

    def get_resource(self):
        return self.resource

    def get_combination_data(self):
        return self.combination_data

    def get_precision(self):
        return self.precision

    def get_logJ(self):
        return self.logJ

    def get_DlogJ(self):
        return self.DlogJ
        
    def get_jnuisance(self):
        return self.jnuisance
        
    def get_DlogJ_comb(self):
        return self.DlogJ_comb
        
    def get_angular_separations(self):
        return self.angular_separations

    def get_logJ_profile(self):
        return self.DlogJ

    def get_DlogJ_profile(self):
        return self.DlogJ_comb

    def get_combination_info(self):
        return self.combination_info
        
    def compute_Jnuisance(self, sigmav_range, ts_values, DlogJ, sigmav_min=1e-28, sigmav_max=1e-18):
        if not np.all(ts_values==ts_values[0]):
            maxdev = 6. * DlogJ
            profiling_steps = 10000
            l = np.linspace(-maxdev, maxdev, num=profiling_steps, endpoint=True)
            lLkl = -2.*np.log((np.exp(-np.power(l,2)/(2.*np.power(DlogJ,2)))/(np.sqrt(2.*np.pi)*DlogJ*np.log(10))))
            lin_interpolation = interp1d(sigmav_range, ts_values, kind='linear', fill_value='extrapolate')
            min_ind = np.min(np.where(ts_values == np.min(ts_values)))
            ts_val = []
            for sv in sigmav_range:
                g = sv/np.power(10,l)
                gSpline = lin_interpolation(g)
                gLkl = np.where(((g>sigmav_min) & (g<sigmav_max)), gSpline,  np.nan)
                totLkl = gLkl + lLkl
                ts_val.append(np.nanmin(totLkl))
            ts_val= np.array(ts_val)
            dis = ts_val[min_ind] - ts_values[min_ind]
            ts_values = ts_val - dis
        return ts_values

    def _get_jfactors(self):
        logJ, DlogJ = {}, {}
        for (source, collaboration) in self.combination_info:
            if source not in logJ and source not in DlogJ:
                logJ[source] = {}
                DlogJ[source] = {}
            angular_separation = self.collaborations[collaboration]
            logJ_profile_interpolation = interp1d(self.angular_separations[source], self.logJ_profile[source], kind='linear', fill_value='extrapolate')
            DlogJ_profile_interpolation = interp1d(self.angular_separations[source], self.DlogJ_profile[source], kind='linear', fill_value='extrapolate')
            logJ[source][collaboration] =  np.around(logJ_profile_interpolation(angular_separation), self.precision)
            DlogJ[source][collaboration] =  np.around(DlogJ_profile_interpolation(angular_separation), self.precision)
        return logJ, DlogJ
        
    
    def _compute_DlogJ_for_combination(self):
        # Compute the actual uncertainties, which is introduced in the combination.
        # Therefore, we have to sort the DlogJ[source] in descending order and compute the
        # actual uncertainty using DlogJ_diff = (DlogJ^2 - DlogJ_next^2)^1/2. For the last
        # element DlogJ_diff is set to (the smallest) DlogJ.
        # After this computation, DlogJ_comb[source] holds the J Factor uncertainty
        # for the combined analysis.
        DlogJ_comb = {}
        for source in self.DlogJ:
            DlogJ_comb[source] = {}
            self.DlogJ[source] = dict(sorted(self.DlogJ[source].items(), key=lambda x: x[1], reverse=True))
            # Sorting the logJ dict accordingly.
            logJ = {}
            for collaboration in self.DlogJ[source]:
                logJ[collaboration] = self.logJ[source][collaboration]
            self.logJ[source] = logJ
            prev_collaboration = None
            for collaboration in self.DlogJ[source]:
                if prev_collaboration:
                    DlogJ_diff = 0.0
                    if DlogJ_comb[source][prev_collaboration] != self.DlogJ[source][collaboration]:
                        DlogJ_diff = np.sqrt(np.power(DlogJ_comb[source][prev_collaboration],2) - np.power(self.DlogJ[source][collaboration],2))
                    DlogJ_comb[source][prev_collaboration] = DlogJ_diff
                prev_collaboration = collaboration
                prev_DlogJ = DlogJ_comb[source][prev_collaboration] = self.DlogJ[source][collaboration]
        return DlogJ_comb

    def _construct_combination_info(self):
        combination_info = []
        if self.combination_data is None:
            # Writing all possible combinations.
            for source in self.sources:
                for collaboration in self.collaborations:
                    combination_info.append((source, collaboration))
        else:
            if self.combination_data.endswith(".h5") or self.combination_data.endswith(".hdf5"):
                # Writing only valid combinations from the lklcom hdf5 file.
                # Opening the hdf5 file.
                h5 = tables.open_file(self.combination_data, "r")
                # Looping over the likelihood or ts tables and append valid combinations.
                for h5_groups in h5.walk_groups("/"):
                    for table in h5.list_nodes(h5_groups, classname='Table'):
                        table_info = table._v_pathname.split("/")
                        if (table_info[2], table_info[1]) not in combination_info and table_info[2] in self.sources and table_info[1] in self.collaborations.keys() and table_info[3] == self.channel:
                            combination_info.append((table_info[2], table_info[1]))

            if os.path.isdir(self.combination_data):
                # Writing only valid combinations from the lklcom hdf5 file.
                txt_files = np.array([x for x in os.listdir(self.combination_data) if x.endswith(".txt")])
                for txt_file in txt_files:
                    file_info = txt_file.replace('.txt','').split("_")
                    if (file_info[1], file_info[2]) not in combination_info and file_info[1] in self.sources and file_info[2] in self.collaborations.keys() and file_info[0] == self.channel:
                        combination_info.append((file_info[1], file_info[2]))
        return combination_info

class Bonnivard(JFactor):
    """
    Jfactor class for handling the J-Factor set from Bonnivard et al. (https://arxiv.org/abs/1504.02048).
    """

    def __init__(self,
                channel,
                sources,
                collaborations,
                resource=None,
                combination_data=None,
                precision=2,
                jnuisance=True):
        """
        Parameters
        ----------
        channel: str
            name of the channel.
        sources: `list of str`
            list of the sources.
        collaborations: dict
            settings of the collaborations, following the skeleton: {'name of collaboration' : angular cut `float`}
        resource : path
            path to any file, which might be needed.
        combination_data : path
            path to the dataset, which will be used in the combination.
        precision: int
            presicion of the J-Factor and it's uncertainty.
        jnuisance: bool
            boolean to enable the J-Factor uncertainty as nuisance parameter in the analysis.
        """

        super().__init__(channel=channel, sources=sources, collaborations=collaborations, resource=resource, combination_data=combination_data, precision=precision, jnuisance=jnuisance)
        
        if resource is None:
            resource = os.path.join(os.path.dirname(__file__), "resources/Bonnivard/")
        self.resource = resource
        self.combination_info = super()._construct_combination_info()
        self.angular_separations, self.logJ_profile, self.DlogJ_profile = self._construct_jprofile()
        self.logJ, self.DlogJ = super()._get_jfactors()
        if self.jnuisance:
            self.DlogJ_comb = super()._compute_DlogJ_for_combination()
                    
    def _construct_jprofile(self):
        # Credits: V.Poireu (LAPP)
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for (source, collaboration) in self.combination_info:
            file = source + "_Jalphaint_cls.output"
            angular_separation, J_array, J_low_array = np.genfromtxt(self.resource + file, skip_header=4, usecols = (0,1,3), unpack=True)
            # The conversion factor is 1 GeV2.cm−5 = 2.25 × 10−7 Msun2.kpc−5
            conversion_factor = 2.25e-7
            logJ = np.log10(J_array / conversion_factor)
            DlogJ = logJ - np.log10(J_low_array / conversion_factor)
            angular_separations[source] = angular_separation
            logJ_profiles[source] = logJ
            DlogJ_profiles[source] = DlogJ
        return angular_separations, logJ_profiles, DlogJ_profiles
        
class Custom(JFactor):
    """
    Custom class to hardcode the J-Factors and it's uncertainties.
    """   

    def __init__(self,
                logJ,
                DlogJ,
                channel=None,
                jnuisance=True):
        """
        Parameters
        ----------
        logJ: dict
            custom dictionary holding the hardcoded log J-Factor to by-pass the build-in J-Factor sets, following the skeleton:
            {'name of the source': {'name of collaboration' : logJ value `float`}}
        DlogJ: dict
            custom dictionary holding the hardcoded log J-Factor uncertainties to by-pass the build-in J-Factor sets, following the skeleton:
            {'name of the source': {'name of collaboration' : DlogJ value `float`}}
        channel: str
            name of the channel.
        jnuisance: bool
            boolean to enable the J-Factor uncertainty as nuisance parameter in the analysis.
        """

        super().__init__(logJ=logJ, DlogJ=DlogJ, channel=channel, jnuisance=jnuisance)
        
        # The arguments logJ and DlogJ should be only used to hardcode the log J-factor
        # and it's uncertainty.
        self.sources = logJ.keys()
        self.collaborations = logJ.values()
        self.logJ = logJ
        self.DlogJ = DlogJ
        self.combination_info = super()._construct_combination_info()
        self.jnuisance = jnuisance
        if self.jnuisance:
            self.DlogJ_comb = super()._compute_DlogJ_for_combination()

class GeringerSameth(JFactor):
    """
    Jfactor class for handling the J-Factor set from Geringer-Sameth et al. (https://arxiv.org/abs/1408.0002).
    """

    def __init__(self,
                channel,
                sources,
                collaborations,
                resource=None,
                combination_data=None,
                precision=2,
                jnuisance=True):
        """
        Parameters
        ----------
        channel: str
            name of the channel.
        sources: `list of str`
            list of the sources.
        collaborations: dict
            settings of the collaborations, following the skeleton: {'name of collaboration' : angular cut `float`}
        resource : path
            path to any file, which might be needed.
        combination_data : path
            path to the dataset, which will be used in the combination.
        precision: int
            presicion of the J-Factor and it's uncertainty.
        jnuisance: bool
            boolean to enable the J-Factor uncertainty as nuisance parameter in the analysis.
        """

        super().__init__(channel=channel, sources=sources, collaborations=collaborations, resource=resource, combination_data=combination_data, precision=precision, jnuisance=jnuisance)
        
        if resource is None:
            resource = os.path.join(os.path.dirname(__file__), "resources/GeringerSameth/intJ_cf.txt")
        self.resource = resource
        self.combination_info = super()._construct_combination_info()
        self.angular_separations, self.logJ_profile, self.DlogJ_profile = self._construct_jprofile()
        self.logJ, self.DlogJ = super()._get_jfactors()
        if self.jnuisance:
            self.DlogJ_comb = super()._compute_DlogJ_for_combination()

    def _construct_jprofile(self):
        angular_separations, logJ_profiles, DlogJ_profiles = {}, {}, {}
        for (source, collaboration) in self.combination_info:
            # Opening txt file.
            file = open(self.resource, "r")
            angular_separation = []
            logJ = []
            DlogJ = []
            for i,line in enumerate(file):
                if i > 41:
                    # The source in the GS file for this line
                    source_GS = line[:18].replace(" ", "")
                    line = line[18:].split()
                    if source == source_GS:
                          angular_separation.append(np.float32(line[0]))
                          logJ.append(np.float32(line[3]))
                          DlogJ.append(np.float32(line[3]) - np.float32(line[2]))
            file.close()
            angular_separations[source] = angular_separation
            logJ_profiles[source] = logJ
            DlogJ_profiles[source] = DlogJ
        return angular_separations, logJ_profiles, DlogJ_profiles
