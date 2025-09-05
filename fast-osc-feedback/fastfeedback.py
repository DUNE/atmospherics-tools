import uproot
import numpy as np
import ROOT
import awkward as ak
import pandas as pd
from enum import IntEnum
import matplotlib.pyplot as plt
import itertools
import os
import sparse
import polars as pl
from iminuit import Minuit
dir_path = os.path.dirname(os.path.realpath(__file__))

lib_oscprob = dir_path + '/build/oscprob-src/lib/libOscProb.so'
lib_oscillogram = dir_path + '/build/src/libOscillogram.so'

prem_default = dir_path + '/build/oscprob-src/PremTables/prem_15layers.txt'


ROOT.gSystem.Load(lib_oscprob)
ROOT.gSystem.Load(lib_oscillogram)

def replace_empty_list(x, replacement):
    if len(x) == 0:
        return [replacement]
    return x

class Flavor(IntEnum):
    NuE = 12
    NuEBar = -12
    NuMu = 14
    NuMuBar = -14
    NuTau = 16
    NuTauBar = -16
    NC = 0

class MH(IntEnum):
    Normal = 1
    Inverted = -1

class NuAxis(IntEnum):
    NuGen = 0
    NuInt = 1
    NuDet = 2

class VarType(IntEnum):
    VarTrue = 0
    VarReco = 1

class Method(IntEnum):
    Truth = 0
    Reco = 1
    Perso = 2
    Efficiency = 3

channels = [ #[ifl, ofl] ; Aranged in a specific order to match the oscillograms,
    (Flavor.NuE, Flavor.NuE),
    (Flavor.NuE, Flavor.NuMu),
    (Flavor.NuE, Flavor.NuTau),
    (Flavor.NuMu, Flavor.NuE),
    (Flavor.NuMu, Flavor.NuMu),
    (Flavor.NuMu, Flavor.NuTau),
    (Flavor.NuEBar, Flavor.NuEBar),
    (Flavor.NuEBar, Flavor.NuMuBar),
    (Flavor.NuEBar, Flavor.NuTauBar),
    (Flavor.NuMuBar, Flavor.NuEBar),
    (Flavor.NuMuBar, Flavor.NuMuBar),
    (Flavor.NuMuBar, Flavor.NuTauBar),
    (Flavor.NuE, Flavor.NC),
    (Flavor.NuMu, Flavor.NC),
    (Flavor.NuEBar, Flavor.NC),
    (Flavor.NuMuBar, Flavor.NC)
]

def get_nufit(mh = MH.Normal):
    """
    Returns the OscPars object with the specified parameters for neutrino oscillation.

    Parameters:
    - mh (MH): The mass hierarchy of the neutrinos. Default is MH.Normal.

    Returns:
    - pars (OscPars): The OscPars object with the specified parameters.
    """
    pars = ROOT.OscPars()
    pars.dm21 = 7.5e-5
    pars.dm31 = 2.457e-3 if mh == MH.Normal else -2.449e-3 + pars.dm21
    pars.th12 = np.arcsin(np.sqrt(0.304))
    pars.th13 = np.arcsin(np.sqrt(0.0218 if mh == MH.Normal else 0.0219))
    pars.th23 = np.arcsin(np.sqrt(0.452 if mh == MH.Normal else 0.579))
    pars.dcp  = (306 if mh == MH.Normal else 254)*np.pi/180
    return pars

def fix_empty_arrays(data):
    """
    Replaces empty lists in the input data with a specified value.

    Args:
        data (list): The input data containing arrays.

    Returns:
        ak.Array: An array with empty lists replaced by the specified value.
    """
    return ak.where(ak.num(data) == 0, ak.Array([[-999]] * len(data)), data)


class DataManager:
    """
    A class for managing and processing CAFs data.

    Args:
        fname (str): The file name or path of the data file.

    Attributes:
        data (pl.DataFrame): A DataFrame containing the loaded data.

    Methods:
        _load_CAFs_data(fname): Load CAFs data from a file.
        set_flavor_discrimination(method, arg=None): Sets the flavor discrimination method.
        set_energy_reco(method, arg=None): Sets the energy reconstruction method.
        set_direc_reco(method, arg=None): Sets the direction reconstruction method.
        set_nunubar_discrimination(method, arg=None): Sets the nunubar discrimination method.
        prepare_data(): Prepares the data by applying various data transformations.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If an invalid method or argument is provided.

    """

    def default_data_selection(self):
        return (pl.col('cvn_numu') > 0) & (pl.col('cvn_numu') > 0) & (pl.col('recoE_numu') > 0) & (pl.col('recoE_nue') > 0) & (pl.col('direc_numu').abs() <= 1) & (pl.col('direc_nue').abs() <= 1)
    
    def default_direc_reco(self):
        return pl.when(
            (pl.col('npfps') < 3) | (pl.col('recoE_nue') > 1.3) #For high energy events or low number of PFPs, we use the hit direc reco
        ).then(
            pl.col('direc_nc')
        ).otherwise(
            pl.when(
                pl.col('reco_pdg') == Flavor.NuMu
            ).then(
                pl.when(pl.col('direc_numu') != 0).then( #Avoiding cases with zero and moving to the hit direc reco
                    pl.col('direc_numu')
                ).otherwise(
                    pl.col('direc_nc')
                )
            ).otherwise(
                pl.when(
                    pl.col('reco_pdg') == Flavor.NuE
                ).then(
                    pl.when(pl.col('direc_nue') != 0).then( #Avoiding cases with zero and moving to the hit direc reco
                        pl.col('direc_nue')
                    ).otherwise(
                        pl.col('direc_nc')
                    )
                ).otherwise(
                    pl.col('direc_nc')
                )
            )
        )


    def __init__(self, fname):
        self._load_CAFs_data(fname)
        self.set_flavor_discrimination(Method.Reco)
        self.set_direc_reco(Method.Reco)
        self.set_energy_reco(Method.Reco)
        self.set_nunubar_discrimination(Method.Reco)
        self.set_data_selection(self.default_data_selection())

    def _load_CAFs_data(self, fname) -> None:
        """
        Load CAFs data from a file.

        Args:
            fname (str): The file name or path of the data file.

        Returns:
            pl.DataFrame: A DataFrame containing the loaded data.

        Raises:
            FileNotFoundError: If the specified file does not exist.

        """
        print("Loading data...")
        with uproot.open(fname) as f:
            # Load data from the file
            weights = f['weights'].arrays(library='pd')
            weights['nuPDG'] = ak.flatten(f['cafTree/rec/mc/mc.nu.pdg'].array())
            weights['Ev'] = ak.flatten(f['cafTree/rec/mc/mc.nu.E'].array())
            weights['isCC'] = ak.flatten(f['cafTree/rec/mc/mc.nu.iscc'].array())
            weights['NuMomY'] = ak.flatten(f['cafTree/rec/mc/mc.nu.momentum.y'].array())
            weights['mode'] = ak.flatten(f['cafTree/rec/mc/mc.nu.mode'].array())

            recoE_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.lep_calo'].array())
            recoE_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.Enu.e_calo'].array())
            weights['recoE_numu'] = ak.flatten(recoE_numu)
            weights['recoE_nue'] = ak.flatten(recoE_nue)

            direc_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.lngtrk.y'].array())
            direc_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.heshw.y'].array())
            direc_nc = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.dir.heshw.y'].array())
            weights['direc_numu'] = -ak.flatten(direc_numu) #Minus sign to have the convention negative=upgoing neutrinos
            weights['direc_nue'] = -ak.flatten(direc_nue) #Minus sign to have the convention negative=upgoing neutrinos
            weights['direc_nc'] = -ak.flatten(direc_nc) #Minus sign to have the convention negative=upgoing neutrinos

            cvn_nue = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.nue'].array())
            cvn_numu = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.numu'].array())
            cvn_nc = fix_empty_arrays(f['cafTree/rec/common/common.ixn.pandora/common.ixn.pandora.nuhyp.cvn.nc'].array())

            weights['cvn_numu'] = ak.flatten(cvn_numu)
            weights['cvn_nue'] = ak.flatten(cvn_nue)
            weights['cvn_nc'] = ak.flatten(cvn_nc)

            weights['npfps'] = ak.flatten(fix_empty_arrays(f['cafTree/rec/fd/fd.hd.pandora/fd.hd.pandora.npfps'].array()))

        weights['direc_true'] = -weights['NuMomY']/weights['Ev'] #Minus sign to have the convention negative=upgoing neutrinos
        weights['nue_w'] *= weights["xsec"]
        weights['numu_w'] *= weights["xsec"]

        self.data = pl.from_pandas(pd.DataFrame(weights))
        print("Finished loading data")

    def set_flavor_discrimination(self, method:Method, arg=None):
        """
        Sets the flavor discrimination method based on the given method.

        Parameters:
            method (Method): The method to use for flavor discrimination.
            arg (pl.Expr, optional): An expression to use for the 'Method.Perso' method.

        Raises:
            ValueError: If an invalid method is provided or if 'arg' is not a polars expression.

        Returns:
            None
        """
        if method == Method.Truth:
            self.flavor_discrimination = lambda: pl.col("nuPDG")
        elif method == Method.Reco:
            self.flavor_discrimination = lambda: pl.when(
                pl.col('cvn_nue') > pl.col('cvn_numu')
            ).then(
                pl.when(
                    pl.col('cvn_nue') > pl.col('cvn_nc')
                ).then(
                    Flavor.NuE
                ).otherwise(
                    Flavor.NC
                )
            ).otherwise(
                pl.when(
                    pl.col('cvn_numu') > pl.col('cvn_nc')
                ).then(
                    Flavor.NuMu
                ).otherwise(
                    Flavor.NC
                )
            )
        elif method == Method.Perso:
            if not isinstance(arg, pl.Expr):
                raise ValueError("A polars expression is expected when using the Method.Perso method")
            self.flavor_discrimination = lambda: arg
        elif method == Method.Efficiency:
            raise ValueError("Method.Efficiency is not implemented for the flavor discrimination")
        else:
            raise ValueError()
        
    def set_energy_reco(self, method:Method, arg=None):
        """
        Sets the energy reconstruction method based on the given method and argument.

        Parameters:
            method (Method): The energy reconstruction method to use.
            arg (optional): An argument specific to the chosen method.

        Raises:
            ValueError: If an invalid method or argument is provided.

        Returns:
            None
        """
        if method == Method.Truth:
            self.energy_reco = lambda: pl.col('Ev')
        elif method == Method.Reco:
            self.energy_reco = lambda: pl.when(
                pl.col('reco_pdg') == Flavor.NuMu
            ).then(
                pl.col('recoE_numu')
            ).otherwise(
                pl.col('recoE_nue')
            )
        elif method == Method.Perso:
            if not isinstance(arg, pl.Expr):
                raise ValueError("A polars expression is expected when using the Method.Perso method")
            self.energy_reco = lambda: arg
        elif method == Method.Efficiency:
            if not isinstance(arg, FakeResolution):
                raise ValueError("A FakeResolution object is expected when using the Method.Efficiency method")
            self.energy_reco = lambda: pl.lit(arg.generate(self.data.select(arg.bin_var).to_series(), self.data['Ev']))
        else:
            raise ValueError()
        
    def set_direc_reco(self, method:Method, arg=None):
        """
        Sets the direction reconstruction method based on the given method and argument.

        Parameters:
            method (Method): The method to use for direction reconstruction.
            arg: The argument to be used with the method.

        Raises:
            ValueError: If an invalid method is provided or if the argument is of the wrong type.

        Returns:
            None
        """

        if method == Method.Truth:
            self.direc_reco = lambda: pl.col('direc_true')
        elif method == Method.Reco:
            self.direc_reco = self.default_direc_reco
        elif method == Method.Perso:
            if not isinstance(arg, pl.Expr):
                raise ValueError("A polars expression is expected when using the Method.Perso method")
            self.direc_reco = lambda: arg
        elif method == Method.Efficiency:
            if not isinstance(arg, FakeResolution):
                raise ValueError("A polars FakeResolution object is expected when using the Method.Efficiency method")
            self.direc_reco = lambda: pl.lit(arg.generate(self.data[arg.bin_var], self.data['direc_true']))
        else:
            raise ValueError()

    def set_nunubar_discrimination(self, method:Method, arg=None):
        """
        Sets the nunubar discrimination method based on the given method and argument.

        Parameters:
            method (Method): The method to use for nunubar discrimination.
            arg (optional): An argument specific to the chosen method.

        Raises:
            ValueError: If an invalid method or argument is provided.

        Returns:
            None
        """
        if method == Method.Truth:
            self.nunubar_discrimination = lambda: pl.col('reco_pdg').abs()*pl.col('nuPDG').sign()
        elif method == Method.Reco:
            self.nunubar_discrimination = lambda: pl.col('reco_pdg')
        elif method == Method.Perso:
            if not isinstance(arg, pl.Expr):
                raise ValueError("A polars expression is expected when using the Method.Perso method")
            self.nunubar_discrimination = lambda: arg
        elif method == Method.Efficiency:
            if not isinstance(arg, FakeEfficiency):
                raise ValueError("A polars FakeEfficiency object is expected when using the Method.Efficiency method")
            self.nunubar_discrimination = lambda: arg.generate(pl.col('nuPDG'), len(self.data)).sign()*pl.col('reco_pdg').abs()
        else:
            raise ValueError()
        
    
    def set_data_selection(self, selection:pl.Expr):
        """
        Sets the data selection based on the given selection expression.

        Parameters:
            selection (pl.Expr): A polars expression representing the selection criteria.

        Raises:
            ValueError: If the selection is not a polars expression.

        Returns:
            None
        """
        if not isinstance(selection, pl.Expr):
            raise ValueError("A polars expression is expected for the selection")
        self.data_selection = selection

    def prepare_data(self):
        """
        Prepares the data by applying various data transformations.

        Returns:
            prepared_data (Table): The prepared data with additional columns.
        """
        prepared_data = self.data.filter(
            self.data_selection
        ).with_columns(
            reco_pdg = self.flavor_discrimination()
        ).with_columns(
            reco_pdg = self.nunubar_discrimination()
        ).with_columns(
            recoE = self.energy_reco()
        ).with_columns(
            direc_reco = self.direc_reco()
        )

        return prepared_data
    

class FakeEfficiency:
    """
    A class representing fake efficiencies for different flavors.

    Parameters:
    - efficiencies (dict[Flavor, float]): A dictionary containing the efficiencies for each flavor.

    Raises:
    - ValueError: If the efficiencies dictionary does not contain entries for all required flavors.

    Methods:
    - generate(true_column: pl.Expr, size: int) -> pl.Expr: Generates fake data based on the true column and size.
    """

    def __init__(self, efficiencies: dict[Flavor, float]):
        """
        Initializes a FakeEfficiency object.

        Parameters:
        - efficiencies (dict[Flavor, float]): A dictionary containing the efficiencies for each flavor.

        Raises:
        - ValueError: If the efficiencies dictionary does not contain entries for all required flavors.
        """
        required_flavors = [Flavor.NuE, Flavor.NuMu, Flavor.NuEBar, Flavor.NuMuBar]
        if not all(flavor in efficiencies for flavor in required_flavors):
            raise ValueError("The efficiencies dictionary must contain entries for all required flavors")

        self._efficiencies = efficiencies

    def generate(self, true_column: pl.Expr, size: int) -> pl.Expr:
        """
        Generates fake data based on the true column and size.

        Parameters:
        - true_column (pl.Expr): The true column to generate fake data from.
        - size (int): The size of the generated data.

        Returns:
        - pl.Expr: The generated fake data.
        """
        random_vals = np.random.rand(size)

        return pl.when(
            pl.lit(random_vals) < true_column.replace_strict(old=pl.Series(self._efficiencies.keys()), new=pl.Series(self._efficiencies.values()))
        ).then(
            true_column
        ).otherwise(
            -true_column
        )


class FakeResolution:
    """
    A class that generates fake values based on given resolutions and binning.

    Parameters:
    - bins (list): A list of bin values.
    - resolutions (list): A list of resolutions corresponding to each bin.
    - bin_var (pl.Expr): The variable used for binning.

    Raises:
    - ValueError: If the length of bins plus 1 is not equal to the length of resolutions.

    Methods:
    - generate(bin_var_values, true_values): Generates fake values based on the given bin variable values and true values.

    """

    def __init__(self, bins, resolutions, bin_var:pl.Expr):
        if len(bins) + 1 != len(resolutions):
            raise ValueError("The resolutions must include the values outside the binning (below lowest bin and above highest bin)")
            
        self._bins = bins
        self._resolutions = resolutions
        self.bin_var = bin_var

    def generate(self, bin_var_values:pl.Series, true_values:pl.Series):
        """
        Generates fake values based on the given bin variable values and true values.

        Parameters:
        - bin_var_values (pl.Series): The bin variable values.
        - true_values (pl.Series): The true values.

        Returns:
        - fake_values (np.ndarray): An array of generated fake values.

        """
        bin_sort = np.digitize(bin_var_values, self._bins)
        fake_values = true_values*np.random.normal(1, self._resolutions[bin_sort])
        return fake_values

class EventDistrib:
    """
    Class for computing event distributions and oscillograms.

    Args:
        events (polars.DataFrame): DataFrame containing event data.
        Ebins (array-like): Energy bin edges for true energy.
        Czbins (array-like): Cosine zenith angle bin edges for true direction.
        Ebins_reco (array-like): Energy bin edges for reconstructed energy.
        Czbins_reco (array-like): Cosine zenith angle bin edges for reconstructed direction.
        fieldIsCC (str, optional): Name of the column indicating whether the event is charged current. Default is 'isCC'.
        fieldNuGen (str, optional): Name of the column indicating the neutrino PDG code. Default is 'nuPDG'.
        fieldNueFlux (str, optional): Name of the column indicating the nue flux weight. Default is 'nue_w'.
        fieldNumuFlux (str, optional): Name of the column indicating the numu flux weight. Default is 'numu_w'.
        fieldTrueE (str, optional): Name of the column indicating the true energy. Default is 'Ev'.
        fieldTrueDir (str, optional): Name of the column indicating the true direction. Default is 'direc_true'.
        fieldRecoE (str, optional): Name of the column indicating the reconstructed energy. Default is 'recoE'.
        fieldRecoDir (str, optional): Name of the column indicating the reconstructed direction. Default is 'direc_reco'.
        fieldRecoFlv (str, optional): Name of the column indicating the reconstructed flavor. Default is 'reco_pdg'.

    Attributes:
        events (polars.DataFrame): DataFrame containing event data.
        fieldIsCC (str): Name of the column indicating whether the event is charged current.
        fieldNuGen (str): Name of the column indicating the neutrino PDG code.
        fieldNueFlux (str): Name of the column indicating the nue flux weight.
        fieldNumuFlux (str): Name of the column indicating the numu flux weight.
        fieldTrueE (str): Name of the column indicating the true energy.
        fieldTrueDir (str): Name of the column indicating the true direction.
        fieldRecoE (str): Name of the column indicating the reconstructed energy.
        fieldRecoDir (str): Name of the column indicating the reconstructed direction.
        fieldRecoFlv (str): Name of the column indicating the reconstructed flavor.
        hists (dict): Dictionary containing the histograms for different channels.
        Ebins (array-like): Energy bin edges for true energy.
        Czbins (array-like): Cosine zenith angle bin edges for true direction.
        Ebins_reco (array-like): Energy bin edges for reconstructed energy.
        Czbins_reco (array-like): Cosine zenith angle bin edges for reconstructed direction.
        osc (ROOT.Oscillogram): Oscillogram object for computing oscillograms.
        oscillograms (np.ndarray): Array containing the computed oscillograms.

    Methods:
        _fill_hists(Ebins, Czbins, Ebins_reco, Czbins_reco): Fills histograms with event data.
        compute_osc(pars): Compute the oscillograms based on the given parameters.
        true_distrib_ch(ch, oscillated=False): Calculate the true distribution for a given channel.
        reco_distrib_ch(ch, oscillated=False): Calculate the reconstructed distribution for a given channel.
        get_distrib(vartype, axis, fl, oscillated=False): Returns the distribution based on the given parameters.
        detected_events(detected_fl): Returns the distribution of detected events for a given flavor.
        get_oscillogram(ifl, ofl): Get the oscillogram for the given input and output channels.
        plot_distrib(vartype, axis, fl, oscillated=False): Plot the distribution of a variable based on the given parameters.
    """

    def __init__(self, events, Ebins, Czbins, Ebins_reco, Czbins_reco,
                 fieldIsCC='isCC', fieldNuGen="nuPDG", fieldNueFlux="nue_w",
                 fieldNumuFlux="numu_w", fieldTrueE="Ev", fieldTrueDir="direc_true",
                 fieldRecoE="recoE", fieldRecoDir="direc_reco", fieldRecoFlv="reco_pdg",
                 prem=prem_default):
        
        self.events = events
        self.fieldIsCC = fieldIsCC
        self.fieldNuGen = fieldNuGen
        self.fieldNueFlux = fieldNueFlux
        self.fieldNumuFlux = fieldNumuFlux
        self.fieldTrueE = fieldTrueE
        self.fieldTrueDir = fieldTrueDir
        self.fieldRecoE = fieldRecoE
        self.fieldRecoDir = fieldRecoDir
        self.fieldRecoFlv = fieldRecoFlv

        self._fill_hists(Ebins, Czbins, Ebins_reco, Czbins_reco)
        self.osc = ROOT.Oscillogram(Ebins, Czbins, prem)

    def _fill_hists(self, Ebins, Czbins, Ebins_reco, Czbins_reco):
        """
        Fills histograms with event data.

        Args:
            Ebins (array-like): Energy bin edges for true energy.
            Czbins (array-like): Cosine zenith angle bin edges for true direction.
            Ebins_reco (array-like): Energy bin edges for reconstructed energy.
            Czbins_reco (array-like): Cosine zenith angle bin edges for reconstructed direction.
        """
        print("Filling histograms...")
        self.hists = {}
        self.Ebins = Ebins
        self.Czbins = Czbins
        self.Ebins_reco = Ebins_reco
        self.Czbins_reco = Czbins_reco
        
        events_cc = self.events.filter(
            pl.col(self.fieldIsCC) == 1
        )

        events_nc = self.events.filter(
            pl.col(self.fieldIsCC) == 0
        )

        detected_channels = self.events[self.fieldRecoFlv].unique()
        self.detected_channels = [Flavor(ch) for ch in detected_channels]
        print("Using the following detected output channels:", [ch.name for ch in self.detected_channels])

        channels_with_detected = itertools.product(channels, self.detected_channels)

        for ch in channels_with_detected:
            (ifl, ofl), detected_fl = ch

            if ofl == Flavor.NC:
                selected = events_nc.filter(
                    pl.col(self.fieldRecoFlv) == detected_fl.value
                )
            else:
                selected = events_cc.filter(
                    pl.col(self.fieldNuGen) == ofl.value,
                    pl.col(self.fieldRecoFlv) == detected_fl.value
                )
            
            
            if (ifl == Flavor.NuE) or (ifl == Flavor.NuEBar):
                weights = selected[self.fieldNueFlux]
            else:
                weights = selected[self.fieldNumuFlux]

            full_ch = (ifl, ofl, detected_fl)

            hist, _ = np.histogramdd(
                (selected[self.fieldTrueDir],
                 selected[self.fieldTrueE],
                 selected[self.fieldRecoDir],
                 selected[self.fieldRecoE]
                ),
                bins=(Czbins, Ebins, Czbins_reco, Ebins_reco),
                weights=weights
                )
            
            self.hists[full_ch] = sparse.COO(hist)
            # self.hists[ch], _ = np.histogramdd((selected['direc_true'], selected['Ev'], selected['direc_true'], selected['Ev']), bins=(Czbins, Ebins, Czbins_reco, Ebins_reco), weights=weights)
        print("Finished filling histograms...")

    def compute_osc(self, pars):
        """
        Compute the oscillograms based on the given parameters.

        Parameters:
        - pars: The parameters used for computing the oscillograms.

        Returns:
        - None

        This method computes the oscillograms using the provided parameters and stores the result in the `oscillograms` attribute of the object.
        The computed oscillograms are reshaped into a 3-dimensional array based on the dimensions of `Czbins` and `Ebins`.
        """
        oscillograms = self.osc.Compute(pars)
        self.oscillograms = np.array(oscillograms).reshape((len(self.Czbins) - 1, len(self.Ebins) - 1, 12))

    def true_distrib_ch(self, ch, oscillated=False):
        """
        Calculate the true distribution for a given channel.

        Parameters:
        - ch (tuple): A tuple representing the channel, where the first element is the channel number and the second element is the flavor.
        - oscillated (bool): A flag indicating whether the distribution should be oscillated or not. Default is False.

        Returns:
        - numpy.ndarray: The calculated true distribution for the given channel.
        """
        if oscillated and ch[1] != Flavor.NC:
            idx = channels.index((ch[0], ch[1]))
            osc = self.oscillograms[:, :, idx, None, None]
        else:
            osc = 1.
        return np.sum(self.hists[ch]*osc, axis=(2, 3))
    
    def reco_distrib_ch(self, ch, oscillated=False):
        """
        Calculate the reconstructed distribution for a given channel.

        Parameters:
        - ch (tuple): A tuple representing the channel (flavor, interaction type).
        - oscillated (bool): Flag indicating whether the channel is oscillated.

        Returns:
        - np.ndarray: The reconstructed distribution for the given channel.
        """
        if oscillated and ch[1] != Flavor.NC:
            idx = channels.index((ch[0], ch[1]))
            osc = self.oscillograms[:, :, idx, None, None]
        else:
            osc = 1
        return np.sum(self.hists[ch]*osc, axis=(0, 1))
    
    def get_distrib(self, vartype:VarType, axis:NuAxis, fl:Flavor, oscillated:bool=False):
        """
        Returns the distribution based on the given parameters.

        Args:
            vartype (VarType): The type of variable (VarTrue or VarReco).
            axis (NuAxis): The axis to consider (NuGen, NuInt, or NuDet).
            fl (Flavor): The flavor to match.
            oscillated (bool, optional): Whether to consider oscillated distribution. Defaults to False.

        Returns:
            final_distrib: The final distribution based on the given parameters.
        """
        final_distrib = None
        if axis == NuAxis.NuGen:
            match = lambda t: t[0] == fl
        elif axis == NuAxis.NuInt:
            match = lambda t: t[1] == fl
        elif axis == NuAxis.NuDet:
            match = lambda t: t[2] == fl
        else:
            raise ValueError()
        
        if vartype == VarType.VarTrue:
            distrib_func = lambda ch: self.true_distrib_ch(ch, oscillated)
        elif vartype == VarType.VarReco:
            distrib_func = lambda ch: self.reco_distrib_ch(ch, oscillated)
        else:
            raise ValueError()
        
        for ch, _ in self.hists.items():
            if not match(ch):
                continue
            if final_distrib is None:
                final_distrib = distrib_func(ch)
            else:
                final_distrib += distrib_func(ch)
        return final_distrib
    
    def detected_events(self, detected_fl: Flavor):
        """
        Returns the distribution of detected events for a given flavor.

        Parameters:
        - detected_fl (Flavor): The flavor of the detected events.

        Returns:
        - Distribution: The distribution of detected events.

        """
        return self.get_distrib(VarType.VarReco, NuAxis.NuDet, detected_fl, True)
    
    def get_oscillogram(self, ifl, ofl):
        """
        Get the oscillogram for the given input and output channels.

        Parameters:
        - ifl (InputChannel): The input channel.
        - ofl (OutputChannel): The output channel.

        Returns:
        - fig (matplotlib.figure.Figure): The generated figure object.

        Raises:
        - ValueError: If no oscillogram is found for the given input and output channels.
        """
        if not (ifl, ofl) in channels:
            raise ValueError(f"No oscillogram found for {ifl.name} -> {ofl.name}")
        osc_id = channels.index((ifl, ofl))

        fig = plt.figure()
        plt.pcolormesh(self.Ebins, self.Czbins, self.oscillograms[:, :, osc_id], cmap='jet')
        plt.xscale('log')
        plt.title(f"{ifl.name} -> {ofl.name}")
        plt.xlabel("E [GeV]")
        plt.ylabel("Zenith angle")
        plt.tight_layout()
        return fig
    
    def plot_distrib(self, vartype:VarType, axis:NuAxis, fl:Flavor, oscillated:bool=False):
        """
        Plot the distribution of a variable based on the given parameters.

        Parameters:
        - vartype (VarType): The type of variable to plot (VarType.VarTrue or VarType.VarReco).
        - axis (NuAxis): The axis to plot (NuAxis.Energy or NuAxis.CosZenith).
        - fl (Flavor): The flavor to plot (Flavor.NuMu, Flavor.NuMuBar, Flavor.NuE, or Flavor.NuEBar).
        - oscillated (bool, optional): Whether to plot the oscillated distribution. Default is False.

        Returns:
        - fig (matplotlib.figure.Figure): The generated figure object.

        """
        distrib = self.get_distrib(vartype, axis, fl, oscillated).todense()
        fig = plt.figure()
        if vartype == VarType.VarTrue:
            Ebins = self.Ebins
            Czbins = self.Czbins
        else:
            Ebins = self.Ebins_reco
            Czbins = self.Czbins_reco

        plt.pcolormesh(Ebins, Czbins, distrib, cmap='jet')
        plt.xscale('log')
        plt.title(f"{vartype.name} ; {axis.name} ; {fl.name}")
        plt.xlabel("E [GeV]")
        plt.ylabel("Zenith angle")
        plt.tight_layout()
        return fig
    
    
def chi2(observed, expected):
    """
    Calculate the chi-square statistic for comparing observed and expected values.

    Parameters:
    observed (array-like): The observed values.
    expected (array-like): The expected values.

    Returns:
    float: The chi-square statistic.

    """
    mask = expected > 0
    return np.sum(((observed[mask] - expected[mask])**2)/expected[mask])

def LnL(data, mc):
    """
    Calculate the log-likelihood (llh) using the given data and model counts.

    Parameters:
    data (numpy.ndarray): The observed data.
    mc (numpy.ndarray): The model counts.

    Returns:
    float: The log-likelihood value.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        llh = np.sum(np.where(mc > 0,
            np.where(data > 0,  
                    2*(mc - data + data*np.log(data/mc)),
                    2*mc
                        ),
            np.where(data > 0, 2*data, 0)
        ))

    return llh

class OscFit:
    def __init__(self, events):
        self.events = events

        self.__setup()
       
    def __setup(self):
        nufit = get_nufit(MH.Normal)

        fitter = Minuit(self.loss, dm21=nufit.dm21, dm31=nufit.dm31, th12=nufit.th12, th13=nufit.th13, th23=nufit.th23, dcp=nufit.dcp)

        self.fitter = fitter

    def set_limits(self, limits):
        """
        Set the limits for the fitting parameters.

        Parameters:
        - limits (dict): A dictionary containing the limits for each parameter.
        """
        self.fitter.limits = limits

    def loss(self, dm21, dm31, th12, th13, th23, dcp):
        pars = ROOT.OscPars()
        pars.dm21 = dm21
        pars.dm31 = dm31
        pars.th12 = th12
        pars.th13 = th13
        pars.th23 = th23
        pars.dcp  = dcp

        self.events.compute_osc(pars)
        new_rates = {fl: self.events.detected_events(fl) for fl in self.events.detected_channels}

        total_loss = 0
        for ch in self.events.detected_channels:
            total_loss += LnL(self.data_rates[ch], new_rates[ch])
        return total_loss

    def fit(self, null_hyp, nuisance_params, asimov_true_params, poisson_throws=False):
        self.events.compute_osc(null_hyp)
        self.data_rates = {fl: self.events.detected_events(fl) for fl in self.events.detected_channels}

        for par in self.fitter.parameters:
            if par in nuisance_params:
                self.fitter.fixed[par] = False
            else:
                self.fitter.fixed[par] = True
            self.fitter.values[par] = getattr(asimov_true_params, par)

        return self.fitter.migrad()

    def plot(self):
        # Plot the results
        pass