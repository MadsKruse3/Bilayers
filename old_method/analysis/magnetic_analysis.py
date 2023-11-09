from abc import ABC, abstractmethod
from asr.core import read_json, ASRResult
from typing import List
from pathlib import Path
import numpy as np
from ase.io import read


def get_monolayer_data(monolayer_folder: str, fname: str, file_must_exist: bool) -> ASRResult:
    path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
    c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
    c2db_path2 = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname

    if Path(c2db_path).is_file():
        return read_json(c2db_path)
    elif Path(c2db_path2).is_file():
        return read_json(c2db_path2)
    elif file_must_exist:
        raise ValueError(f"{fname} did not exist for {monolayer_folder}")
    else:
        return None

                        
def get_bilayer_data(bilayer_folder: str, fname: str, file_must_exist: bool) -> ASRResult:
    path = bilayer_folder + "/" + fname

    
    if Path(path).is_file():
        return read_json(path)
    elif file_must_exist:
        raise ValueError(f"{fname} did not exist for {monolayer_folder}")
    else:
        return None


class Analysis(ABC):
    def __init__(self, savelocation):
        self.savelocation = savelocation

    @abstractmethod
    def run_monolayer(self, monolayer_folder):
        pass

    @abstractmethod
    def run_bilayers(self, bilayer_folders):
        pass

    @abstractmethod
    def run_switchables(self, switchables):
        pass

    def run_invertables(self, invertables):
        pass

    @abstractmethod
    def finalize_monolayer(self):
        pass

    @abstractmethod
    def save(self):
        pass

    @abstractmethod
    def plot(self, save):
        pass


class MagstateEmergence(Analysis):
    """Check if BL becomes magnetic or changes magnetic state."""
    def __init__(self, savelocation):
        super().__init__(savelocation)

        self.magstate_mlbl = []
        self.magstate_ml = None

    def run_monolayer(self, monolayer_folder):
        ml_data = get_monolayer_data(monolayer_folder, "results-asr.magstate.json", False)
        if ml_data is None:
            return

        self.magstate_ml = ml_data["magstate"]

    def run_bilayers(self, bilayer_folders):
        for bilayer_folder in bilayer_folders:
            bl_data = get_bilayer_data(bilayer_folder, "results-asr.magstate.json", False)
            if bl_data is None:
                continue

            magstate_bl = bl_data["magstate"]

            self.magstate_mlbl.append((self.magstate_ml, magstate_bl))

    def run_switchables(self, switchables):
        pass

    def finalize_monolayer(self):
        pass

    def save(self):
        data_dct = dict(magstate_mlbl=self.magstate_mlbl)
        np.save(self.savelocation / "magstateemergence.npy", data_dct)

    def get_changes(self, d):
        states = ["NM", "FM", "AFM"]
        changes = {}

        for state1, state2 in d["magstate_mlbl"]:
            if state1 == "NM" and state2 == "NM":
                continue
            if state1 == state2:
                key = "No change"
            else:
                key = state1 + " to " + state2
            if key not in changes:
                changes[key] = 0
            changes[key] += 1

        return changes

    def plot(self, save):
        import matplotlib.pyplot as plt
        import seaborn as sns

        d = np.load(self.savelocation / "magstateemergence.npy", allow_pickle=True).item()

        changes = self.get_changes(d)

        labels = sorted(list(changes.keys()))[::-1]
        values = np.array([changes[label] for label in labels])

        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

        plot_piechart(labels, values, ax)
        
        if save:
            plt.savefig("figures/magstateemergence.pdf", bbox_inches="tight")
        else:
            plt.show()


def plot_piechart(labels, values, ax):
    wedges, texts, autotexts = ax.pie(values, autopct=lambda val: absolute_value(values, val),
                                      shadow=False)
    for autotext in autotexts:
        autotext.set_color("white")

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2 * np.sign(x), 1.2 * y),
                    horizontalalignment=horizontalalignment, **kw)


def absolute_value(arr, val):
    a  = arr[ np.abs(arr - val/100.*arr.sum()).argmin() ]
    return a


def get_interlayer_magstate(bl_data, bilayer_folder):
    """Get magstate based on interlayer coupling.

    First check if proper FM and AFM states exist.
    If they do, return magstate based on relative energies.
    If only one state exists, return that.
    If none exists, return NM.
    """
    from mads_mag import get_bl_magstate
    return get_bl_magstate(bilayer_folder)
    # afm_exists = abs(bl_data["M_AFM"]) < 0.01

    # initial_magmoms = read_json(f"{bilayer_folder}/../structure.json")[1].get("magmoms")
    # bl_atoms = read(f"{bilayer_folder}/structure.json")
    # TM3d_atoms = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]

    # bl_magmoms = bl_data["magmoms_fm"]
    # fm_exists = True
    # for i, atom in enumerate(bl_atoms):
    #     if atom.symbol not in TM3d_atoms:
    #         continue
    #     if np.sign(bl_magmoms[i]) != np.sign(initial_magmoms[i % len(initial_magmoms)]):
    #         fm_exists = False
    #         break

    # if not afm_exists and not fm_exists:
    #     return "NM"
    # elif not afm_exists and fm_exists:
    #     return "FM"
    # elif afm_exists and not fm_exists:
    #     return "AFM"
    # else:
    #     ediff = bl_data["eDIFF"]
    #     if ediff > 0.0:
    #         return "AFM"
    #     else:
    #         return "FM"


class MagneticEmergence(Analysis):
    """Check if BL becomes magnetic or changes magnetic state."""
    def __init__(self, savelocation):
        super().__init__(savelocation)

        self.magstate_mlbl = []
        self.magstate_ml = None

    def run_monolayer(self, monolayer_folder):
        # ml_data = get_monolayer_data(monolayer_folder, "results-asr.magstate.json", False)
        ml_data = get_bilayer_data(monolayer_folder, "results-asr.monolayer_magnetism.json", False)
        if ml_data is None:
            return
        if np.allclose(ml_data["magmoms_FM"], 0.0):
            self.magstate_ml = "NM"
        elif not np.allclose(ml_data["magmoms_FM"], 0.0) and np.allclose(ml_data["M_FM"], 0.0):
            self.magstate_ml = "AFM"
        elif not np.allclose(ml_data["M_FM"], 0.0):
            self.magstate_ml = "FM"
        else:
            raise ValueError("Magstate not determinable. Magmoms: {ml_data['magmoms_FM']}")

        # self.magstate_ml = ml_data["magstate"]

    def run_bilayers(self, bilayer_folders):
        if self.magstate_ml is None:
            return

        for bilayer_folder in bilayer_folders:
            bl_data = get_bilayer_data(bilayer_folder,
                                       "results-asr.interlayer_magnetic_exchange.json", False)
            if bl_data is None:
                continue

            interlayer_magstate = get_interlayer_magstate(bl_data, bilayer_folder)
            if interlayer_magstate is None:
                continue

            self.magstate_mlbl.append((self.magstate_ml, interlayer_magstate, bilayer_folder))

    def run_switchables(self, switchables):
        pass

    def finalize_monolayer(self):
        pass

    def save(self):
        data_dct = dict(magstate_mlbl=self.magstate_mlbl)
        np.save(self.savelocation / "magneticemergence.npy", data_dct)

    def get_changes(self, d):
        changes = {}
        mls = set()

        for state1, state2, bilayer_folder in d["magstate_mlbl"]:
            if state2 == "NM":
                print("This became NM:", bilayer_folder)
            if state1 == state2:
                assert state1 != "NM"
                key = "No change"
            else:
                key = state1 + " to " + state2

            if key not in changes:
                changes[key] = 0
            changes[key] += 1

            ml = bilayer_folder.split("/")[-2]
            mls.add(ml)

        return changes, mls

    def plot(self, save):
        import matplotlib.pyplot as plt
        import seaborn as sns

        d = np.load(self.savelocation / "magneticemergence.npy", allow_pickle=True).item()
        changes, mls = self.get_changes(d)
        print(f"{len(mls)} magnetic monolayers analysed")
        print(next(x for x in mls))

        labels = sorted(list(changes.keys()))[::-1]
        values = np.array([changes[label] for label in labels])

        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

        plot_piechart(labels, values, ax)
        
        if save:
            plt.savefig("figures/magneticeemergence.pdf", bbox_inches="tight")
        else:
            plt.show()


class MagneticSwitching(Analysis):
    """Check if BL becomes magnetic or changes magnetic state."""
    def __init__(self, savelocation):
        super().__init__(savelocation)

        self.magstate_blbl = []

    def run_monolayer(self, monolayer_folder):
        pass

    def run_bilayers(self, bilayer_folders):
        pass

    def run_switchables(self, switchables):
        for mat1, mat2 in switchables:
            data1 = get_bilayer_data(mat1, "results-asr.interlayer_magnetic_exchange.json", False)
            data2 = get_bilayer_data(mat2, "results-asr.interlayer_magnetic_exchange.json", False)

            if data1 is None or data2 is None:
                continue
                
            
            ddip = abs(get_dipole_change(mat1, mat2))
            interlayer_magstate1 = get_interlayer_magstate(data1, mat1)
            interlayer_magstate2 = get_interlayer_magstate(data2, mat2)

            self.magstate_blbl.append((interlayer_magstate1, interlayer_magstate2, ddip,
                                       mat1, mat2))

    def finalize_monolayer(self):
        pass

    def save(self):
        data_dct = dict(magstate_blbl=self.magstate_blbl)
        np.save(self.savelocation / "magneticswitch.npy", data_dct)

    def get_changes(self, d, dip):
        changes = []
        change_counts = {}
        for state1, state2, ddip, mat1, mat2 in d["magstate_blbl"]:
            if abs(ddip) <= dip:
                continue
            if state1 == "NM" or state2 == "NM":
                    continue

            if not state1 in ["FM", "AFM"] or not state2 in ["FM", "AFM"]:
                continue

            ml = mat1.split("/")[-2]
            changes.append((ml, state1, state2, ddip, mat1, mat2))


        return changes

    def plot(self, save):
        import matplotlib.pyplot as plt
        import seaborn as sns

        d = np.load(self.savelocation / "magneticswitch.npy", allow_pickle=True).item()
        changes = self.get_changes(d, 1e-13)
        changes2 = self.get_changes(d, -1)

        ml = set([t[0] for t in changes])
        ml2 = set([t[0] for t in changes2])


        if save:
            pass
        else:
            from ase.db import connect
            c2db = connect("/home/niflheim2/cmr/databases/c2db/c2db.db")
            print(len(changes))
            print(len(changes2))
            print(len(ml))
            print(len(ml2))

            for uid in ml:
                row = c2db.get(f"uid={uid}")
                
                if row.ehull < 0.05:
                    dbid = next((x for x in ["doi", "icsd_id", "cod_id"] if hasattr(row, x)), "none")
                    print(uid, dbid)
            print("-----")
            c = 0
            for uid in ml2:
                row = c2db.get(f"uid={uid}")
                
                if row.ehull < 0.05:
                    dbid = next((x for x in ["doi", "icsd_id", "cod_id"] if hasattr(row, x)), "none")
                    print(uid, dbid)
                    c += 1

            print(c)


        # labels = sorted(list(changes.keys()))[::-1]
        # values = np.array([changes[label] for label in labels])

        # fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

        # plot_piechart(labels, values, ax)
        
        # if save:
        #     plt.savefig("figures/magneticswitch.pdf", bbox_inches="tight")
        # else:
        #     plt.show()






analysis_objs = [EmassEmergence(Path("results")), EmassSwitching(Path("results")),
                 MagstateEmergence(Path("results")), MagneticEmergence(Path("results")),
                 MagneticSwitching(Path("results"))]
# analysis_objs = analysis_objs[-2:] # Do something like this if you want to only run a few analyses
# Can later be changed to a function that returns the analyses based on inputs
