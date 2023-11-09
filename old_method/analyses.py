from abc import ABC, abstractmethod
from asr.core import read_json, ASRResult
from typing import List
from pathlib import Path
import numpy as np
from ase.io import read


def get_dipole_change(folder1: str, folder2: str, blbl=True) -> float:
    if not blbl:
        raise NotImplementedError

    atoms = read(f"{folder1}/structure.json")
    area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))

    dip1 = get_dip(folder1, area)
    dip2 = get_dip(folder2, area)
    
    return dip1 - dip2


def get_dip(folder, area, ml=False):
    """Returns dipole moment in units of C/m.

    Technically dipole moment density."""
    if not ml:
        dpath = f"{folder}/results-asr.gs.json"
        if not Path(dpath).is_file():
            raise ValueError(f"asr.gs does not exist for {folder}")
        data = read_json(dpath)
    else:
        data = get_monolayer_data(folder)
        
    dip = data["dipz"] / area * 1.602 * 1e-19 / 1e-10
    return dip


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


def get_masses(data, bt):
    masses = [data[key]
              for key in [f"emass_{bt}_dir" + str(i)
                          for i in [1,2,3]]
              if data[key] is not None]

    return masses

def get_holes_elecs(data):
    holes = get_masses(data, "vb")
    elecs = get_masses(data, "cb")
    return holes, elecs


class EmassEmergence(Analysis):
    """Calculate changes in effective masses from ML to BL.

    Compares heavy masses and light masses."""
    def __init__(self, savelocation):
        super().__init__(savelocation)
        self.light_hole_mlbl = []
        self.heavy_hole_mlbl = []
        self.light_elec_mlbl = []
        self.heavy_elec_mlbl = []

        self.light_hole_ml = []
        self.heavy_hole_ml = []
        self.light_elec_ml = []
        self.heavy_elec_ml = []

        self.light_hole_bl = []
        self.heavy_hole_bl = []
        self.light_elec_bl = []
        self.heavy_elec_bl = []

        self.light_hole = None
        self.heavy_hole = None
        self.light_elec = None
        self.heavy_elec = None

    def run_monolayer(self, monolayer_folder):
        """Get monolayer masses for later comparison."""
        data = get_monolayer_data(monolayer_folder, "results-asr.emasses@validate.json", False)
        if data is None:
            return
        else:
            try:
                # Some data is partially missing or otherwise corrupt
                holes = get_masses(data, "vb")
                elecs = get_masses(data, "cb")
            except KeyError as e:
                return

            if any(hole > 0.0 for hole in holes) or any(elec < 0.0 for elec in elecs):
                return
            elif any(abs(hole) > 100 for hole in holes) or any(abs(elec) > 100 for elec in elecs):
                return
            else:
                self.light_hole = min(holes, key=lambda t:abs(t))
                self.heavy_hole = max(holes, key=lambda t:abs(t))
                self.light_elec = min(elecs, key=lambda t:abs(t))
                self.heavy_elec = max(elecs, key=lambda t:abs(t))
            
    def run_bilayers(self, bilayer_folders):
        """Get bilayer masses and calculate changes."""
        if self.light_hole is None:
            # The monolayer is not insulating
            # We don't need to include these here because they are
            # already captured in the bandstructure analysis
            # Perhaps it would be interesting to see how large the masses
            # get when the monolayer is metallic but that is for another time
            return
        else:
            for bilayer_folder in bilayer_folders:
                data = get_bilayer_data(bilayer_folder, "results-asr.emasses@validate.json", False)
                if data is None:
                    # Bilayer is not insulating or data is missing
                    continue
                else:
                    try:
                        holes = get_masses(data, "vb")
                        elecs = get_masses(data, "cb")
                    except KeyError as e:
                        continue

                    if (any(abs(hole) > 100 for hole in holes)
                        or any(abs(elec) > 100 for elec in elecs)):
                        continue

                    light_hole = min(holes, key=lambda t:abs(t))
                    heavy_hole = max(holes, key=lambda t:abs(t))
                    light_elec = min(elecs, key=lambda t:abs(t))
                    heavy_elec = max(elecs, key=lambda t:abs(t))

                    self.light_hole_mlbl.append((self.light_hole, light_hole))
                    self.heavy_hole_mlbl.append((self.heavy_hole, heavy_hole))
                    self.light_elec_mlbl.append((self.light_elec, light_elec))
                    self.heavy_elec_mlbl.append((self.heavy_elec, heavy_elec))

                    self.light_hole_bl.append((bilayer_folder, light_hole))
                    self.heavy_hole_bl.append((bilayer_folder, heavy_hole))
                    self.light_elec_bl.append((bilayer_folder, light_elec))
                    self.heavy_elec_bl.append((bilayer_folder, heavy_elec))
            
    def run_switchables(self, switchables):
        pass

    def finalize_monolayer(self):
        """Finalize calculation.
        
        Do any final steps necessary for analysis."""
        pass

    def save(self):
        data_dct = dict(light_hole_mlbl=self.light_hole_mlbl,
                        heavy_hole_mlbl=self.heavy_hole_mlbl,
                        light_elec_mlbl=self.light_elec_mlbl,
                        heavy_elec_mlbl=self.heavy_elec_mlbl,
                        light_hole_bl=self.light_hole_bl,
                        heavy_hole_bl=self.heavy_hole_bl,
                        light_elec_bl=self.light_elec_bl,
                        heavy_elec_bl=self.heavy_elec_bl)
        np.save(self.savelocation / "emassemergence.npy", data_dct)

    def find_small_masses(self, d):
        keys = ["light_hole_bl", "heavy_hole_bl",
                "light_elec_bl", "heavy_hole_bl"]
        mats = []
        for k in keys:
            for folder, mass in d[k]:
                if abs(mass) < 0.1 and abs(mass) > 0.02:
                    mats.append(folder)

        with open("figures/smallmasses.txt", "w+") as f:
            f.write("\n".join(mats))

    def plot(self, save):
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set_context("paper")

        d = np.load(self.savelocation / "emassemergence.npy", allow_pickle=True).item()

        self.find_small_masses(d)

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        light_hole_mlbl = [(a, b) for a, b in d["light_hole_mlbl"]]
        heavy_hole_mlbl = [(a, b) for a, b in d["heavy_hole_mlbl"]]
        light_elec_mlbl = [(a, b) for a, b in d["light_elec_mlbl"]]
        heavy_elec_mlbl = [(a, b) for a, b in d["heavy_elec_mlbl"]]
        light_hole_changes = [b - a for a, b in light_hole_mlbl]
        heavy_hole_changes = [b - a for a, b in heavy_hole_mlbl]
        light_elec_changes = [b - a for a, b in light_elec_mlbl]
        heavy_elec_changes = [b - a for a, b in heavy_elec_mlbl]
        _, bins1 = np.histogram(light_hole_changes, bins=40)
        _, bins2 = np.histogram(light_elec_changes, bins=40)
        axes[0].hist(light_hole_changes, bins1, alpha=0.75, label="Light hole $\Delta m$")
        axes[0].hist(heavy_hole_changes, bins1, alpha=0.75, label="Heavy hole $\Delta m$")
        C = np.corrcoef(light_hole_changes, heavy_hole_changes)[0, 1]
        axes[0].annotate("$C=" + str(round(C, 4)) + "$", (0.75, 0.7), xycoords="axes fraction")
        axes[1].hist(light_elec_changes, bins1, alpha=0.75, label="Light elec. $\Delta m$")
        axes[1].hist(heavy_elec_changes, bins1, alpha=0.75, label="Heavy elec. $\Delta m$")
        C = np.corrcoef(light_elec_changes, heavy_elec_changes)[0, 1]
        axes[1].annotate("$C=" + str(round(C, 4)) + "$", (0.75, 0.7), xycoords="axes fraction")
        axes[0].legend()
        axes[1].legend()
        axes[0].set_xlabel("Hole mass change [$m_e$]")
        axes[1].set_xlabel("Electron mass change [$m_e$]")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")
        
        if save:
            plt.savefig("figures/emassemergence.pdf", bbox_inches="tight")
        else:
            plt.show()


def get_light_heavy(masses):
    light = min(masses, key=lambda t:abs(t))
    heavy = max(masses, key=lambda t:abs(t))
    return light, heavy


class EmassSwitching(Analysis):
    def __init__(self, savelocation):
        super().__init__(savelocation)
        self.light_hole_pairs = []
        self.heavy_hole_pairs = []
        self.light_elec_pairs = []
        self.heavy_elec_pairs = []

        self.candidates = []

    def run_monolayer(self, monolayer_folder):
        pass

    def run_bilayers(self, bilayer_folders):
        pass

    def run_switchables(self, switchables):
        for mat1, mat2 in switchables:
            data1 = get_bilayer_data(mat1, "results-asr.emasses@validate.json", False)
            data2 = get_bilayer_data(mat2, "results-asr.emasses@validate.json", False)
            if data1 is None or data2 is None:
                continue

            try:
                holes1, elecs1 = get_holes_elecs(data1)
                holes2, elecs2 = get_holes_elecs(data2)
            except KeyError as e:
                continue

            if any(hole > 0.0 for hole in holes1) or any(hole > 0.0 for hole in holes2):
                continue

            if any(elec < 0.0 for elec in elecs1) or any(elec < 0.0 for elec in elecs2):
                continue
                
            if (any(abs(hole) > 100 for hole in holes1)
                or any(abs(elec) > 100 for elec in elecs1)):
                continue
            if (any(abs(hole) > 100 for hole in holes2)
                or any(abs(elec) > 100 for elec in elecs2)):
                continue

            ddip = abs(get_dipole_change(mat1, mat2))
            
                
            light_hole1, heavy_hole1 = get_light_heavy(holes1)
            light_elec1, heavy_elec1 = get_light_heavy(elecs1)
            light_hole2, heavy_hole2 = get_light_heavy(holes2)
            light_elec2, heavy_elec2 = get_light_heavy(elecs2)

            self.light_hole_pairs.append((light_hole1, light_hole2, ddip, mat1, mat2))
            self.heavy_hole_pairs.append((heavy_hole1, heavy_hole2, ddip, mat1, mat2))
            self.light_elec_pairs.append((light_elec1, light_elec2, ddip, mat1, mat2))
            self.heavy_elec_pairs.append((heavy_elec1, heavy_elec2, ddip, mat1, mat2))

    def find_candidates(self):
        from ase.db import connect
        c2db = connect("/home/niflheim2/cmr/databases/c2db/c2db.db")
        light_hole_pairs = self.light_hole_pairs
        heavy_hole_pairs = self.heavy_hole_pairs
        light_elec_pairs = self.light_elec_pairs
        heavy_elec_pairs = self.heavy_elec_pairs

        def check_for_ehull(pairs):
            candidates = []
            for m1, m2, ddip, mat1, mat2 in pairs:
                ml_uid = mat1.split("/")[-2]
                row = c2db.get(f"uid={ml_uid}")
            
                if row.ehull < 0.05:
                    dbid = next((x for x in ["doi", "icsd_id", "cod_id"] if hasattr(row, x)), "none")
                    candidates.append((ml_uid, dbid, m1, m2, ddip, mat1, mat2))
            return candidates

        candidates = check_for_ehull(light_hole_pairs) + check_for_ehull(heavy_hole_pairs)
        candidates += check_for_ehull(light_elec_pairs) + check_for_ehull(heavy_elec_pairs)

        self.candidates = candidates

    def finalize_monolayer(self):
        pass

    def save(self):
        self.find_candidates()

        data_dct = dict(light_hole_pairs=self.light_hole_pairs,
                        heavy_hole_pairs=self.heavy_hole_pairs,
                        light_elec_pairs=self.light_elec_pairs,
                        heavy_elec_pairs=self.heavy_elec_pairs,
                        candidates=self.candidates)
        np.save(self.savelocation / "emassswitch.npy", data_dct)

    def get_by_dip(self, d, dip):
        light_hole_pairs = [(a, b) for a, b, ddip, _, _ in d["light_hole_pairs"]
                            if ddip > dip]
        heavy_hole_pairs = [(a, b) for a, b, ddip, _, _ in d["heavy_hole_pairs"]
                            if ddip > dip]
        light_elec_pairs = [(a, b) for a, b, ddip, _, _ in d["light_elec_pairs"]
                            if ddip > dip]
        heavy_elec_pairs = [(a, b) for a, b, ddip, _, _ in d["heavy_elec_pairs"]
                            if ddip > dip]

        return light_hole_pairs, heavy_hole_pairs, light_elec_pairs, heavy_elec_pairs

    def plot_it(self, save, pairs, savename):
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set_context("paper")
        light_hole_pairs, heavy_hole_pairs, light_elec_pairs, heavy_elec_pairs = pairs
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

        assert len(light_hole_pairs) == len(heavy_hole_pairs)
        assert len(light_elec_pairs) == len(heavy_elec_pairs)
        light_hole_changes = [abs(b - a) for a, b in light_hole_pairs]
        heavy_hole_changes = [abs(b - a) for a, b in heavy_hole_pairs]
        light_elec_changes = [abs(b - a) for a, b in light_elec_pairs]
        heavy_elec_changes = [abs(b - a) for a, b in heavy_elec_pairs]
        _, bins1 = np.histogram(light_hole_changes, bins=40)
        _, bins2 = np.histogram(light_elec_changes, bins=40)
        axes[0].hist(light_hole_changes, bins1, alpha=0.75, label="Light hole $\Delta m$")
        axes[0].hist(heavy_hole_changes, bins1, alpha=0.75, label="Heavy hole $\Delta m$")
        C = np.corrcoef(light_hole_changes, heavy_hole_changes)[0, 1]
        axes[0].annotate("$C=" + str(round(C, 4)) + "$", (0.65, 0.7), xycoords="axes fraction")
        axes[1].hist(light_elec_changes, bins1, alpha=0.75, label="Light elec. $\Delta m$")
        axes[1].hist(heavy_elec_changes, bins1, alpha=0.75, label="Heavy elec. $\Delta m$")
        C = np.corrcoef(light_elec_changes, heavy_elec_changes)[0, 1]
        axes[1].annotate("$C=" + str(round(C, 4)) + "$", (0.65, 0.7), xycoords="axes fraction")
        axes[0].legend()
        axes[1].legend()
        axes[0].set_xlabel("Hole mass change [$m_e$]")
        axes[1].set_xlabel("Electron mass change [$m_e$]")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")

        if save:
            plt.savefig(f"figures/{savename}.pdf", bbox_inches="tight")
        else:
            plt.show()

    def extract_it(self, save, pairs, savename):
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set_context("paper")

        light_hole_pairs, heavy_hole_pairs, light_elec_pairs, heavy_elec_pairs = pairs
        light_hole_changes = [abs(b - a)/min(abs(a), abs(b))*100 for a, b in light_hole_pairs]
        heavy_hole_changes = [abs(b - a)/min(abs(a), abs(b))*100 for a, b in heavy_hole_pairs]
        light_elec_changes = [abs(b - a)/min(a, b)*100 for a, b in light_elec_pairs]
        heavy_elec_changes = [abs(b - a)/min(a, b)*100 for a, b in heavy_elec_pairs]

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        
        _, bins1 = np.histogram(light_hole_changes)
        bins1 = np.linspace(min(min(light_hole_changes), min(heavy_hole_changes)),
                            max(max(light_hole_changes), max(heavy_hole_changes)), 40)
        bins2 = np.linspace(min(min(light_elec_changes), min(heavy_elec_changes)),
                            max(max(light_elec_changes), max(heavy_elec_changes)), 40)
        # _, bins2 = np.histogram(light_elec_changes)
        axes[0].hist(light_hole_changes, bins1, alpha=0.75, label="Light hole $\Delta m$")
        axes[0].hist(heavy_hole_changes, bins1, alpha=0.75, label="Heavy hole $\Delta m$")
        C = np.corrcoef(light_hole_changes, heavy_hole_changes)[0, 1]
        axes[0].annotate("$C=" + str(round(C, 4)) + "$", (0.65, 0.7), xycoords="axes fraction")
        axes[1].hist(light_elec_changes, bins2, alpha=0.75, label="Light elec. $\Delta m$")
        axes[1].hist(heavy_elec_changes, bins2, alpha=0.75, label="Heavy elec. $\Delta m$")
        C = np.corrcoef(light_elec_changes, heavy_elec_changes)[0, 1]
        axes[1].annotate("$C=" + str(round(C, 4)) + "$", (0.65, 0.7), xycoords="axes fraction")

        axes[0].legend()
        axes[1].legend()

        axes[0].set_xlabel("Hole mass change [%]")
        axes[1].set_xlabel("Electron mass change [%]")
        axes[0].set_ylabel("Count")
        axes[1].set_ylabel("Count")


        if save:
            plt.savefig(f"figures/{savename}.pdf", bbox_inches="tight")
        else:
            plt.show()

    def analyse_candidates(self, d, save, dip):
        candidates = d["candidates"]
        cand_dip = [(a, c, m1, m2, ddip, mat1, mat2)
                    for a, c, m1, m2, ddip, mat1, mat2 in candidates
                    if abs(ddip) > dip and abs(m1 - m2)/min(abs(m1), abs(m2)) > 1.0]
        mls = set([t[0] for t in cand_dip])
        if save:
            np.save("figures/emassswitch_candidates.txt", cand_dip)
        else:
            for x in cand_dip:
                print(x)

            print("###############")
            print("###############")
            for x in mls:
                print(r"\item ", x)
            print(f"Number of candidates: {len(cand_dip)}")
            print(f"Number of monolayers: {len(mls)}")

    def plot(self, save):
        d = np.load(self.savelocation / "emassswitch.npy", allow_pickle=True).item()

        self.analyse_candidates(d, save, 1e-13)
        return

        pairs = self.get_by_dip(d, -1)
        self.plot_it(save, pairs, "emassswitch_nodip")

        pairs = self.get_by_dip(d, 1e-13)
        # self.plot_it(save, pairs, "emassswitch_dip")
        self.extract_it(save, pairs, "emassswitch_dip")


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
