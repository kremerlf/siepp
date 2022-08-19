# import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as cp


class Bands:
    """
    Siesta band object, uses the SystemName.bands.

    Functionalities:
    - Plot band -> self.plot()
    - Get gap info -> self.gap()
    """

    def __init__(self, system_name: str) -> None:
        self.system_name = system_name

        # bands file to a list
        with open(self.system_name + '.bands', 'r') as f:
            bands_file = f.read().split()

        # get last element of the band energies
        for i in range(len(bands_file)):
            if len(bands_file[i]) == 1:
                end_bands = i

        # band information
        e_fermi = float(bands_file[0])
        x_range = [float(bands_file[1]), float(bands_file[2])]
        n_bnd = [int(bands_file[5]), int(bands_file[6]), int(bands_file[7])]
        kpt_point = [float(i) for i in bands_file[end_bands + 1 :: 2]]
        kpt_symbol = [
            str(i).replace('Gamma', "'$\mathrm{\Gamma}$'").replace("'", '')
            for i in bands_file[end_bands + 2 :: 2]
        ]

        # returns the band info in one list
        self.info = [e_fermi, x_range, n_bnd, kpt_point, kpt_symbol]

        bands = bands_file[8:end_bands]

        k_path, k_path_e = [], []
        for j in range(0, n_bnd[0] * n_bnd[1]):
            for i in range(0, len(bands), n_bnd[0] * n_bnd[1] + 1):
                k_path.append(float(bands[i]))
            k_path.append(None)
            for i in range(j + 1, len(bands), n_bnd[0] * n_bnd[1] + 1):
                k_path_e.append(float(bands[i]))
            k_path_e.append(None)

        # returns the band in one list
        if n_bnd[1] == 1:
            # x_up, y_up
            self.band_points = [k_path[: len(k_path)], k_path_e[: len(k_path)]]
        else:
            # x_up, y_up, x_dn, y_dn
            self.band_points = [
                k_path[: int(len(k_path) / 2)],
                k_path_e[: int(len(k_path) / 2)],
                k_path_e[int(len(k_path) / 2) :],
            ]

        # removes all auxiliar variabels from the namespace
        del (
            e_fermi,
            x_range,
            n_bnd,
            kpt_point,
            kpt_symbol,
            k_path,
            k_path_e,
            end_bands,
        )

    def gap(self, info: str = 'gap') -> float:
        """
        Get the band gap of the structure and its vbm and cbm.

        If you want all the info pass info='all' to the function,
        otherwise it will give only the gap
        """
        vb_up, cb_up, vb_dn, cb_dn = [], [], [], []

        for i in range(len(self.band_points[0])):
            if (
                self.band_points[1][i] != None
                and self.band_points[1][i] < self.info[0]
            ):
                vb_up += [[self.band_points[0][i], self.band_points[1][i]]]

            if (
                self.band_points[1][i] != None
                and self.band_points[1][i] < self.info[0]
                and self.info[2][1] == 2
            ):
                vb_dn += [[self.band_points[0][i], self.band_points[2][i]]]

        if (
            self.band_points[1][i] != None
            and self.band_points[1][i] > self.info[0]
        ):
            cb_up += [[self.band_points[0][i], self.band_points[1][i]]]

        if (
            self.band_points[1][i] != None
            and self.band_points[1][i] > self.info[0]
            and self.info[2][1] == 2
        ):
            cb_dn += [[self.band_points[0][i], self.band_points[2][i]]]

        if self.info[2][1] == 2:
            vbm = [
                sorted(vb_up, key=lambda x: x[1])[-1],
                sorted(vb_dn, key=lambda x: x[1])[-1],
            ]

            cbm = [
                sorted(cb_up, key=lambda x: x[1])[0],
                sorted(cb_dn, key=lambda x: x[1])[0],
            ]

            gap_info = [
                f'up: {abs(vbm[0][1] - cbm[0][1]):6f}',
                'Direct' if vbm[0][0] - cbm[0][0] == 0 else 'Indirect',
                f'dn: {abs(vbm[1][1] - cbm[1][1]):6f}',
                'Direct' if vbm[1][0] - cbm[1][0] == 0 else 'Indirect',
            ]
        else:
            vbm = sorted(vb_up, key=lambda x: x[1])[-1]
            cbm = sorted(cb_up, key=lambda x: x[1])[0]
            gap_info = [
                f'{abs(vbm[1] - cbm[1]):6f}',
                'Direct' if vbm[0] - cbm[0] == 0 else 'Indirect',
            ]

        if info == 'gap':
            return gap_info
        elif info == 'all':
            return [vbm, cbm, gap_info]

    def plot(self, fermi_level=False, yFromEF=5.5, save=False, spin_pol=True):
        """
        Plot the band structure of the SystemName.
        args:
                fermi_level (False):
                        if True will shit the plot to the Fermi level.
                yFromEF (5.5): takes a float in eV units that is the yrange of the plot.
                        The standart yrange is (e_fermi-5.5, e_fermi+5.5).
                save (False):
                        You can save an eps figure if set to True
                spin_pol (True):
                        If set to False will only plot the spin upp channel,
                        useful if your ran a spin polarized calculation but the system is not polarized

        """
        eF = self.info[0]
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.size'] = 18
        plt.figure(figsize=[6, 5])
        plt.axis([0, float(self.info[1][1]), eF - yFromEF, eF + yFromEF])
        plt.xticks(self.info[3], self.info[4])
        for hS in self.info[3]:
            plt.axvline(hS, color='black', lw=0.5)
        plt.axhline(eF, color='red', ls='--', lw=0.75)
        plt.ylabel('Energy (eV)')
        plt.tight_layout()

        bnd_p = cp(self.band_points)
        if fermi_level == True:
            plt.axis([0, float(self.info[1][1]), eF, abs(eF)])
            plt.axhline(0, color='red', ls='--', lw=0.75)
            for i in range(len(bnd_p[0])):
                for j in 1, 2:
                    if bnd_p[j][i] == None:
                        pass
                    else:
                        bnd_p[j][i] = bnd_p[j][i] - eF

        if self.info[2][1] == 2 and spin_pol == True:
            plt.plot(bnd_p[0], bnd_p[1], label=r'$\uparrow$')
            plt.plot(bnd_p[0], bnd_p[2], label=r'$\downarrow$')
            plt.legend(loc='upper right')
        elif spin_pol == False:
            plt.plot(bnd_p[0], bnd_p[1])

        if save == True:
            plt.savefig(f'band-{self.system_name}.eps', dpi=95)

    def __str__(self):

        to_print = f'E_F(eV) = {self.info[0]:6f}\n'
        to_print += f'Band x_max = {self.info[1][1]:6f}\n'
        to_print += f'nBands = {self.info[2][0]}\n'
        to_print += f'nSpin = {self.info[2][1]}\n'
        to_print += f'nKpts (path) = {self.info[2][2]}\n'
        to_print += f'Gap = {self.gap()[0]} ({self.gap()[1]})\n\t dn: {self.gap()[0]} ({self.gap()[1]})'

        return to_print
