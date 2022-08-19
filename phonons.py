import matplotlib.pyplot as plt

from bands import Bands


class Phonons(Bands):
    """
    Reads ands plot phonons from the siesta calculation
    """
	
    def plot(self, save: bool = False) -> None:

        # figure params
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.size'] = 18
        plt.figure(figsize=[10, 9])

        plt.axis(
            [0, float(self.info[1][1]), -50, (self.band_points[1][-2] + 50)]
        )
        plt.xticks(self.info[3], self.info[4])
        plt.ylabel('cm$^{-1}$')
        plt.tight_layout()

        # create vertical lines for the high-symm points
        for hS in self.info[3]:
            plt.axvline(hS, color='black', lw=0.5)

        # plot
        plt.plot(self.band_points[0], self.band_points[1])

        # save file
        if save == True:
            plt.savefig(f'phonon-{self.system_name}.pdf', dpi=95)
