import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import re


class Aimd:
	"""
	The Aimd class creates a pandas DataFrame with the data
	contained in the MD run of SIESTA (system_name.MDE).
	Atributes:
		system_name: srt
			use the prefix of the MDE file
		
		sim_step: float
			time step used in the simulation
	Methods:
		self.moving_avg(col_to_move: str, window: int = 500) -> None
			calculates the simple moving average (SMA) from a chosen column
			(col_to_move) from the self.df for a given window.
			
		self.plot_moving_avg_T_and_E(self, window: int = 500) -> None
			calculate and plot the SMA for the Temperature and
			Kohn-Sham energy.
	"""

	def __init__(self, system_name: str, sim_step: float = 1.0) -> None:
		"""
		Sets self.system_name and self.step.
		Then reads the system_name.MDE file and creates the DataFrame.
		"""
		self.system_name = system_name
		self.step = sim_step*(1/1000)
		
		with open(self.system_name+'.MDE', 'r') as f:
			file_info = f.read().split('\n')
		
		# read siesta .out file to add FreeEnergy
		with open(self.system_name+'.out', 'r') as f:
			out_file = f.read().split('\n')
		
		# get FreeEnergy values from .out
		free_energy = []
		for i,line in enumerate(out_file):
			if "SCF Convergence by DM+H criterion" in line:
				free_energy.append(float(list(filter(None,out_file[i-3].split()))[4]))
		
		# remove units and # from header
		file_info[0] = re.sub(r'\([^()]*\)','', file_info[0])
		file_info[0] = re.sub(r'#','', file_info[0])
		
		# organizing the file to create the df
		# siesta throws a empty line at the end, 
		# so we ignore the last item of the file
		temp_list = [file_info[idx].split() for idx,item in enumerate(file_info[:-1])]
		
		# create df and change all to numeric
		self.df = pd.DataFrame(temp_list[1:], columns = temp_list[0])
		for i in self.df.columns:
			self.df[i] = pd.to_numeric(self.df[i])
		
		self.df['FreeE'] = free_energy
		self.df['dFreeE'] = self.df['FreeE'] - self.df['FreeE'][0]
		# scaling some columns
		self.df['Step'] = self.df['Step']*self.step
		self.df['dE_KS'] = self.df['E_KS']-self.df['E_KS'][0]
	
	def moving_avg(self, col_to_move: str, window: int = 500) -> None:
		"""
		Computes the simple moving average (SMA) 
		for a given window of a column of the DataFrame.        
		"""
		self.df[col_to_move+'_m'] = self.df[col_to_move].rolling(window).mean()
    
	def plot_moving_avg_T_and_E(self, window: int = 500) -> None:
		"""
		Plots the values of Temperature and Kohn-Shan energy
		gven by the calculation and it respective SMA.
		"""
        
		# SMA computing
		self.moving_avg('T', window)
		# self.moving_avg('dE_KS', window)
		self.moving_avg('dFreeE', window)
        
		# figure parameters
		# fig = plt.figure(figsize = [10*(1/2.54), 8*(1/2.54)])        
		matplotlib.use('pgf')
		plt.rcParams.update({
		    "font.family": "serif",  # use serif/main font for text elements
		    "text.usetex": True,     # use inline math for ticks
		    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
		    "font.size" : 12,
			"figure.figsize" : [10*(1/2.54), 8*(1/2.54)]
		    # "pgf.preamble": "\n".join([
				# r"\usepackage{url}",            # load additional packages
				# r"\usepackage{unicode-math}",   # unicode math setup
				# r"\setmainfont{Times}",  # serif font via preamble
	    	# ])
		})
		plt.tight_layout()

		# plot ranges
		x_max = self.df['Step'].iloc[-1]
		t_ymin = self.df['T'].min()*0.95
		t_ymax = self.df['T'].max()*1.05
		e_ymin = self.df['dFreeE'].mean()*0.25
		e_ymax = self.df['dFreeE'].mean()*2
		
		# T subplot
		ax1 = plt.subplot(211)
		ax1.axis([0, x_max, t_ymin, t_ymax])
		ax1.plot(self.df['Step'],self.df['T'], color = '#c4c4c4')
		ax1.plot(self.df['Step'],self.df['T_m'],linewidth = 2, 
				 color = 'black', label=f"SMA {window/1000} ps")
		plt.ylabel('T (K)')
		ax1.set_yticks(np.arange(self.df['T'][0]-250, self.df['T'][0]+300, 250, dtype = int))
		ax1.legend()
		ax1.tick_params('x', labelbottom=False)

		# Free Energy subplot
		ax2 = plt.subplot(212)
		ax2.axis([0, x_max, e_ymin, e_ymax])
		ax2.plot(self.df['Step'],self.df['dFreeE'], color = '#c4c4c4')
		ax2.plot(self.df['Step'],self.df['dFreeE_m'], linewidth = 2, color = 'black')
		# ax2.set_yticks(np.arange(0, e_ymax, 0.5))
		plt.ylabel(r'$\Delta$Free E (eV)')
		plt.xlabel("Time (ps)")
		plt.subplots_adjust(hspace = 0)
		# fig.align_ylabels([ax1,ax2])
		plt.savefig(f'aimd_{self.system_name}.pgf',bbox_inches = "tight")
