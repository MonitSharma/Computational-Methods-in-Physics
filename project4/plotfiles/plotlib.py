import numpy as np
import matplotlib.pyplot as plt

#Fetches files from datafiles folder
#Saves plots to Figures with '_figure.pdf' surfix
def makeplots(plotfor, x_label, y_label,plot_type):

	if plot_type == 'same':
		for filename in plotfor:
			data = np.loadtxt('../datafiles/'+filename, skiprows=1)
			x = np.array(data[:,0])
			y = np.array(data[:,1])


			fig, ax = plt.subplots(figsize = (6, 5))
			plt.subplots_adjust(
			top=0.95,
			bottom=0.15,
			left=0.15,
			right=0.98,
			hspace=0.2,
			wspace=0.2
			)
			plt.plot(x,y,label=filename.split('.')[0])
			#ax.set_title("Training mse")
			ax.set_ylabel(y_label)
			ax.set_xlabel(x_label)
			plt.grid()
			plot_name = filename.split('.')[0] +"_figure.pdf"
			plt.legend()
			plt.savefig('../figures/'+plot_name)
		
	if plot_type =='seperate':
		fig, ax = plt.subplots(figsize = (6, 5))
		plt.subplots_adjust(
		top=0.95,
		bottom=0.15,
		left=0.15,
		right=0.98,
		hspace=0.2,
		wspace=0.2
		)
		for filename in plotfor:
				data = np.loadtxt('../datafiles/'+filename, skiprows=1)
				x = np.array(data[:,0])
				y = np.array(data[:,1])
				plt.plot(x,y,label=filename.split('.')[0])

		ax.set_ylabel(y_label)
		ax.set_xlabel(x_label)
		plt.grid()
		plot_name = 'combined_'
		for filename in plotfor:
			plot_name += filename.split('.')[0] + '_'
		plot_name += "_figure.pdf"

		plt.legend()
		plt.savefig('../figures/'+plot_name)