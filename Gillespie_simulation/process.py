import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import scipy

'''
This file contains functions to process parallelized Gillespie simulations or
sensitivity analyses. See Gillespie method for information about the simulation.
'''

F_RANDOM_FLANK = 2.23e-1
F_REPEAT_FLANK = 7.89e-1
F_MOTIF = 8.46e-1


def main():
	parser = argparse.ArgumentParser(
		description='Get number of jobs and sensitivity analysis target')
	parser.add_argument('run_num', type=int, help='number of runs')
	parser.add_argument('target', type=str, help='sensitivity analysis target variable')
	parser.add_argument("-log", type=bool, default=True,
						help='if sensitivity heatmap color scale should be log or not')
	args = parser.parse_args()
	factor = np.geomspace(1e-2, 1e2, 10)

	if args.target == 'simulation':
		factor = np.geomspace(1e-2, 1e2)
		simulation_plot(args.target, args.run_num, factor)
	elif args.target == 'sensitivity_summary':
		plot_sensitivity_summary(args.run_num)
	else:
		first_passage, mean_occupancy_mot, mean_occupancy_flanks, mean_occupancy_local = load_arrays(
			args.run_num, args.target)
		run_vars = get_run_vars(args.target)
		sensitivity_plot(args.target, args.run_num, factor, run_vars, first_passage,
						 mean_occupancy_mot, mean_occupancy_flanks, args.log)
		print('done')


# plots relevant sensitivity plots (heatmaps and correlations)
def sensitivity_plot(target, num_jobs, factor, run_vars, first_passage,
					 mean_occupancy_mot, mean_occupancy_flanks, log):
	# sensitivity heatmaps
	fig = plt.figure(figsize=(20, 6))
	factor = np.around(factor, 3)
	run_vars = np.asarray(run_vars)
	ax = fig.add_subplot(131)
	plot_heatmap(ax, np.mean(first_passage, axis=0), target, run_vars, factor, log,
				 'First passage time')
	ax = fig.add_subplot(132)
	plot_heatmap(ax, np.mean(mean_occupancy_mot, axis=0), target, run_vars, factor, log,
				 'Mean motif occupancy')
	ax = fig.add_subplot(133)
	plot_heatmap(ax, np.mean(mean_occupancy_flanks, axis=0), target, run_vars, factor, log,
				 'Mean flanks occupancy')
	plt.tight_layout()
	plt.savefig(target + '/' + target + '_heatmaps.pdf', dpi=300)

	# correlation plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	spearman_arr = np.zeros(len(run_vars))
	for i, n in enumerate(run_vars):
		r, p = scipy.stats.spearmanr(-np.hstack([factor] * num_jobs), first_passage[:, :, i].flat)
		spearman_arr[i] = r
	np.save(target + '/' + target + '_correlation.npy', spearman_arr)
	ax.plot(run_vars, spearman_arr, marker='o', ls='')
	set_value = get_set_value(target)
	ax.axvline(x=set_value)
	ax.legend(['Spearman rho', 'Set value'])
	if target != 'koff_intercept':
		ax.set_xscale('log')
	ax.set_ylim(0, 1.0)
	ax.set_xlabel(target, fontsize=18)
	ax.set_ylabel('correlation', fontsize=18)
	ax.set_title('Correlation between ' + target + ' and time to first passage')
	ax.tick_params()
	fig.tight_layout()
	plt.savefig(target + '/' + target + '_correlation.pdf', dpi=300)


# plots simulation results
def simulation_plot(target, run_num, factor):
	first_passage = np.zeros((run_num, len(factor)))
	mean_occupancy_mot = np.zeros((run_num, len(factor)))
	mean_occupancy_flanks = np.zeros((run_num, len(factor)))
	mean_occupancy_local = np.zeros((run_num, len(factor)))

	for i in range(run_num):
		first_passage[i, :] = np.load(
			target + '/simulation_output/first_passage_' + str(i) + '.npy')

	fig = plt.figure(figsize=(18, 10))
	ax = fig.add_subplot(111)
	l1 = ax.errorbar(x=factor, y=np.median(first_passage, axis=0),
					 yerr=np.std(first_passage, axis=0) / np.sqrt(run_num), marker='o',
					 color='k', capsize=5)
	ax.set_xscale('log')
	ax.set_ylabel('time (sec)', fontsize=20)
	ax.set_xlabel('ratio of f_flank to f_motif', fontsize=20)
	ax.tick_params(labelsize=14)

	# highlight repeat and random regions for affinity ratios
	# l5 = plt.axvspan(0.2, 0.3, color='black', alpha=0.3)
	# l6 = plt.axvspan(0.9, 1, color='red', alpha=0.3)

	ax.legend((l1),
			   ['mfpt'], fontsize=20, bbox_to_anchor=(1.1, 1),
			   loc='upper left')

	fig.tight_layout()
	plt.savefig(target + '/' + target + '_results_' + str(run_num) + '.pdf', dpi=300)


# plots heatmap for data across affinity ratios
def plot_heatmap(ax, data, target, run_vars, factor, log, title):
	if log:
		sns.heatmap(data, xticklabels=run_vars, norm=LogNorm())
	else:
		sns.heatmap(data, xticklabels=run_vars)
	ax.set_ylabel('affinity ratio\n(flanks/core)', fontsize=18)
	ax.set_xlabel(target, fontsize=18)
	run_vars_ticks = [format(x, ".1e") for x in run_vars]
	factor_ticks = [format(y, ".1e") for y in factor]
	ax.set_xticklabels(run_vars_ticks, rotation=30)
	ax.set_yticklabels(factor_ticks, rotation=0)
	ax.tick_params(labelsize=14)
	ax.set_title(title, fontsize=20)

def plot_sensitivity_summary(run_num):
	if run_num == 0:
		fig, axs = plt.subplots(4, 1, figsize=(15,7))
		titles = ['TF diffusion constant (um^2/s)', 'Motif affinity (M)', 'Number of flanks', 'k_off slope']
		factors = ['diffusion', 'core_affinity', 'n_flanks', 'koff_slope']
		save_as = 'sensitivity_summary_small.pdf'
	else:
		fig, axs = plt.subplots(9, 1, figsize=(15, 16))
		titles = ['TF diffusion constant (um^2/s)', 'Motif affinity (M)', 'Number of flanks',
				  'k_off slope', 'k_off intercept', 'Number of TFs', 'Switching ratio',
				  'DNA concentration (M)', 'Local volume (um^3)']
		factors = 	['diffusion', 'core_affinity', 'n_flanks', 'koff_slope',
					'koff_intercept', 'n_TF', 'switching_ratio', 'DNA_concentration',
					'segment_len']
		save_as = 'sensitivity_summary.pdf'

	cbar_ax = fig.add_axes([0.91, .15, .03, .7])
	for i, target in enumerate(factors):
		correlation = np.load(target + '_correlation.npy')
		sns.heatmap(correlation.reshape(1,10), vmin=0, vmax=0.65, ax=axs[i], cbar_ax=cbar_ax, cbar_kws={'label': 'spearman correlation'})
		run_vars = get_run_vars(target)
		run_vars_ticks = [format(x, ".1e") for x in run_vars]
		set_value = get_set_value(target)
		if target == 'koff_intercept':
			mi = run_vars[0]
			ma = run_vars[-1]
			arrow = (9*(set_value-mi) / (ma - mi)) + 0.5
		else:
			mi = np.log10(run_vars[0])
			ma = np.log10(run_vars[-1])
			arrow = (9*(np.log10(set_value)-mi)/(ma-mi)) + 0.5
		axs[i].axvline(arrow, color='black')
		axs[i].set_xticklabels(run_vars_ticks, fontsize=13)
		axs[i].set_yticks([])
		axs[i].set_title(titles[i], fontsize=18)
	fig.subplots_adjust(hspace=.7)
	plt.savefig(save_as)


# returns variables for sensitivity analyses
def get_run_vars(target):
	if target == 'n_TF':
		return np.geomspace(100, 10000, 10)
	if target == 'k_on':
		return np.geomspace(1e-6, 1e-1, 10)


# returns set variable for simulation
def get_set_value(target):
	if target == 'n_TF':
		return 2600
	if target == 'k_on':
		return 2.67e-4


# initializes storage for first passage and occupancy arrays
def initialize_storage(n_jobs, n_factor, n_var):
    return np.zeros((n_jobs, n_factor, n_var)), np.zeros((n_jobs, n_factor, n_var)), np.zeros((n_jobs, n_factor, n_var)), np.zeros((n_jobs, n_factor, n_var))


# loads and aggregates first passage and occupancy arrays for the given target from each run
def load_arrays(num_jobs, target):
	example = np.load(target + '/simulation_output/' + target + '_sensitivity_first_passage_1.npy')
	first_passage, mean_occupancy_mot, mean_occupancy_flanks, mean_occupancy_local = initialize_storage(
		num_jobs, example.shape[0], example.shape[1])
	for run_num in range(num_jobs):
		first_passage[run_num, :, :] = np.load(
			target + '/simulation_output/' + target + '_sensitivity_first_passage_' + str(run_num) + '.npy')
		mean_occupancy_mot[run_num, :, :] = np.load(
			target + '/simulation_output/' + target + '_sensitivity_mean_occupancy_mot_' + str(run_num) + '.npy')
		mean_occupancy_flanks[run_num, :, :] = np.load(
			target + '/simulation_output/' + target + '_sensitivity_mean_occupancy_flanks_' + str(run_num) + '.npy')
		mean_occupancy_local[run_num, :, :] = np.load(
			target + '/simulation_output/' + target + '_sensitivity_mean_occupancy_local_' + str(run_num) + '.npy')
	return first_passage, mean_occupancy_mot, mean_occupancy_flanks, mean_occupancy_local


if __name__ == "__main__":
    main()