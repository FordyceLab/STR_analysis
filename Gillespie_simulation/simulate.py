import numpy as np
import argparse
import tqdm
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.stats

'''
This file contains functions to run a gillespie simulation for TF search with or 
without repeats flanking a motif. Sensitivity analysis functions are also included.
The simulation and the sensitivity analysis functions are intended to be parallelized
to speed up computation. See Gillespie method for information about states and 
equations defining rates, affinities, and reaction likelihoods.
'''

# Global variables
DT_INDEX = 0
TIME_INDEX = 1
LOCAL_INDEX = 2
TEST_INDEX = 3
MOTIF_INDEX = 4
FLANK_INDEX = 5
DNA_INDEX = 6
RXN_INDEX = 7

SIM_TIME = 100000
MAX_TIME = 100000
PLOT_TIME = 100

F_RANDOM_FLANK = 2.23e-1
F_REPEAT_FLANK = 7.89e-1
F_MOTIF = 8.46e-1

NM_per_BP = 0.34 # nm / bp
AVOGADROS_NUMBER = 6.022e23 # molecules/mol
NUC_VOL = 3	#um^3
PERSISTENCE_LEN = 0.05 # um
R_G = PERSISTENCE_LEN/2 # um
TF_DIFFUSION = 1 # um^2/s
N_TF = 2600
SEGMENT_LEN = 150

REPEAT_FACTOR = 0.1
RANDOM_FACTOR = 0.01


# parses input and runs appropriate simulation or sensitivity analysis
def main():
	parser = argparse.ArgumentParser(description='Get run number and sensitivity analysis target')
	parser.add_argument('run_num', type=str, help='run number')
	parser.add_argument('target', type=str, help='sensitivity analysis target variable')
	parser.add_argument("-y0", type=int, nargs="+", default=[N_TF, 0, 0, 0, 1])
	args = parser.parse_args()
	factors = np.geomspace(1e-2, 1e2, 10)

	# these targets are all set up for parallelization
	if args.target == 'n_TF':
		n_tf_sensitivity(args.target, args.run_num, args.y0, factors)
	if args.target == 'k_on':
		k_on_sensitivity(args.target, args.run_num, args.y0, factors)
	if args.target == 'simulation':
		factors = np.geomspace(1e-2, 1e2)
		simulation(factors, args.run_num, args.y0)

	# these targets are set up to run locally for quick testing
	if args.target == 'one_simulation':
		one_simulation(args.target, args.y0)
	if args.target == 'mfpt_simulation':
		mfpt_simulation(args.target, args.run_num, args.y0, F_REPEAT_FLANK/F_MOTIF, verbose=True)
	if args.target == 'steady_state':
		steady_state(args.target, args.y0)
	if args.target == 'f_flank':
		f_flank(args.target, args.run_num, args.y0)

	print('Run completed')

def f_flank(target, run_num, y0):
	ratio = np.geomspace(1e-2, 1e2, 20)
	stats = np.zeros((len(ratio), 4))
	for i, r in enumerate(ratio):
		stats[i,:] = mfpt_simulation(target, run_num, y0, r, ifReturn=True)

	fig = plt.figure(figsize=(18, 10))
	ax = fig.add_subplot()
	ax.errorbar(x=ratio, y=stats[:, 1],
				yerr=stats[:, 3], marker='o', capsize=5)
	ax.set_xscale('log')
	ax.set_ylabel('time (sec)', fontsize=20)
	ax.set_xlabel('Ratio of f_flank to f_motif', fontsize=20)
	ax.tick_params(labelsize=14)
	# plt.legend(['Random flanks', 'Repeat flanks'])
	plt.savefig('f_flank.pdf')
	np.save('mfpt_')

# returns kinetic parameters for gillespie simulation
def get_k_array(ratio, k_on_MAX=2.67e-4):
	k_on_max = k_on_MAX*1e9 # 1/Ms
	k_off_u_m = 3.17e-1
	k_off_u_f = 3.73e-1
	f_motif = 8.46e-1
	k_off_M = 1e6

	kLT = k_on_max
	kTL = k_off_M

	kMT = k_off_u_m # this is koff, 1/s
	kTM = f_motif*k_off_M  # detailed balance - this is kon, 1/Ms

	kFT = k_off_u_f # this is koff, 1/s
	kTF = f_motif*ratio*k_off_M  # detailed balance - this is kon, 1/Ms

	return kLT, kTL, kTM, kTF, kMT, kFT


# Gillespie simulation of TFsearch
def simulate_tf_search(sim_T, max_T, y0, k_array, mfpt_only=False):
	# initialization
	n = AVOGADROS_NUMBER
	kLT, kTL, kTM, kTF, kMT, kFT = k_array

	stored_data = False
	first_passage = False
	first_passage_time = -1
	rxn = -1
	t = 0
	i = 1
	y = y0.copy()
	n_rows = 100000
	sim_data = np.zeros((n_rows, len(y0) + 3))
	sim_data[0] = np.hstack([0, t, y, rxn])
	sim_data_occ = np.asarray(sim_data)

	local_index = LOCAL_INDEX - 2
	test_index = TEST_INDEX - 2
	tf_motif_index = MOTIF_INDEX - 2
	tf_flank_index = FLANK_INDEX - 2
	dna_index = DNA_INDEX - 2

	# maps a given reaction to how molecules should be updated (subtract from, add to)
	a_mapping = [([local_index, dna_index], [test_index]),		# entering NS state
				 ([test_index], [local_index, dna_index]),		# back to local state
				 ([test_index], [tf_motif_index]),	# binding to motif
				 ([test_index], [tf_flank_index]),	# binding to flanks
				 ([tf_motif_index], [test_index]),	# unbinding from motif
				 ([tf_flank_index], [test_index]),	# unbinding from flanks
				 ]

	while t < max_T:
		# calculates likelihood of each reaction
		aLT = kLT / (NUC_VOL*1e-15) / n * y[local_index] * y[dna_index]
		aTL = kTL * y[test_index]
		aTM = kTM * y[test_index]
		aTF = kTF * y[test_index]
		aMT = kMT * y[tf_motif_index]
		aFT = kFT * y[tf_flank_index]

		# calculates tau
		a = [aLT, aTL, aTM, aTF, aMT, aFT]
		a_tot = np.sum(a)
		tau = -np.log(np.random.rand()) / a_tot
		t = t + tau

		# chooses next reaction
		rand = np.random.rand()
		for j in range(len(a)):
			if rand <= (np.sum(a[:j+1]) / a_tot):
				rxn = j
				# updates molecule counts
				idx_from, idx_to = a_mapping[j]
				for idx in idx_from:
					y[idx] -= 1
				for idx in idx_to:
					y[idx] += 1
				break

		# checks to see if dna is bound for the first time
		if (y[tf_motif_index] == 1 or y[tf_flank_index] == 1) and not first_passage:
			first_passage_time = t
			first_passage = True

		# allocates more space if needed so that sim_data is not stored dynamically
		if i >= n_rows:
			sim_data = np.vstack((sim_data, np.zeros((n_rows, len(y0) + 3))))
			n_rows += n_rows

		# updates sim_data
		sim_data[i] = np.hstack([tau, t, y, rxn])
		i += 1

		# stores data for the amount of simulation time (for calculating occupancy)
		if t >= sim_T and not stored_data:
			sim_data_occ = np.asarray(sim_data)
			sim_data_occ = sim_data_occ[:np.argmax(sim_data_occ[:, 1]) + 1, :]
			stored_data = True

		# returns if only mfpt is needed
		if mfpt_only and first_passage:
			sim_data_occ = np.asarray(sim_data)
			sim_data_occ = sim_data_occ[:np.argmax(sim_data_occ[:, 1]) + 1, :]
			return sim_data_occ, first_passage_time

		# returns if mfpt and simulation data is stored
		if stored_data and first_passage:
			return sim_data_occ, first_passage_time

	# returns an MFPT 10 times the maximum if it did not occur within maximum time
	return sim_data_occ, max_T*10


# sensitivity analysis for number of TFs
def n_tf_sensitivity(target, run_num, y0, factors):
	TF_numbers = get_run_vars(target)
	first_passage, mean_occupancy_mot, mean_occupancy_flanks, mean_occupancy_local = initialize_storage(
		len(factors), len(TF_numbers))
	for i, factor in enumerate(factors):
		for j, n_TF in enumerate(TF_numbers):
			k_array = get_k_array(factor)
			y0[0] = n_TF
			sim_data, first_passage_time = simulate_tf_search(SIM_TIME, MAX_TIME, y0, k_array)
			first_passage[i, j] = first_passage_time
			mean_occupancy_mot[i, j] = compute_mean_occupancy(sim_data, MOTIF_INDEX)
			mean_occupancy_flanks[i, j] = compute_mean_occupancy(sim_data, FLANK_INDEX)
			mean_occupancy_local[i, j] = compute_mean_occupancy(sim_data, MOTIF_INDEX, FLANK_INDEX)
	save_output(target, run_num, first_passage, mean_occupancy_mot, mean_occupancy_flanks,
					mean_occupancy_local)

# sensitivity analysis for number of TFs
def k_on_sensitivity(target, run_num, y0, factors):
	k_on_numbers = get_run_vars(target)
	first_passage, mean_occupancy_mot, mean_occupancy_flanks, mean_occupancy_local = initialize_storage(
		len(factors), len(k_on_numbers))
	for i, factor in enumerate(factors):
		for j, k_on in enumerate(k_on_numbers):
			k_array = get_k_array(factor, k_on_MAX=k_on)
			sim_data, first_passage_time = simulate_tf_search(SIM_TIME, MAX_TIME, y0, k_array)
			first_passage[i, j] = first_passage_time
			mean_occupancy_mot[i, j] = compute_mean_occupancy(sim_data, MOTIF_INDEX)
			mean_occupancy_flanks[i, j] = compute_mean_occupancy(sim_data, FLANK_INDEX)
			mean_occupancy_local[i, j] = compute_mean_occupancy(sim_data, MOTIF_INDEX, FLANK_INDEX)
	save_output(target, run_num, first_passage, mean_occupancy_mot, mean_occupancy_flanks,
					mean_occupancy_local)


# test for the occupancy function to make sure it is working correctly
def test_occupancy_function():
	sim_data = np.asarray(
		[[0, 0, 0, 0],[4000, 1, 0, 0], [0.25, 0, 1, 0], [0.25, 1, 0, 0], [0.25, 0, 1, 0], [0.25, 0, 0, 1], [0.25, 0, 1, 0]])
	print(compute_mean_occupancy(sim_data, 2, 0))
	print(compute_mean_occupancy(sim_data, 3, 0))
	print(compute_mean_occupancy(sim_data, 2, 0, 1))
	print(compute_mean_occupancy(sim_data, 1, 0))


# runs baseline gillespie simulation
def simulation(factors, run_num, y0, core_affinity=1e-7):
	# initialize
	first_passage = np.zeros(len(factors))
	mean_occupancy_mot = np.zeros(len(factors))
	mean_occupancy_flanks = np.zeros(len(factors))
	mean_occupancy_local = np.zeros(len(factors))

	# loop through all factors and run simulation
	for j, factor in enumerate(factors):
		k_array = get_k_array(factor)
		sim_data, first_passage_time = simulate_tf_search(SIM_TIME, MAX_TIME, y0, k_array)
		first_passage[j] = first_passage_time
		print('factor: ', factor, ' first passage: ', first_passage_time)
		mean_occupancy_mot[j] = compute_mean_occupancy(sim_data, MOTIF_INDEX)
		mean_occupancy_flanks[j] = compute_mean_occupancy(sim_data, FLANK_INDEX)
		mean_occupancy_local[j] = compute_mean_occupancy(sim_data, MOTIF_INDEX, FLANK_INDEX)

		np.save('simulation_output/first_passage_' + run_num + '.npy', first_passage)
		np.save('simulation_output/mean_occupancy_mot_' + run_num + '.npy', mean_occupancy_mot)
		np.save('simulation_output/mean_occupancy_flanks_' + run_num + '.npy', mean_occupancy_flanks)
		np.save('simulation_output/mean_occupancy_local_' + run_num + '.npy', mean_occupancy_local)


# runs one simulation, prints information and plots local, flank and motif bound tfs over time
def one_simulation(target, y0):
	k_array_rpt = get_k_array(F_REPEAT_FLANK/F_MOTIF)
	print('k_array_rpt: ', k_array_rpt)
	sim_data_rpt, first_passage = simulate_tf_search(SIM_TIME, MAX_TIME, y0, k_array_rpt)
	print('rpt first passage: ', first_passage)

	k_array_rand = get_k_array(F_RANDOM_FLANK/F_MOTIF)
	print('k_array_rand: ', k_array_rand)
	sim_data_rand, first_passage = simulate_tf_search(SIM_TIME, MAX_TIME, y0, k_array_rand)
	print('rand first passage: ', first_passage)
	print_occupancy_info(sim_data_rand, sim_data_rpt)

	plot_time = PLOT_TIME
	# plot_single_molecule_trace(sim_data_rpt[:, TIME_INDEX], sim_data_rpt[:, MOTIF_INDEX],
	# 						   sim_data_rand[:, TIME_INDEX], sim_data_rand[:, MOTIF_INDEX],
	# 						   plot_time, 'bound to motif', 'simulation_tfs_motif', target)
	# plot_single_molecule_trace(sim_data_rpt[:, TIME_INDEX], sim_data_rpt[:, FLANK_INDEX],
	# 						   sim_data_rand[:, TIME_INDEX], sim_data_rand[:, FLANK_INDEX],
	# 						   plot_time, 'bound to flanks', 'simulation_tfs_flanks', target)
	# plot_single_molecule_trace(sim_data_rpt[:, TIME_INDEX], sim_data_rpt[:, LOCAL_INDEX] + sim_data_rpt[:, TEST_INDEX],
	# 						   sim_data_rand[:, TIME_INDEX], sim_data_rand[:, LOCAL_INDEX] + sim_data_rand[:, TEST_INDEX],
	# 						   plot_time, 'in local volume', 'simulation_tfs_local', target)

	plot_trace_all_states(sim_data_rpt, sim_data_rand, plot_time, target)

# runs simulation of only mfpt
def mfpt_simulation(target, run_num, y0, ratio, plot=False, verbose=True, ifReturn = False):
	# initialization
	run_num = int(run_num)
	rpt_fpt = np.zeros(run_num)
	rand_fpt = np.zeros(run_num)
	rpt_mfpt_mode = np.zeros(2)
	rand_mfpt_mode = np.zeros(2)

	for i in tqdm.tqdm(range(int(run_num))):
		k_array_rpt = get_k_array(ratio)
		sim_data_rpt, first_passage = simulate_tf_search(1e5, 1e6, y0, k_array_rpt, mfpt_only=True)
		rpt_fpt[i] = first_passage

		k_array_rand = get_k_array(1)
		sim_data_rand, first_passage = simulate_tf_search(1e5, 1e6, y0, k_array_rand, mfpt_only=True)
		rand_fpt[i] = first_passage

		# keeps track of mode for first passage (L-->F-->M or L-->M)
		if sim_data_rpt[:, -1][-1] == 2:
			rpt_mfpt_mode[1] += 1
		else:
			rpt_mfpt_mode[0] += 1
		if sim_data_rand[:, -1][-1] == 2:
			rand_mfpt_mode[1] += 1
		else:
			rand_mfpt_mode[0] += 1

	# prints summary statistics
	if verbose:
		print('rpt mfpt: ', np.round(np.average(rpt_fpt), 3), ', stdev: ', np.round(np.std(rpt_fpt), 3),
			  ', mode: ', rpt_mfpt_mode/np.sum(rpt_mfpt_mode))
		print('rand mfpt: ', np.round(np.average(rand_fpt), 3), ', stdev: ', np.round(np.std(rand_fpt), 3),
			  ', mode: ', rand_mfpt_mode/np.sum(rand_mfpt_mode))
		print('p-value: ', scipy.stats.ttest_ind(rpt_fpt, rand_fpt, equal_var=False)[1])

	# # save data
	# np.save(target + '/rpt_fpt_' + str(run_num) + '.npy', rpt_fpt)
	# np.save(target + '/rand_fpt_' + str(run_num) + '.npy', rand_fpt)

	if plot:
		plot_mfpt(rand_fpt, rpt_fpt, run_num, target)

	if ifReturn:
		return np.average(rand_fpt), np.average(rpt_fpt), np.std(rand_fpt)/np.sqrt(run_num), np.std(rpt_fpt)/np.sqrt(run_num)


def steady_state(target, y0):
	k_array_rpt = get_k_array(REPEAT_FACTOR, y0[-1])
	k_array_rand = get_k_array(RANDOM_FACTOR, y0[-1])

	for k_array in [k_array_rand, k_array_rpt]:
		kNL, kLN, kLM, kLF, kML, kMF, kFL, kFM = k_array
		ss_mat = np.array([
			[-kNL, kLN, 0, 0],
			[0, -kLM*5e-5, kML, 0],
			[0, -kLF*5e-5*100, 0, kFL],
			[0, 0, -kMF, kFM]
		])
		ss_vec = scipy.linalg.null_space(ss_mat) / np.sum(scipy.linalg.null_space(ss_mat))
		print()
		print('matrix: ', ss_mat)
		print('equilibrium vector: ', np.round(ss_vec, 2))


### HELPER FUNCTIONS


# computes the fraction of time that the target (or optionally two targets) is occupied


def compute_fraction_time_occupied(simulation_data, target_idx, target_idx_2=None):
    target_data = simulation_data[:-1, target_idx]
    if target_idx_2 is not None:
        target_data = np.add(target_data, simulation_data[:-1, target_idx_2])
    tot_occupied_time = np.sum(simulation_data[(np.where(target_data > 0)[0] + 1), DT_INDEX])
    return tot_occupied_time/simulation_data[-1, TIME_INDEX]


# computes the average occupancy of the target (or optionally two targets)
def compute_mean_occupancy(simulation_data, target_idx, target_idx_2=None):
	if target_idx_2 is not None:
		return np.average(np.add(simulation_data[:-1, target_idx],
								 simulation_data[:-1, target_idx_2]),
						  weights=simulation_data[1:, DT_INDEX])
	return np.average(simulation_data[:-1, target_idx], weights=simulation_data[1:, DT_INDEX])


# saves simulation results to .npy files
def save_output(target, run_num, first_passage, mean_occupancy_mot,
				mean_occupancy_flanks, mean_occupancy_local):
    np.save('simulation_output/' + target + '_sensitivity_first_passage_' + run_num + '.npy', first_passage)
    np.save('simulation_output/' + target + '_sensitivity_mean_occupancy_mot_' + run_num + '.npy', mean_occupancy_mot)
    np.save('simulation_output/' + target + '_sensitivity_mean_occupancy_flanks_' + run_num + '.npy', mean_occupancy_flanks)
    np.save('simulation_output/' + target + '_sensitivity_mean_occupancy_local_' + run_num + '.npy', mean_occupancy_local)


# creates storage variables for simulation results
def initialize_storage(n_factors, n_var):
    return np.zeros((n_factors, n_var)), np.zeros((n_factors, n_var)), np.zeros((n_factors, n_var)), np.zeros((n_factors, n_var))


# returns variables for sensitivity analyses
def get_run_vars(target):
	if target == 'n_TF':
		return np.geomspace(100, 10000, 10)
	if target == 'k_on':
		return np.geomspace(1e-6, 1e-1, 10)


# prints occupancy information about one simulation
def print_occupancy_info(sim_data_rand, sim_data_rpt):
	print(
		'                               [cellular milleu, local, motif bound, flanks bound, free]')
	print('')
	print('repeat fraction time occupied: ',
		  [compute_fraction_time_occupied(sim_data_rpt, x) for x in
		   list(np.arange(2, 7))])
	print('random fraction time occupied: ',
		  [compute_fraction_time_occupied(sim_data_rand, x) for x in
		   list(np.arange(2, 7))])
	print('repeat mean occupancy: ',
		  [compute_mean_occupancy(sim_data_rpt, x) for x in np.arange(2, 7)])
	print('random mean occupancy: ',
		  [compute_mean_occupancy(sim_data_rand, x) for x in np.arange(2, 7)])
	print('')
	print('repeat mean number of TFs bound in antennae: ',
		  compute_mean_occupancy(sim_data_rpt, 4, 5))
	print('random mean number of TFs bound in antennae: ',
		  compute_mean_occupancy(sim_data_rand, 4, 5))
	print('')
	print('repeat fraction of time with TF bound in antennae: ',
		  compute_fraction_time_occupied(sim_data_rpt, 4, 5))
	print('random fraction of time with TF bound in antennae: ',
		  compute_fraction_time_occupied(sim_data_rand, 4, 5))


# plots example single molecule trace of gillespie algorithm for one state
def plot_single_molecule_trace(x1, y1, x2, y2, plot_time, location, title, target):
	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax1.plot(x1, y1, alpha=0.5, drawstyle='steps-post')
	ax1.set_xlim([0, plot_time])
	ax1.set_title('TFs ' + location + ' with repeat flanks')
	ax1.set_ylabel('# of molecules')
	ax1.set_yticks([0, 1])
	ax1.set_xlabel('time (s)')
	ax2 = fig.add_subplot(212)
	ax2.plot(x2, y2, alpha=0.5, drawstyle='steps-post')
	ax2.set_title('TFs ' + location + ' to motif with random flanks')
	ax2.set_xlim([0, plot_time])
	ax2.set_yticks([0, 1])
	ax2.set_ylabel('# of molecules')
	ax2.set_xlabel('time (s)')
	plt.tight_layout()
	plt.savefig(target + '/' + title + '.pdf')


# plots example single molecule trace of all four states together
def plot_trace_all_states(sim_data_rpt, sim_data_rand, plot_time, target):
	# plot of all four states
	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	state_data = sim_data_rpt[:,
				 DNA_INDEX] + (sim_data_rpt[:,
								  DNA_INDEX]==0 - sim_data_rpt[:,
													 FLANK_INDEX] - sim_data_rpt[:,
																		MOTIF_INDEX]) * 2 + sim_data_rpt[:,
												   FLANK_INDEX]*3 + sim_data_rpt[:,
																	MOTIF_INDEX]*4
	ax1.plot(sim_data_rpt[:, TIME_INDEX], state_data, alpha=0.5, drawstyle='steps-post')
	ax1.set_xlim([0, plot_time])
	ax1.set_title('TFs in local volume with repeat flanks')
	ax1.set_ylabel('# of molecules')
	ax1.set_yticks([1, 2, 3, 4])
	ax1.set_yticklabels(['local', 'test', 'flanks', 'motif'])
	ax1.set_xlabel('time (s)')
	ax2 = fig.add_subplot(212)
	state_data = sim_data_rand[:,
				 DNA_INDEX] + (sim_data_rand[:,
								  DNA_INDEX]==0 - sim_data_rand[:,
													 FLANK_INDEX] - sim_data_rand[:,
																		MOTIF_INDEX]) * 2 + sim_data_rand[:,
													 FLANK_INDEX] * 3 + sim_data_rand[:,
																		MOTIF_INDEX] * 4
	ax2.plot(sim_data_rand[:, TIME_INDEX], state_data, alpha=0.5, drawstyle='steps-post')
	ax2.set_xlim([0, plot_time])
	ax2.set_title('TFs in local volume with random flanks')
	ax2.set_ylabel('# of molecules')
	ax2.set_yticks([1, 2, 3, 4])
	ax2.set_yticklabels(['local', 'test', 'flanks', 'motif'])
	ax2.set_xlabel('time (s)')
	plt.tight_layout()
	plt.savefig(target + '/simulation_summary.pdf')


# plots MFPT for repeat and random with SEM
def plot_mfpt(rand_fpt, rpt_fpt, run_num, target):
	fig, ax = plt.subplots()
	ax.bar([0, 1], [np.average(rand_fpt), np.average(rpt_fpt)], yerr=[np.std(rand_fpt)/np.sqrt(run_num), np.std(rpt_fpt)/np.sqrt(run_num)],
		   align='center', alpha=0.5, ecolor='black', capsize=10)
	# ax.violinplot([rand_fpt, rpt_fpt], [0, 1], showmeans=True, showextrema=True)
	ax.set_ylabel('Mean first passage time (s)')
	ax.set_xticks([0, 1])
	ax.set_xticklabels(['random', 'repeat'])
	ax.set_title('MFPT for random and repeat flanks (N = %d)'%run_num)
	ax.yaxis.grid(True)

	# Save the figure and show
	plt.tight_layout()
	plt.savefig(target + '/MFPT_plot.pdf')


if __name__ == "__main__":
	main()
