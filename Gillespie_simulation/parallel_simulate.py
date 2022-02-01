import os
import argparse

def run_jobs(num_jobs, target):
	for i in range(num_jobs):
		os.environ["TARGET"] = target
		os.environ["RUN_NUM"] = str(i)
		run_command = 'sbatch simulate.sh ' + str(i)
		os.system(run_command)

def main():
	parser = argparse.ArgumentParser(description='Get number of jobs and sensitivity analysis target')
	parser.add_argument('num_jobs', type=int, help='number of jobs')
	parser.add_argument('target', type=str, help='sensitivity analysis target variable')
	args = parser.parse_args()
	run_jobs(args.num_jobs, args.target)

if __name__ == "__main__":
    main()
