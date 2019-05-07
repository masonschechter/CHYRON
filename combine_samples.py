import os

filepath1 = '/data/users/mscheche/data/Jan_2019/orange_rerun/results/separated_by_barcode/'
filepath2 = '/pub/mscheche/data/Mar_2019/new/run_1/orange_16nt/results/separated_by_barcode/'
outdir = '/pub/mscheche/data/combo/orange/Jan+Run_1/results/separated_by_barcode/'
for filename in os.listdir(filepath1):
	if filename in os.listdir(filepath2):
		with open(f"{outdir+filename}", "w") as file:
			os.system(f"cat {filepath1+filename} {filepath2+filename} >> {outdir+filename}")