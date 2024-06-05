#!/usr/bin/python3

from pymol import cmd
import pymol
import glob

SIM_DIR = os.getcwd()
FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]
color_map = {
    'model_01': 'red',
    'model_02': 'blue',
    'model_03': 'green',
    'model_04': 'purple',
    'model_05': 'orange'
}


pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
cmd.set("ray_opaque_background", 1)

for i in pdb_data[:25]:
	rep_name = i.split('/')[-3]
	forcefield = i.split('/')[-2]
	selected = color_map.get(rep_name, 'black')

	md = glob.glob(i + 'md*smooth*xtc')
	temp = glob.glob(i + 'temp*gro')

	if len(md) > 0 and len(temp) > 0:
		cmd.load(temp[0], 'structure')
		cmd.load_traj(md[0], 'structure', state=1, interval=1000)
		cmd.set('all_states', 'on')
		cmd.color(selected, 'structure')
		cmd.ray(300, 300)
		cmd.png(i + 'Ensemble_' + rep_name + '_' + forcefield + '.png')
		cmd.delete('all')

cmd.quit()

fig, axs = plt.subplots(1, 5, figsize=(15, 6))
current_ax = 0

for item in FORCEFIELDS:
	relax_data = glob.glob(SIM_DIR + "/model*/" + item)
	cmd.set("ray_opaque_background", 1)
	for data in sorted(relax_data):
		rep_name = data.split("/")[-2]

		md = glob.glob(data + '/md*smooth*xtc')
		temp = glob.glob(data + '/temp*gro')

		if len(md) > 0 and len(temp) > 0:
			cmd.load(temp[0], rep_name)
			cmd.load_traj(md[0], rep_name, state=1, interval=2000)
			cmd.color('green', rep_name)
			cmd.set('all_states', 'on')
			cmd.ray(300, 300)
	obj = cmd.get_object_list('all')
	for i in obj[1:]:
		cmd.align(obj[0], i, object='aln', transform=0)
	cmd.png(relax_folder + item + '_aligned_fig.png')
	axs[current_ax].imshow(plt.imread(relax_folder + item + '_aligned_fig.png'), aspect='auto')
	axs[current_ax].set_title(item)	
	current_ax += 1
	cmd.delete('all')
                
cmd.quit()

plt.tight_layout()
plt.savefig(relax_folder + 'aligned_fig.png')

pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
cmd.set("ray_opaque_background", 1)

for i in pdb_data[:25]:
	rep_name = i.split('/')[-3]
	forcefield = i.split('/')[-2]
	selected = color_map.get(rep_name, 'black')

	md = glob.glob(i + 'md*smooth*xtc')
	temp = glob.glob(i + 'temp*gro')

	if len(md) > 0 and len(temp) > 0:
		cmd.load(temp[0], rep_name)
		cmd.load_traj(md[0], rep_name, state=1, interval=2000)
		cmd.color('green', rep_name)
		cmd.set('all_states', 'on')
		cmd.ray(300, 300)
	obj = cmd.get_object_list('all')
	for i in obj[1:]:
		cmd.align(obj[0], i, object='aln', transform=0)
	cmd.png(i + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png')	
	cmd.delete('all')

cmd.quit()
