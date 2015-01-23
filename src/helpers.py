import os

# Links to directory on TSCC
tscc_path = "/projects/ps-jcvi/projects/Lasken_Single_Cell/"
tscc_data = os.path.join(tscc_path, 'data')
tscc_code = os.path.join(tscc_path, 'code')

# Links to local directories
local_path = "/home/kunal/tscc_projects/lasken/"
local_data = os.path.join(local_path, 'data')
local_code = os.path.join(local_path, 'code')

def local_to_tscc(path):
    return path.replace(local_path, tscc_path)

def tsccify_command(command):
    '''
    Currently assumes that a command is a string, but can be changed later.
    '''
    return command.replace(local_path, tscc_path)

def get_sample_map():
    original_samples = """1E1HCNP111314
1F1HCNP111314
1G1HCNP111314
1H1HCNP111314
2A1HCNPN111314
2B1HCNPN111314
2C1HCNPN111314
2D1HCNPN111314
2G1HCNPN111314
2H1HCNPN111314
3AHCN100pg111314
3BHCN10pg111314
3DHCN100pg111314
3EHCN10pg111314""".split('\n')

    new_sample_names = """Non-neuronal nucleus-1
Non-neuronal nucleus-2
Non-neuronal nucleus-3
Non-neuronal nucleus-4
Neuronal nucleus-1
Neuronal nucleus-2
Neuronal nucleus-3
Neuronal nucleus-4
Neuronal nucleus-5
Neuronal nucleus-6
Total RNA-100pg-1
Total RNA-10pg-1
Total RNA-100pg-2
Total RNA-10pg-2""".split('\n')

    sample_map = dict((sample[0], sample[1]) for sample in zip(original_samples, new_sample_names))

    return sample_map


color_map = {'Neuronal': 'red', 'Non-neuronal':'green', 'Pooled':'blue'}

def get_sorted_color_map():
    return [color_map['Neuronal']] * 6 + [color_map['Non-neuronal']] * 4 + [color_map['Pooled']] * 4

def get_sorted_order():
    index = ['Total RNA-100pg-2',
 'Total RNA-100pg-1',
 'Total RNA-10pg-2',
 'Total RNA-10pg-1',
 'Non-neuronal nucleus-4',
 'Non-neuronal nucleus-3',
 'Non-neuronal nucleus-2',
 'Non-neuronal nucleus-1',
 'Neuronal nucleus-6',
 'Neuronal nucleus-5',
 'Neuronal nucleus-4',
 'Neuronal nucleus-3',
 'Neuronal nucleus-2',
 'Neuronal nucleus-1']
    return index