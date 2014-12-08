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