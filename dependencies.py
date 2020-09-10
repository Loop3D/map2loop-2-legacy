import subprocess
import platform

deps = []
with open('./dependencies.txt', 'r') as f:
    for line in f:
        deps.append(line.strip())

command = 'conda install -c loop3d -c conda-forge -y python=3.7'.split() + deps

if platform.system() == "Windows":
    subprocess.run(command, shell=True)
else:
    # Unix like
    subprocess.run(command)