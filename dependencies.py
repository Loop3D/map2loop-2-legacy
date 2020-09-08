import subprocess

deps = []
with open('./dependencies.txt', 'r') as f:
    for line in f:
        deps.append(line.strip())

command = 'conda install -c loop3d -y python=3.7'.split() + deps
subprocess.run(command)
