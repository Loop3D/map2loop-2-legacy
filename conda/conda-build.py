from subprocess import run, Popen, PIPE
from concurrent.futures import ThreadPoolExecutor
import tqdm
import sys
import os


def process(python_version):
    try:

        command = 'conda build --py {} .'.format(
            python_version)
        p = Popen(
            command.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    except Exception as e:
        return e

    output, err = p.communicate()
    return str(err.decode('ascii'))


def build(packages):
    for dirname in packages:
        if not os.path.exists(dirname):
            sys.exit(
                "Make sure each package directory is at the same level as this script.")

    python_versions = [3.6, 3.7, 3.8, 3.9]
    print("Building " + " ".join(packages),
          "for python versions {}".format(python_versions))

    with ThreadPoolExecutor(15) as executor:
        output = list(
            tqdm.tqdm(executor.map(process, python_versions), total=len(packages)))

        for message in output:
            lines = message.split('\\n')
            for line in lines:
                print(line)


def upload(root_buildpath):

    # Create OSX packages from linux versions
    output = ''
    command = "conda convert -p osx-64 {}*tar.bz2 -o {} -f".format(
        root_buildpath + 'linux-64/', root_buildpath)
    p = Popen(
        command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output += p.communicate()[0].decode('ascii')
    output += p.communicate()[1].decode('ascii')
    command = "conda convert -p win-64 {}*tar.bz2 -o {} -f".format(
        root_buildpath + 'linux-64/', root_buildpath)
    p = Popen(
        command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output += p.communicate()[0].decode('ascii')
    output += p.communicate()[1].decode('ascii')

    command = "anaconda upload {} {} {} --force".format(
        root_buildpath + 'linux-64/*.tar.bz2',
        root_buildpath + 'osx-64/*.tar.bz2',
        root_buildpath + 'win-64/*tar.bz2',
    )
    p = Popen(
        command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output += p.communicate()[0].decode('ascii')
    output += p.communicate()[1].decode('ascii')
    print(output)


if __name__ == "__main__":

    packages = [
        '.'
    ]

    command = 'which conda'
    p = Popen(
        command.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    std, err = p.communicate()
    root_buildpath = "/".join(std.decode('ascii').split('/')
                              [:-2]) + '/conda-bld/'

    build(packages)
    # upload(root_buildpath)
