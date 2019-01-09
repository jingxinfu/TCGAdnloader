from setuptools import setup
import os,re

NAME = 'dataCollection'
PACKAGE = [NAME]
VERSION = __import__(NAME).__version__
try:
    f = open("requirements.txt", "rb")
    REQUIRES = [i.strip() for i in f.read().decode("utf-8").split("\n")]
except:
    print("'requirements.txt' not found!")
    REQUIRES = []

DESCRIPTION='''
Data Collection
'''


def path_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if filename[0] is not '.':  # filter hidden files
                paths.append(os.path.join(
                    re.sub(NAME+'/', '', path), filename))
    return paths


package_files = path_files(NAME+'/data')

def main():
    setup(name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        url='https://github.com/jingxinfu/dataCollection',
        author='Jingxin Fu',
        author_email='jingxinfu.tj@gmail.com',
        packages=PACKAGE,
        package_data={NAME: package_files},
        install_requires=REQUIRES,
        license=open('LICENSE').read(),
        scripts=['bin/collect'],
        include_package_data=True,
        zip_safe=False)

if __name__ == '__main__':
    main()
