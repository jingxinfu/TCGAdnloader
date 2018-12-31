from setuptools import setup


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


def main():
    setup(name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        url='https://github.com/jingxinfu/dataCollection',
        author='Jingxin Fu',
        author_email='jingxinfu.tj@gmail.com',
        packages=PACKAGE,
        install_requires=REQUIRES,
        license=open('LICENSE').read(),
        scripts=['bin/collect'],
        include_package_data=True,
        zip_safe=False)

if __name__ == '__main__':
    main()
