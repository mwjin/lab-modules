from distutils.core import setup

setup(
    name='lab-modules',
    version='1.0.0',
    packages=['lab', 'lab.vcf', 'lab.gene', 'lab.peak', 'lab.genome', 'tests', 'tests.vcf', 'tests.gene', 'tests.peak',
              'tests.genome'],
    url='https://github.com/jeongmw/lab-modules',
    license='MIT',
    author='Minwoo Jeong',
    author_email='jmw95217@gmail.com',
    description='Modules for the biological research',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.13.3",
        "matplotlib >= 2.1.0",
        "pytest >= 3.2.1",
    ],
)
