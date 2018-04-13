from setuptools import setup
import lab

setup(
    name='lab-modules',
    version=lab.__version__,
    packages=['lab', 'lab.vcf', 'lab.gene', 'lab.peak', 'lab.genome', 'tests', 'tests.vcf', 'tests.gene', 'tests.peak',
              'tests.genome'],
    url='https://github.com/jeongmw/lab-modules',
    license='MIT',
    author='Minwoo Jeong',
    author_email='jmw95217@gmail.com',
    description='Modules for the biological research',
    long_description=open('README.md').read(),
    setup_requires=[
        "pytest >= 3.2.1",
        "numpy >= 1.13.3",
        "matplotlib >= 2.1.0"
    ],
)
