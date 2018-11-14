from setuptools import setup
import lab

setup(
    name='lab-modules',
    version=lab.__version__,
    packages=['lab', 'lab.vcf', 'lab.gene', 'lab.bed', 'lab.genome', 'tests', 'tests.vcf', 'tests.gene', 'tests.bed',
              'tests.genome'],
    url='https://github.com/jeongmw/lab-modules',
    license='MIT',
    author='Minwoo Jeong',
    author_email='jmw95217@gmail.com',
    description='Modules for the biological research',
    long_description=open('README.md').read(),
    install_requires=[
        "pytest >= 3.2.1",
        "numpy >= 1.13.3",
        "matplotlib >= 2.1.0",
        "scipy >= 1.0.0",
        "statsmodels >= 0.9.0"
    ],
)
