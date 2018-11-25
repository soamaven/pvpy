from setuptools import setup

setup(
    name='pvpy',
    version='1.04.1.alpha',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        'Topic :: Scientific/Engineering :: Physics',
    ],
    description='Some routines relating photovoltaic solar cells',
    author='Colton Bukowsky',
    author_email='cb@caltech.edu',
    url='https://github.com/soamaven/pvpy',
    download_url='https://github.com/soamaven/pvpy/tarball/1.0.alpha',
    keywords=['solar',
              'photovoltaics',
              'shockley',
              'queisser',
              'yablonovitch',
              'short-circuit current',
              'open-circuit voltage'],
    license='MIT',
    packages=[
        'pvpy',
    ],
    install_requires=[
        'numpy',
        'scipy'
    ],
    include_package_data=True,
    zip_safe=False
)
