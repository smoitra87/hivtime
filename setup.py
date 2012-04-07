try : 
    from setuptools import setup
except ImportError :
    from distutils.core import setup


config = {
    'description' : 'HIV time sequence model',
    'author' : 'Subhodeep Moitra',
    'url' : 'https://github.com/smoitra87/...' ,
    'download_url' : 'https://github.com/smoitra87/.../downloads',
    'author_email' : 'subho@cmu.edu',
    'version' : '0.1',
    'install_requires' : ['nose'],
    'packages' : ['hiv-time'],
    'scripts' : [],
    'name' : 'hiv-time'
}

setup(**config)

