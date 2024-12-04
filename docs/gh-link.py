# Based on https://github.com/sphinx-doc/sphinx/issues/823
# Taken from https://github.com/GalSim-developers/GalSim/blob/2bf7969b88d8c2e9913338afb7fb703b386f0c9d/docs/gh-link.py

#import euclidlike

# The short X.Y version
#version = '.'.join(map(str,euclidlike.__version_info__[:2]))

#blob_url = 'https://github.com/GalSim-developers/GalSim-Euclid-Like/blob/releases/' + version + '/'
blob_url = 'https://github.com/GalSim-developers/GalSim-Euclid-Like/blob/main/'

def gh_link_role(rolename, rawtext, text, lineno, inliner,
                 options={}, content=()):
    from docutils import nodes, utils
    name, path = text.split('<')
    path = path.split('>')[0]
    full_url = blob_url + path
    pnode = nodes.reference(internal=False, refuri=full_url)
    pnode += nodes.literal(name, name, classes=['file'])
    return [pnode], []


def setup(app):
    app.add_role('gh-link', gh_link_role)