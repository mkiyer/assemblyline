'''
Created on Nov 29, 2010

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# ---- Extension Modules ----------------------------------------------------
def get_extension_modules():
    extensions = []
    # Graph node
    #extensions.append( Extension( "lib.cnode", ["lib/cnode.pyx"] ) )
    # Interval clustering                
    extensions.append( Extension( "assemblyline.lib.bx.cluster", [ "assemblyline/lib/bx/cluster.pyx", "assemblyline/lib/bx/intervalcluster.c"], 
                                  include_dirs=["assemblyline/lib/bx"]) )
    # Interval intersection
    extensions.append( Extension( "assemblyline.lib.bx.intersection", [ "assemblyline/lib/bx/intersection.pyx" ] ) )
    return extensions

def main():
    setup(name = "assemblyline",
          ext_modules = get_extension_modules(),
          author = "Matthew Iyer",
          author_email = "mkiyer@umich.edu",
          description = "Transcriptome meta-assembly pipeline",
          url = "http://assemblyline.googlecode.com",
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__': main()