#
# happyapj.py -- Copy tex files to a new directory with names to make apj happy.
#
# Written by Jonathan Coles <jonathan@physik.uzh.ch>
# This software is in the public domain.
#
# This program will first attempt to create a directory named 'apj' in the
# current directory. If it fails it will try 'apj-1', 'apj-2', etc., up to
# 'apj-10'. Then it will copy the files listed below into that directory,
# renaming them to conform with the apj guidelines. It will then modify
# the copied .tex file by replacing references to the old filenames with
# the new filenames.
#
# List the files to copy here. Figures should be listed in the order they
# appear in the document. If multiple files appear in a single figure they
# should be grouped as a sublist. For example myfig1.eps appears in a single
# figure, while myfig2.eps and myfig3.eps are part of the same figure:
#
# files = ['myfig1.eps', ['myfig2.eps', 'myfig3.eps']]
#
# There are a certain number of assumptions made: 
#   (1) There is only one .tex file 
#   (2) There is only one .bbl file (and it is already made)
#   (3) Files referenced in the tex file include the complete name with the
#       extension.
#
# Files ending in .tex, .bbl, .eps, .ps, .pdf, .png, .jpg, and .jpeg are
# considered figures and a renamed accordingly.  Any files listed below that
# are not recognized will be copied without being renamed.

files = ['ms.tex',
         'ms.bbl',
         ['recon.3153.14.0.eps', 'recon.3153.14.1.eps'],
         'likelihood.eps',
         'fullsim-t2.eps',
         ['recons-ppm.eps', 'recons-mass.eps'],
         ['recons-area.eps', 'recons-area-v-mass.eps'],
         'apj-jour.bib',
         'apj.bst'
        ]

#-------------------------------------------------------------------------------
# DO NOT EDIT BELOW THIS SECTION.
#-------------------------------------------------------------------------------

conversions = {}

def newname(file, ext, mod):
    """Returns an apj-acceptable name for a given file."""
    if ext == '.tex':
        return 'ms', mod
    elif ext == '.bbl':
        return 'ms', mod
    elif ext in ['.eps', '.ps', '.pdf', '.png', '.jpg', '.jpeg']:
        return 'figure%s' % mod[0], mod[1:]
    else:
        return file, mod

def convert_files(files, mod=[]):
    """Copy files to dest and rename them appropriately for apj."""
    for f in files:
        if isinstance(f, list):
            convert_files(f, map(lambda x: mod[0]+x, ascii_letters[:26]))
            mod = mod[1:]
        else:
            file, ext = os.path.splitext(f)
            name, mod = newname(file, ext, mod)
            conversions[f] = name + ext
    return mod

def flat_count(lst):
    """Count (recursively) the total number of non-list items in a list."""
    count = 0
    for item in lst:
        if isinstance(item, list):
            count += flat_count(item)
        else:
            count += 1
    return count

if __name__ == "__main__":

    import sys, os, shutil, subprocess
    from string import ascii_letters

    nfiles = flat_count(files)
    convert_files(files, map(str, range(1, nfiles)))

    dest = None
    tag = ''
    for i in xrange(1, 11):
        try:
            dest = 'apj%s' % tag
            os.mkdir(dest)
            break
        except:
            dest = None
            tag = '-%i' % i
            
    if not dest:
        assert 0, "Can't create a spare apj directory."

    has_tex = False

    for k,v in conversions.iteritems():
        print "Copying %s to %s/%s" % (k, dest, v)
        shutil.copy(k, '%s/%s' % (dest,v))
        file, ext = os.path.splitext(v)
        if v == '.tex': has_tex = True

    if has_tex:
        cmd = 'cd %s && cat ms.tex | sed "' % dest
        for k,v in conversions.iteritems():
            cmd += 's/{%s}/{%s}/g; ' % (k,v)
        cmd += '" > happyapj.tex && mv happyapj.tex ms.tex'
        subprocess.call(cmd, shell=True)

