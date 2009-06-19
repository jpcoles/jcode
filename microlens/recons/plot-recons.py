from __future__ import division, with_statement
from numpy import dtype, loadtxt, fromstring

fmt = dtype({'CNS NAME':            ('S15',   5)})
#            '# comp':              ('S1',   21),
#            'LHS':                 ('S4',   25),
#            'RA':                  ('S10',  32),
#            'DEC':                 ('S9',   43),
#            'position ref':        ('S1',   53),
#            'proper motion "/yr':  ('S5',   57),
#            'proper motion angle': ('S5',   63),
#            'pm ref':              ('S1',   69),
#            'parallax':            ('S7',   73),
#            'parallax err':        ('S7',   81),
#            'Spectral':            ('S4',   90),
#            'Spectral ref':        ('S1',  104),
#            'V':                   ('S6',  107),
#            'Mv':                  ('S5',  116),
#            'mass':                ('S4',  124),
#            'notes':               ('S17', 132),
#            'Common Name':         ('S36', 152)})

fmt = [['CNS NAME',           15,    5, None],
       ['# comp',              1,   21, None],
       ['LHS',                 4,   25, None],
       ['RA',                 10,   32, None],
       ['DEC',                 9,   43, None],
       ['position ref',        1,   53, None],
       ['proper motion "/yr',  5,   57, float],
       ['proper motion angle', 5,   63, float],
       ['pm ref',              1,   69, None],
       ['parallax',            7,   73, float],
       ['parallax err',        7,   81, float],
       ['Spectral',            4,   90, None],
       ['Spectral ref',        1,  104, None],
       ['V',                   6,  107, None],
       ['Mv',                  5,  116, None],
       ['mass',                4,  124, float],
       ['notes',              17,  132, None],
       ['Common Name',        36,  152, None]]

reconsSOA = []
reconsAOS = {}
for col,length,offs,func in fmt:
    reconsAOS[col] = []


with open('RECONS-100-nearest-systems.X', 'r') as f:
    for lineno,line in enumerate(f):
        line = line.rstrip()
        print line
        if not 9 <= lineno <= 257: continue
        if len(line) == 0: continue
    
        l = {}
        for col,length,offs,func in fmt:
            val = line[offs:offs+length].strip()
            print col, val
            if val and func: val = func(val)
            reconsAOS[col].append(val)
            l[col] = val
        reconsSOA.append(l)

