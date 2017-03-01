#!/usr/bin/env python

from __future__ import division

import numpy as np
from scipy.interpolate import griddata
import argparse
import os
import h5py


################################################################
# function to compute sliced posterior
def Do_Slicing(data, delta):
    condition=[]
    for v in VARS:
        if argsd[v] is not None:
            if argsd['frac_'+v]:
                # fractional width
                condition.append(np.abs((argsd[v]-data[v])/argsd[v])<delta*argsd['mult_'+v])
            else:
                # absolute width
                condition.append(np.abs(argsd[v]-data[v])<delta*argsd['mult_'+v])
    out=data[np.where(np.all(condition,axis=0))]
    return out
################################################################


def PrintOnlyExisting(data, members):
    """Helper function for printing only those members that exist in the array 'data'"""
    q=[]
    for m in members:
        if m in data.dtype.names: q.append(m)
    print "# {}={}".format(q,data[q][0])


################################################################
def Analyse(data, delta, argsd, h5_file):
    """Given an imported posterior array 'data', parameters stored in
    argsd, and a width 'delta', analyse the posterior array and 
    print results to screen"""
   
################################################################
    # do slicing
    if delta is not None:
        out=Do_Slicing(data, delta)
    else:
        print "# Determine delta by bisection.  Goal Nremaining={}".format(Nremaining)
        # need to bisect
        delta_min=0.
        delta_max=1.
        out_min=Do_Slicing(data, delta_min)
        if len(out_min) > Nremaining:
            raise Exception('logic error in algorithm')
        out_max=Do_Slicing(data, delta_max)
        while len(out_max)<Nremaining:
            delta_max=delta_max*8.
            out_max=Do_Slicing(data, delta_max)
            if delta>1000.:
                raise Exception('failed to bound Nremaining samples')
        while(len(out_max) != len(out_min)):
            print "# delta-range=[{:.6f}, {:.6f}], Nsamples=[{}, {}]".format(delta_min, delta_max, len(out_min), len(out_max))
            # bifurcate
            delta=(delta_min+delta_max)/2.
            out=Do_Slicing(data, delta)
            if len(out)==Nremaining: break
            if len(out)<Nremaining:
                delta_min=delta
                out_min=out
            else:
                delta_max=delta
                out_max=out

    ################
    # generate nice output
    s=""
    for v in VARS:
        if argsd[v] is None: continue
        range=delta*argsd['mult_'+v]
        if argsd['frac_'+v]: range=range*argsd[v]
        s=s+", {}={}+-{:.4f}".format(v, argsd[v],range)
    s=s[2:] # eliminate leading comma


    if len(out) == 0:
        print "No sample point found within "+s
        return

    print "# Highest "+ranking+" point found within "+s
    PrintOnlyExisting(data, ['q', 'a1z', 'a2z', 'tilt1', 'tilt2'])
    print "# ['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']={}".format(out[['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']][0])
    print "# ['time','h1_end_time','l1_end_time','h1l1_delay']={}".format(out[['time','h1_end_time','l1_end_time','h1l1_delay']][0])
    
    if args.likelihood:
        maxLpost=(data['logl'])[0]
        maxLslice=(out['logl'])[0]
    else:
        maxLpost=(data['logprior']+data['logl'])[0]
        maxLslice=(out['logprior']+out['logl'])[0]
        print "# overall highest log"+ranking+" ={}".format(maxLpost)
        print "# highest log"+ranking+" in slice={}".format(maxLslice)
        print "# difference in log"+ranking+"={}".format(maxLpost-maxLslice)


    if h5_file!="":
        argsd['mtotal']=out['mtotal'][0]
        print "{: <25}  {:2.5f}  {:1.4f}  {mtotal:4.2f}  {q:.5f}  {a1z: .4f}  {a2z: .4f}  {tilt1:.3f}  {tilt2:.3f}".format(h5_file,maxLpost-maxLslice, delta, **argsd)

    if argsd['output']:
        print "# Writing "+str(len(out))+" points to file "+argsd['output']
        np.savetxt(argsd['output'],out,header='\t'.join(out.dtype.names))




################################################################
################################################################
#### MAIN
################################################################
################################################################


parser = argparse.ArgumentParser(description=

    """Slices and interpolate a LIGO CBC-PE posterior file to give the
parameters (total mass and extrinsic parameters) fitting best specific
mass-ratio, a1z, a2z, and (if desired) tilt1 and tilt2 values.
""")


# list of slice-able variables.  This list
# can be extended as desired
VARS=['q', 'a1z', 'a2z', 'tilt1', 'tilt2']



parser.add_argument('-p', '--posterior', type=str, help='Posterior file (ASCII table with one-line header).')

parser.add_argument('-o', '--output', type=str, nargs='?', help='Output file.',
                    default=None)

parser.add_argument('-l', '--likelihood', action='store_true', help='Use the likelihood instead of the posterior to rank points.')

parser.add_argument('-c', '--credible-interval', type=float, help='Credible interval to use as subset of the posterior (defaul: %(default)s).',
                    default=0.9)

parser.add_argument('--params_from_h5_files', type=str, default='', nargs='*',
                    help="Read target parameters for q, a1z, a2z, tilt1, tilt2 from the metadata of this LVC-NR waveform file.")

group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--delta', type=float, default=None,
                   help='Allowed +/- range around the given mass-ratio, spin1, spin2 values.')
group.add_argument('--Nremaining', type=int, default=None,
                    help='Choose a delta such that Nremaining posterior samples remain in sliced posterior.')

for v in VARS:
    parser.add_argument('-{}'.format(v), type=float, default=None,
                        help="specify a value for {} at which to slice".format(v))
for v in VARS:
    parser.add_argument('-frac_{}'.format(v), default=False, action='store_true',
                        help='enforce a fractional width for variable {}'.format(v))
for v in VARS:
    parser.add_argument('-mult_{}'.format(v), default=1., type=float,
                        help="width multiplier of this variable in the posterior slice (default: %(default)s)")


args = parser.parse_args()
delta=args.delta
Nremaining=args.Nremaining
argsd=vars(args) # turn args into dictionary for use in Do_Slicing




# Import data
org_data=np.genfromtxt(args.posterior, names=True)

# Sort the data in order of increasing logposterior=logprior+loglikelihood
if args.likelihood:
    ranking='likelihood'
    sorted_data=org_data[np.argsort(org_data['logl'])][::-1]
else:
    ranking='posterior'
    sorted_data=org_data[np.argsort(org_data['logprior']+org_data['logl'])][::-1]
# Take the top (credible-interval)
cred_data=sorted_data[0:int(args.credible_interval*len(sorted_data))]

data=cred_data

#print args.params_from_h5_files
#os.sys.exit()
if len(args.params_from_h5_files)==0:
    Analyse(data, delta, argsd, "")

else:

    # LOOP OVER ALL GIVEN h5-FILES
    for h5_file in args.params_from_h5_files:
        F=h5py.File(h5_file, 'r')
        A=F.attrs
        q=A['mass2']/A['mass1']
        vec1=[A['spin1x'], A['spin1y'], A['spin1z']]
        vec2=[A['spin2x'], A['spin2y'], A['spin2z']]
        a1=(vec1[0]**2+vec1[1]**2+vec1[2]**2)**0.5
        a2=(vec2[0]**2+vec2[1]**2+vec2[2]**2)**0.5
        a1z=A['spin1z']
        a2z=A['spin2z']
        tilt1=np.arctan2((vec1[0]**2+vec1[1]**2)**0.5, a1z)
        tilt2=np.arctan2((vec2[0]**2+vec2[1]**2)**0.5, a2z)
        print "# Read from {}:  q={}".format(h5_file, q)
        print '#  vec{{chi}}_1={}, a1={}, tilt1={}'.format(vec1, a1, tilt1)
        print '#  vec{{chi}}_2={}, a2={}, tilt2={}'.format(vec2, a2, tilt2)
        argsd['q']=q
        argsd['a1z']=a1z
        argsd['tilt1']=tilt1
        argsd['a2z']=a2z
        argsd['tilt2']=tilt2
        Analyse(data, delta, argsd, h5_file)



#pts=np.array([[x,y,z] for x,y,z in cred_data[['q','a1z','a2z']]])
#mtot=griddata(pts, cred_data['mtotal'], (q, a1z, a2z), method='linear')
#print "Interpolated total mass: "+str(mtot)
