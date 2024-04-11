
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from os.path import join,exists
import os
from multiprocessing import Pool
import matplotlib.pyplot as plt
import json
import argparse


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument(
                '-i', '--pathToData', 
                help = 'File path to data table (must be a tab-delimited file).',
                required=True      
        ) 
        parser.add_argument(
                '-o', '--outdir', 
                help = 'Path to write all the outputs to.',
                required=True
        )
        parser.add_argument(
                '-c', '--label_combo', 
                help = 'comma delimited list of labels(celltypes/annotations) to compare',
                required=True
        )
        parser.add_argument(
                '-rf', '--roi_field', 
                help = 'The name of column which spcifies ROI. Default is sample_id',
                default="sample_id"  
        ),
        parser.add_argument(
                '-pf', '--position_fields', 
                help = 'The name of column which spcifies ROI. Default is sample_id',
                default="x,y" 
        ),
        parser.add_argument(
                '-lf', '--label_field', 
                help = 'The name of column which specifies cell type/annotation. Default is annotation',
                default="annotation",
        )
        parser.add_argument(
                '-pt', '--plot_type', 
                help = 'the type of plot pdf,png or svg. Default is png',
                default="png"
        )
        parser.add_argument(
                '-ri', '--rois', 
                help = 'A comma delimited list of rois to analyse. Default is all rois',   
        ),      
        parser.add_argument(
                '-rm', '--roi_metadata', 
                help = 'a tab delimited text file with info about each roi i.e. wdith , height, disease state',
        ),
        parser.add_argument(
                '-g', '--groups', 
                help = 'comma delimited list of states to group by , each state be a column in the roi metadata e.g. disease state/condition',
        )
        parser.add_argument(
                '-bm', '--bootstrap_method', 
                help = 'How the bootstrap should be done.Either random position (position) or label shuffling (shuffle) Default is position',
                default="position",
        )
        parser.add_argument(
                '-bn', '--bootstrap_n', 
                help = 'Number of bootstrap iterations. Default is 1000',
                default=1000,
                type=int
        )
        parser.add_argument(
                '-t', '--threads', 
                help = 'Number of threads to use. Default is 1',
                default=1,
                type=int
        )
        parser.add_argument(
                '-m', '--maxR', 
                help = 'Maximum radius to calculate NCF. Default is 150',
                default=150,
                type=int
        )
       
        parser.add_argument(
                '-st', '--step', 
                help = 'Step size for NCF. Default is 5',
                default=5,
                type=int
        )
        args = parser.parse_args()
        position_fields = args.position_fields.split(",")
        roi_sizes={}
        conditions={}
        if args.roi_metadata:
            df =pd.read_csv(args.roi_metadata,sep="\t")
            if "width" in df.columns and "height" in df.columns:
                roi_sizes= {x[0]:[x[1],x[2]] for x in zip(df[df.columns[0]],df["width"],df["height"])}
            if args.groups:
                gps  = args.groups.split(",")
                for g in gps:
                     conditions[g]= {x[0]:x[1] for x in zip(df[df.columns[0]],df[g])}
        rois = None
        if args.rois:
            rois = args.rois.split(",")

        run_ncf(args.pathToData,args.label_combo.split(","),args.outdir,
            rois=rois,
            roi_sizes=roi_sizes,
            plot_type=args.plot_type,
            position_fields=position_fields,
            threads=args.threads,
            roi_field=args.roi_field,
            label_field=args.label_field,
            maxR=args.maxR,
            step=args.step,
            bootstrap_n=args.bootstrap_n,
            bootstrap_method=args.bootstrap_method)
        for c in conditions:
            group_by_condition(args.outdir,conditions[c],c,args.plot_type)


def run_ncf(data,label_combo,outdir,
            rois=None,roi_sizes={},
            position_fields=["x","y"],
			plot_type="png",
            threads=1,roi_field="sample_id",label_field="annotation",
            maxR=150,step=5,bootstrap_n=1000,bootstrap_method="position"):
    df= pd.read_csv(data,sep="\t")
    #get all rois
    all_rois= df[roi_field].unique()

    #use all rois if a list is not supplied
    if not rois:
        rois=list(all_rois)
    if not exists(outdir):
        os.makedirs(outdir)
    bins = np.arange(0,maxR+step,step)
   
    #is the run halfway through?
    sum_file = join(outdir,"summary.json")
    if exists(sum_file):
        with open(sum_file) as o1:
            summary = json.loads(o1.read())
            maxR=summary["maxR"]
            step=summary["step"]
            bootstrap_n=summary["bootstrap_n"]
            bootstrap_method=summary["bootstrap_method"]
            label_combo=summary["combination"]
        
    else:
        summary = {
            "combination":label_combo,
            "maxR":maxR,
            "step":step,
            "bootstrap_n":bootstrap_n,
            "bootstrap_method":bootstrap_method,
            "numbers":{}
    }
    for roi in rois:
        points = np.array([[x[0],x[1]] for x in zip(df[position_fields[0]],df[position_fields[1]],df[roi_field]) if x[2]==roi])
        #get point co-ords plus labels for the roi
        labels= [x[0] for x in zip(df[label_field],df[roi_field]) if x[1]==roi]
        #get number of points for each label specified
        numbers = [labels.count(x) for x in label_combo]
        #can't run if any of the labels are missing
        if 0 in numbers:
            continue
        #allready been run
        if summary["numbers"].get(roi):
            continue
        circles,triplets = neighbourhoodCorrelationFunction(points,labels,label_combo,maxR=maxR)
        np.save(join(outdir,f"{roi}_circles.npy"),circles)
        np.save(join(outdir,f"{roi}_triplets.npy"),triplets)
        vals_observed, rs = np.histogram(np.array(circles).transpose()[2],bins=bins)
        circles = None
        triplets=None
        size  = roi_sizes.get(roi)
        if not size:
            size = points.max(0)
        arr =  [{"points":points,
            "roi_size":size,
            "labels":labels,
            "maxR":maxR,
            "bootstrap_method":bootstrap_method,
            "label_combo":label_combo,
            "bins":bins} for x in range(bootstrap_n)]
        with Pool(threads) as p:
            bcurve = p.map(bootstrap,arr)
        roi
        np.save(join(outdir,f"{roi}_bscurves.npy"),np.asarray(bcurve))
        write_file(outdir,roi,vals_observed,bcurve,rs)
        summary["numbers"][roi]=numbers
        with open(sum_file,"w") as o1:
            o1.write(json.dumps(summary,indent=4))

    draw_plots(outdir,plot_type)



def group_by_condition(indir,roi_to_condition,name,plot_type="png"):
    summary = json.loads(open(join(indir,"summary.json")).read())
    out_dir  = join(indir,name)
    if not exists(out_dir):
        os.makedirs(out_dir)
    conditions={}
    for k,v in roi_to_condition.items():
        i = conditions.get(v)
        if not i:
            i=[k]
        else:
            i.append(k)
        conditions[v]=i
    bins = np.arange(0,summary["maxR"]+summary["step"],summary["step"])
 
    for condition,rois in conditions.items():
        observed  = []
        bscurves = []
        for roi in rois:
            fl = join(indir,f"{roi}_circles.npy")
            if not exists(fl):
                print ("missing",fl)
                continue
            arr = np.load(fl).transpose()[2]
            if len(observed)==0:
                observed=arr
            else:
                observed=np.concatenate((observed,arr))
            f2 = join(indir,f"{roi}_bscurves.npy")
            arr= np.load(f2)
            if len(bscurves)==0:
                bscurves=arr
            else:
                bscurves=np.add(bscurves,arr)
        vals, rs  =np.histogram(observed,bins=bins)
        write_file(out_dir,condition,vals,bscurves,rs)
    draw_plots(out_dir,plot_type)


def write_file(outdir,roi,vals_observed,bcurve,rs):
    per5 = np.percentile(bcurve,5,0)
    per95= np.percentile(bcurve,95,0)
    average = np.mean(bcurve,axis=0)
    with open(join(outdir,f"{roi}_curves.txt"),"w") as o1:
        o1.write("\t"+"\t".join([str(x) for x in rs[1:]])+"\n")
        o1.write("observed\t"+"\t".join([str(x) for x in vals_observed])+"\n")
        o1.write("bootstrap_per5\t"+"\t".join([str(x) for x in per5])+"\n")
        o1.write("bootstrap_per95\t"+"\t".join([str(x) for x in per95])+"\n")
        o1.write("bootstrap\t"+"\t".join([str(x) for x in average])+"\n")

def draw_plots(dir,image_format="svg"):
    for curve_file in [join(dir,x) for x in os.listdir(dir) if "_curves.txt" in x]:
        name = curve_file.split("/")[-1].replace("_curves.txt","")
        data = pd.read_csv(curve_file,sep="\t")
        data.set_index("Unnamed: 0",inplace=True)
        data.columns= data.columns.astype(int)

        plt.rcParams['font.family'] = 'helvetica'
        #set font size 12pt for all text
        plt.rcParams['font.size'] = 35
        #figure size
        plt.figure(figsize=(12,9))
        #plot line thickness 1pt
        plt.rcParams['lines.linewidth'] = 4
        #line color black
        plt.rcParams['lines.color'] = 'k'
        #set y scale form 0 to 5
        #plt.ylim([0,6])
        #label x and y axis
        plt.xlabel('r ($\\mu$m)')
        plt.ylabel('g (r)')
        #draw gray line at y=1

        #plot the lines
        plt.plot(data.loc["observed"],color="#946BA9",lw=5)
        plt.plot(data.loc["bootstrap"],color="gray",linewidth=3,linestyle=(0,(5,3)))
        plt.fill_between(data.columns ,data.loc["bootstrap_per95"],data.loc["bootstrap_per5"],alpha=0.7,color="lightgray")

        #remove top and right axis
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        #increase thickness of bottom and left axis
        plt.gca().spines['bottom'].set_linewidth(3)
        plt.gca().spines['left'].set_linewidth(3)
        plt.locator_params(axis='x', nbins=7)
        #thick ticks
        plt.gca().tick_params(width=3)
        #ticks face inwards
        plt.gca().tick_params(direction='in')
        plt.savefig(join(dir,f"{name}_ncf_1.{image_format}"),dpi=300,format=image_format,bbox_inches="tight")
        plt.close()
        plt.figure(figsize=(12,9))
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.ylim([0,6])
        #increase thickness of bottom and left axis
        plt.gca().spines['bottom'].set_linewidth(3)
        plt.gca().spines['left'].set_linewidth(3)
        plt.plot(data.loc["observed"]/data.loc["bootstrap"],lw=5,color="#946BA9")
        plt.gca().axhline(1,color="lightgray",lw=2)
        plt.xlabel('$r$ ($\mu$m)')
        plt.ylabel('NCF$_{C_1 C_2 C_3}(r)$')
        plt.locator_params(axis='x', nbins=7)
        plt.savefig(join(dir,f"{name}_ncf_2.{image_format}"),dpi=300,format=image_format,bbox_inches="tight")
        plt.close()


def bootstrap(args):
    np.random.seed(random.randint(0,10000000))
    x = args["roi_size"][0]
    y = args["roi_size"][1]
    if args["bootstrap_method"]=="position":
        np1 = np.random.rand(len(args["points"]),2)
        np1[:,0]*=x
        np1 [:,1]*= y
        labels=args["labels"]
    else:
        np1= args["points"]
        labels = args["labels"].copy()
        random.shuffle(labels)
    circles,triplets = neighbourhoodCorrelationFunction(np1,labels,args["label_combo"],maxR=args["maxR"])
    radii = np.array(circles)
    
    if len(radii)==0:
        radii=np.array([])
    else:
        radii = np.array(radii).transpose()[2]
    h= np.histogram(radii,bins=args["bins"])[0]
    return h


# 
# Smallest enclosing circle - Library (Python)
# 
# Copyright (c) 2020 Project Nayuki
# https://www.nayuki.io/page/smallest-enclosing-circle
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING.txt and COPYING.LESSER.txt).
# If not, see <http://www.gnu.org/licenses/>.
# 
# This version from https://www.nayuki.io/res/smallest-enclosing-circle/smallestenclosingcircle.py dated 27/07/2022

import math, random

# Data conventions: A point is a pair of floats (x, y). A circle is a triple of floats (center x, center y, radius).

# Returns the smallest circle that encloses all the given points. Runs in expected O(n) time, randomized.
# Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
# Output: A triple of floats representing a circle.
# Note: If 0 points are given, None is returned. If 1 point is given, a circle of radius 0 is returned.
# 
# Initially: No boundary points known
def make_circle(points):
	# Convert to float and randomize order
	shuffled = [(float(x), float(y)) for (x, y) in points]
	random.shuffle(shuffled)
	
	# Progressively add points to circle or recompute circle
	c = None
	for (i, p) in enumerate(shuffled):
		if c is None or not is_in_circle(c, p):
			c = _make_circle_one_point(shuffled[ : i + 1], p)
	return c


# One boundary point known
def _make_circle_one_point(points, p):
	c = (p[0], p[1], 0.0)
	for (i, q) in enumerate(points):
		if not is_in_circle(c, q):
			if c[2] == 0.0:
				c = make_diameter(p, q)
			else:
				c = _make_circle_two_points(points[ : i + 1], p, q)
	return c


# Two boundary points known
def _make_circle_two_points(points, p, q):
	circ = make_diameter(p, q)
	left  = None
	right = None
	px, py = p
	qx, qy = q
	
	# For each point not in the two-point circle
	for r in points:
		if is_in_circle(circ, r):
			continue
		
		# Form a circumcircle and classify it on left or right side
		cross = _cross_product(px, py, qx, qy, r[0], r[1])
		c = make_circumcircle(p, q, r)
		if c is None:
			continue
		elif cross > 0.0 and (left is None or _cross_product(px, py, qx, qy, c[0], c[1]) > _cross_product(px, py, qx, qy, left[0], left[1])):
			left = c
		elif cross < 0.0 and (right is None or _cross_product(px, py, qx, qy, c[0], c[1]) < _cross_product(px, py, qx, qy, right[0], right[1])):
			right = c
	
	# Select which circle to return
	if left is None and right is None:
		return circ
	elif left is None:
		return right
	elif right is None:
		return left
	else:
		return left if (left[2] <= right[2]) else right


def make_diameter(a, b):
	cx = (a[0] + b[0]) / 2
	cy = (a[1] + b[1]) / 2
	r0 = math.hypot(cx - a[0], cy - a[1])
	r1 = math.hypot(cx - b[0], cy - b[1])
	return (cx, cy, max(r0, r1))


def make_circumcircle(a, b, c):
	# Mathematical algorithm from Wikipedia: Circumscribed circle
	ox = (min(a[0], b[0], c[0]) + max(a[0], b[0], c[0])) / 2
	oy = (min(a[1], b[1], c[1]) + max(a[1], b[1], c[1])) / 2
	ax = a[0] - ox;  ay = a[1] - oy
	bx = b[0] - ox;  by = b[1] - oy
	cx = c[0] - ox;  cy = c[1] - oy
	d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
	if d == 0.0:
		return None
	x = ox + ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
	y = oy + ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d
	ra = math.hypot(x - a[0], y - a[1])
	rb = math.hypot(x - b[0], y - b[1])
	rc = math.hypot(x - c[0], y - c[1])
	return (x, y, max(ra, rb, rc))


_MULTIPLICATIVE_EPSILON = 1 + 1e-14

def is_in_circle(c, p):
	return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * _MULTIPLICATIVE_EPSILON


# Returns twice the signed area of the triangle defined by (x0, y0), (x1, y1), (x2, y2).
def _cross_product(x0, y0, x1, y1, x2, y2):
	return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)


def neighbourhoodCorrelationFunction(points,labels,categoriesToPlot,maxR=0.5):
    nCategories = len(categoriesToPlot)
    assert(nCategories > 1)
    pointsToCompare = []
    for category in categoriesToPlot:
        p = np.array([[float(x[0][0]),float(x[0][1])] for x in zip(points,labels) if x[1] == category])
        pointsToCompare.append(p)
    # Prefiltering
    AB = cdist(pointsToCompare[0],pointsToCompare[1]) < maxR*2
    AC = cdist(pointsToCompare[0],pointsToCompare[2]) < maxR*2
    BC = cdist(pointsToCompare[1],pointsToCompare[2]) < maxR*2
    ab_candidates = np.where(AB)
    allCandidates = []
    for i in range(np.shape(ab_candidates)[1]):
        target = [ab_candidates[0][i], ab_candidates[1][i]]
        # Check this against BC and AC to see if we have a candidate
        successes = np.where(AC[target[0],:] & BC[target[1],:])[0]
        if len(successes) > 0:
            allCandidates.extend([[target[0],target[1],v] for v in successes])
    allCandidates = np.asarray(allCandidates)
    # allCandidates is now a nCandidates x 3 array, where nCandidates is the number of triplets which are pairwise within a distance of 2*maxR
    # For each triplet, get the points and calculate the smallest enclosing circle
    nCandidates = np.shape(allCandidates)[0]
    circles = []
    triplets = []
    for i in range(nCandidates):
        (center_x, center_y, radius) = make_circle([pointsToCompare[v][allCandidates[i,v]] for v in range(nCategories)])
        if radius < maxR:
            circles.append([center_x, center_y, radius])
            triplets.append([pointsToCompare[v][allCandidates[i,v]] for v in range(nCategories)])

    return circles, triplets





'''
run_ncf("/ceph/project/IPF_UKRMP/shared/annotations_2/cells.txt",
        ["ABI_2_DC_ADJ (10)","CD206hi_mac (1)","Prolif_mono (18)"],
        "mytest",
        bootstrap_n=5,
        threads=2,
		plot_type="png",
        label_field="annotation_2")
'''
#input_file  = "/ceph/project/IPF_UKRMP/shared/annotations_2/disease_states.txt"
#df = pd.read_csv(input_file,sep="\t")

#di = {x[0]:x[1] for x in zip(df["sample_id"],df["disease_state_v3"])}


#collapse_to_condition("mytest",di,"disease_state",plot_type="png")

if __name__ == "__main__":
    main()