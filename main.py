# GroupID-22 (22114030_22114057_22114073)
# 24th Sept, 2025
#main.py File calling all the functions and the GUI

from numpy.random import randint
from easygui import multenterbox
from numpy import lexsort,asarray,append
from matplotlib.pyplot import figure,show
import matplotlib.pyplot as plt
from ThreeColoring import *
from DCEL import *
from MonotonePartitioning import *
from Triangulation import *
DEBUG = False
from EdgeGuardSolver import mobile_edge_guards, plot_mobile_guards
import copy
import math

msg = "Enter the value of 'n'"
title = "Polygon with 'n' sides"
fieldNames = ["  n"]
fieldValues = multenterbox(msg,title,fieldNames)
while 1:
    if fieldValues == None: break
    errmsg = ""
    try:
        int(fieldValues[0])
    except:
        errmsg = errmsg + ('%s is a required to be a number.\n\n' % fieldNames[0])
    if fieldValues[0].strip() == "":
        errmsg = errmsg + ('%s is a required field.\n\n' % fieldNames[0])
    if errmsg == "": break
    fieldValues = multenterbox(errmsg, title, fieldNames, fieldValues)
    


n = int(fieldValues[0])
from numpy.random import randint
from numpy import lexsort,asarray,append
coords = randint(0,90000, size=(2,n))
x = coords[0] 
y = coords[1]
ind = lexsort((y,x))
coords = [(x[i],y[i]) for i in ind] 
x = asarray([c[0] for c in coords])
y = asarray([c[1] for c in coords])
pivot = coords[0]

y_diff_pivot = y-pivot[1]
x_diff_pivot = x-pivot[0]
tan = (y_diff_pivot[1:]+0.0)/x_diff_pivot[1:]

pairs = zip(tan,coords[1:])
# pairs = sorted(pairs, key = lambda (x,y): x)
pairs = sorted(pairs, key=lambda pair: pair[0])

coords = asarray([pivot])
coords = append(coords, asarray([c for _, c in pairs]), axis=0)
coords = append(coords, asarray([pivot]), axis=0 )

p=[(c[1],c[0]) for c in coords][::-1][1:]
d = buildSimplePolygon(p)
print("################################ DCEL POLYGON FORMED ###################################")
map_points ={x.coords:x for x in d.getVertices()}
for i in range(1,len(p)-1):
    map_points[p[i]].next = map_points[p[i+1]]
    map_points[p[i+1]].prev = map_points[p[i]]
map_points[p[0]].prev = map_points[p[-1]]
map_points[p[0]].next = map_points[p[1]]
map_points[p[1]].prev = map_points[p[0]]
map_points[p[-1]].next = map_points[p[0]]

d1 = d
def newline(p, q):
    X = np.linspace(p[0], q[0], endpoint=True)
    Y = np.linspace(p[1], q[1], endpoint=True)
    return X, Y

def nowDraw(toDraw):
    for x in toDraw:
        plt.plot(x[0], x[1], x[2])
    plt.show()

frames = [] 

def push_frame(toDraw_list, title=""):
    frames.append({"toDraw": copy.deepcopy(toDraw_list), "title": title})

def nowDrawMulti(frames_list=None, ncols=3, figsize_per_panel=(4,3)):
    if frames_list is None:
        fl = frames
    else:
        fl = frames_list

    if not fl:
        print("no frames to draw")
        return

    m = len(fl)
    cols = ncols
    rows = math.ceil(m / cols)
    figsize = (figsize_per_panel[0] * cols, figsize_per_panel[1] * rows)
    fig, axes = plt.subplots(rows, cols, figsize=figsize, squeeze=False)
    fig.tight_layout(pad=3.0)

    for idx, frame in enumerate(fl):
        r = idx // cols
        c = idx % cols
        ax = axes[r][c]

        for item in frame["toDraw"]:
            X, Y, style = item
   
            try:
                ax.plot(X, Y, style)
            except Exception:
                try:
                    ax.plot([X], [Y], style)
                except Exception:
                    pass
        ax.set_title(f"Step {idx+1}: {frame.get('title','')}")
        ax.set_aspect('equal', 'box')
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # turn off empty subplots if any
    for j in range(m, rows * cols):
        r = j // cols
        c = j % cols
        axes[r][c].axis('off')

    plt.show()

toDraw = []
for e in d.getEdges():
    p1,q1 = list(e.origin.coords),list(e.getTwin().origin.coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'r'])

push_frame(toDraw, title="Original polygon")
DEBUG = False
ret = getTrapEdges(d)
for r in ret:
    X,Y = newline(r.left,r.right)
    toDraw.append([X,Y,'b'])
    toDraw.append([r.pivot.coords[0],r.pivot.coords[1],'go'])

push_frame(toDraw, title="Trapezoidal edges and pivots")
DEBUG = False
diagnls = monotonePartitioningDgnls(d)
for dg in  diagnls:
    X,Y = newline(dg[0].coords,dg[1].coords)
    toDraw.append([X,Y,'g'])

push_frame(toDraw, title="Partition diagonals (monotone split)")
listOfMonos = insertDgnls(d,[(x[0].coords,x[1].coords) for x in diagnls])
import os

listOfMonos = insertDgnls(d, [(x[0].coords, x[1].coords) for x in diagnls])

toDraw = []
for m in listOfMonos:
    for e in m.getEdges():
        p1,q1 = list(e.origin.coords),list(e.getTwin().origin.coords)
        X,Y = newline(p1,q1)
        toDraw.append([X,Y,''])
        toDraw.append([e.origin.coords[0],e.origin.coords[1],'o'])
push_frame(toDraw, title="Monotone polygons")

DEBUG = True
toDraw = []
listOfTriangles = []
tmp = -1
# print len(listOfMonos)
for mono in listOfMonos:
    diagnls = triangulateMonotonePolygon(mono)
    vv = [(x[0].coords,x[1].coords) for x in diagnls]
    listOfTriangles += insertDgnls(mono,vv)
    tmp+=len(diagnls)+1

listOfTriangles = [[t.getFaces()[1].getOuterBoundary()[0].origin,
                    t.getFaces()[1].getOuterBoundary()[1].origin,
                    t.getFaces()[1].getOuterBoundary()[2].origin
                   ] for t in listOfTriangles]

for t in listOfTriangles:
    p1,q1 = list(t[0].coords),list(t[1].coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'r'])
    p1,q1 = list(t[1].coords),list(t[2].coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'b'])
    p1,q1 = list(t[2].coords),list(t[0].coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'g'])


push_frame(toDraw, title="Triangulated polygons")

colorizer = Colorizer(d, listOfTriangles)
x = colorizer.colorize()

# ---------------- 3-coloring visualization ----------------
colors = colorizer.colors  
counts = {0:0, 1:0, 2:0}

toDraw = []
# draw polygon edges
for e in d1.getEdges():
    p1,q1 = list(e.origin.coords), list(e.getTwin().origin.coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'k-'])

# draw vertices with color markers
for coords, c in colors.items():
    if c == 0:
        toDraw.append([[coords[0]], [coords[1]], 'ro'])  # red
    elif c == 1:
        toDraw.append([[coords[0]], [coords[1]], 'go'])  # green
    elif c == 2:
        toDraw.append([[coords[0]], [coords[1]], 'bo'])  # blue
    counts[c] += 1

push_frame(
    toDraw,
    title=f"3-coloring (red={counts[0]}, green={counts[1]}, blue={counts[2]})"
)

toDraw = []
for e in d1.getEdges():
    p1,q1 = list(e.origin.coords),list(e.getTwin().origin.coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'r'])
for g in x[0]:
    toDraw.append([g.coords[0],g.coords[1],'bo'])
push_frame(toDraw, title=f"Vertex guards (count {x[1]})")

# # ---------------- Mobile edge guard solution ----------------
# guards, reduced_poly = mobile_edge_guards(d, listOfTriangles, debug=True)
# print(f"Mobile edge guards required: {len(guards)}")
# plot_mobile_guards(d, guards, nowDraw, newline)

guards, reduced_poly = mobile_edge_guards(d, listOfTriangles, debug=True)
print(f"Mobile edge guards required: {len(guards)}")

toDraw = []
# draw polygon edges
for e in d.getEdges():
    p1,q1 = list(e.origin.coords),list(e.getTwin().origin.coords)
    X,Y = newline(p1,q1)
    toDraw.append([X,Y,'k-'])
# draw guards
for g in guards:
    p1,p2 = g
    X,Y = newline(list(p1), list(p2))
    toDraw.append([X,Y,'y-'])
    mx, my = (p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0
    toDraw.append([[mx],[my],'ro'])

push_frame(toDraw, title=f"Mobile edge guards (count={len(guards)})")


nowDrawMulti(frames, ncols=3, figsize_per_panel=(4,3))
