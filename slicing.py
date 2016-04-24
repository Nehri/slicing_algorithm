# -*- coding: utf-8 -*-

import math
import sys
import string

#printerbed size: 500mm side cube

class Point:
    def __init__(self, x_, y_, z_):
        self.x = x_
        self.y = y_
        self.z = z_

    def dotProduct(self, p):
        return self.x*p.x + self.y*p.y + self.z*p.z

    def normalize(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def toString(self):
        return "("+str(self.x)+","+str(self.y)+","+str(self.z)+")"


class Line:
    def __init__(self, p0_, p1_):
        self.p0 = p0_
        self.p1 = p1_

def lineEqual(L1,L2):
    if (((L1.p0.x == L2.p0.x) and (L1.p0.y == L2.p0.y) and
        (L1.p1.x == L2.p1.x) and (L1.p1.y == L2.p1.y))
    or ((L1.p0.x == L2.p1.x) and (L1.p0.y == L2.p1.y) and
        (L1.p1.x == L2.p0.x) and (L1.p1.y == L2.p0.y))):
        return True
    else:
        return False

class Triangle:
    def __init__(self, p0_, p1_, p2_, norm_):
        self.p1 = p0_
        self.p2 = p1_
        self.p3 = p2_
        self.norm = norm_

class Slice:
    def __init__(self, perimeter_, isSurface_):
        self.perimeter = perimeter_
        self.isSurface = isSurface_
        self.support = list()
        self.infill = list()

# given an stl file of standard format,
# returns a list of the triangles
def fileToTriangles(filename):
    with open(filename, 'r') as f:
        next(f)
        counter = 0
        triangles = list()
        triangle = Triangle(p0_=None, p1_=None, p2_=None, norm_=None)
        for line in f:
            l_ = line.split(" ")
            l = [value for value in l_ if value != '']
            if counter in [1,5,6]:
                if counter == 6:
                    counter = 0
                else:
                    counter += 1
                continue
            elif counter == 0:
                if l[0] == 'endsolid':
                    break
                triangle.norm=Point(float(l[2]), float(l[3]), float(l[4][:-1]))
                counter += 1
            elif counter == 2:
                triangle.p0 = Point(float(l[1]), float(l[2]), float(l[3]))
                counter+=1
            elif counter == 3:
                triangle.p1 = Point(float(l[1]), float(l[2]), float(l[3]))
                counter+=1
            elif counter == 4:
                triangle.p2 = Point(float(l[1]), float(l[2]), float(l[3]))
                triangles.append(triangle)
                counter += 1
        return triangles
         
# given a line segment and a plane,
# returns the point at which those two
# planes intersect. Will return the first point
# in the line if the whole line is in the plane
def intersectSlice(line, plane):
    if line.p0.z == line.p1.z and line.p1.z == plane:
        return line.p0
    elif line.p0.z == line.p1.z:
        return None
    else:
        slope = Point(x_=line.p1.x-line.p0.x, y_=line.p1.y-line.p0.y, z_=line.p1.z-line.p0.z)
        t = float(plane-line.p0.z)/float(slope.z)

        if t >= 0 and t <= 1:
            testZ = line.p0.z+t*slope.z
            if testZ <= max(line.p0.z, line.p1.z) and testZ >= min(line.p0.z, line.p1.z):
                testP = Point(x_=line.p0.x+t*slope.x, y_=line.p0.y+t*slope.y, z_=line.p0.z+t*slope.z)
                
            else: 
                return None
        else:
            return None
         
# given two lines on the same z plane,
# returns the point at which they intersect,
# or None if there is no intersection
def intersection(L1,L2):

    #make sure all lines are on the same z plane
    assert (L1.p0.z == L1.p1.z)
    assert (L2.p0.z == L2.p1.z)
    assert (L1.p0.z == L2.p0.z)

    x1 = L1.p0.x
    y1 = L1.p0.y
    x2 = L1.p1.x
    y2 = L1.p1.y
    x3 = L2.p0.x
    y3 = L2.p0.y
    x4 = L2.p1.x
    y4 = L2.p1.y

    xnum = (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)
    xden = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    ynum = (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)
    yden = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)

    try:
        intersect = Point(xnum/xden,ynum/yden,L1.p0.z)
        return intersect
    except:
        return None

# given a list of triangles in 3D space,
# returns a tuple of the highest and lowest Z values
def findBoundaries(triangles):
    bottomZ = 500
    topZ = -500

    for triangle in triangles:
        maximum = max(triangle.p0.z, triangle.p1.z, triangle.p2.z)
        minimum = min(triangle.p0.z, triangle.p1.z, triangle.p2.z)

        if maximum > topZ:
            topZ = maximum
        if minimum < bottomZ:
            bottomZ = minimum

    return (bottomZ, topZ)

# given a list of triangles and a thickness per layer,
# computes the total number of layers and checks which
# line segments need to be drawn in each layer,
# returns the list of slices with the list of line segments
# to draw per slice, as a tuple with a bool for if the slice 
# is a bottom or top
def separateSlices(triangles, layerThickness):
    bounds = findBoundaries(triangles)
    numSlices = int((bounds[1]-bounds[0])/layerThickness)
    slices = [bounds[0]+z*layerThickness for z in range(0, numSlices+1)]
    segments = list()
    for s in slices:
        currentSegment = list()
        currentSegmentSurface = False
        for triangle in triangles:
            point1 = intersectSlice(Line(p0_=triangle.p0, p1_=triangle.p1), s)
            point2 = intersectSlice(Line(p0_=triangle.p1, p1_=triangle.p2), s)
            point3 = intersectSlice(Line(p0_=triangle.p2, p1_=triangle.p0), s)

            points = list(set([point1, point2, point3]))

            for point in points:
                if point == None:
                    points.remove(None)
                    break
            if s <= (bounds[0]+layerThickness) or s >= (bounds[1]-layerThickness):
                currentSegmentSurface = True

            if len(points) == 1:
                currentSegment.append(Line(points[0], points[0]))
            elif len(points) == 2:
                currentSegment.append(Line(points[0], points[1]))
            elif len(points) == 3:
                segment1 = Line(points[0], points[1])
                segment2 = Line(points[1], points[2])
                segment3 = Line(points[2], points[0])
                currentSegmentSurface = True

                currentSegment.append(segment1)
                currentSegment.append(segment2)
                currentSegment.append(segment3)
         
        segments.append(Slice(perimeter_=currentSegment,isSurface_=currentSegmentSurface))

    return segments

# given two lines on the same z plane, 
# returns the point at which they intersect,
# or None if there is no intersection
def intersection(L1,L2):

    #make sure all lines are on the same z plane
    #assert (math.isclose(L1.p0.z, L1.p1.z, abs_tol=0.0001))
    #assert (L2.p0.z == L2.p1.z)
    #assert (L1.p0.z == L2.p0.z)

    x1 = L1.p0.x
    y1 = L1.p0.y
    x2 = L1.p1.x
    y2 = L1.p1.y
    x3 = L2.p0.x
    y3 = L2.p0.y
    x4 = L2.p1.x
    y4 = L2.p1.y

    xnum = (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)
    xden = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    ynum = (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)
    yden = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)

    try:
        intersect = Point(xnum/xden,ynum/yden,L1.p0.z) 

        if ((intersect.x >= min(x1,x2)) and (intersect.x <= max(x1,x2)) and
            (intersect.y >= min(y1,y2)) and (intersect.y <= max(y1,y2)) and
            (intersect.x >= min(x3,x4)) and (intersect.x <= max(x3,x4)) and
            (intersect.y >= min(y3,y4)) and (intersect.y <= max(y3,y4))):
            return intersect
        else:
            return None
    except:
        return None

#given a list of lines that make a manifold perimeter on a slice,
#and a percentage of space that should be infill,
#returns a list of infill lines (grid pattern) for that slice
#assumes print bed area is a square
def infill(perimeter,percent):

    assert (percent >= 0)
    assert (percent <= 1)

    if (len(perimeter) == 0):
        return []
    Z = perimeter[0].p0.z #should be the same across all lines

    #printer specific constants, should be suplied as args
    bedWidth = 500.0#mm
    extrudeWidth = 1.0#mm

    linesPerSide = round((bedWidth*percent)/(extrudeWidth*2))
    gap = bedWidth/linesPerSide
    infill = []

    #horizontal lines
    for y in range(int(linesPerSide)):
        
        #start with full line
        fullLine = Line(Point(0,y*gap,Z),Point(bedWidth,y*gap,Z))
        inters = []

        #find intersections
        for line in perimeter:
            i = intersection(line,fullLine)
            if (i != None):
                inters.append(i)

        #sort by x to get matching pairs for internal lines
        inters.sort(key=lambda point: point.x)
        assert(len(inters)%2 == 0) #if not even, then perimeter was not manifold
        for i in range(int(round(len(inters)/2))):
            overlap = False;
            newLine = Line(inters[i*2],inters[i*2+1])
            for l in perimeter:
                if lineEqual(l,newLine):
                    overlap = True;
            if not overlap:
                infill.append(newLine)

    #vertical lines
    for x in range(int(linesPerSide)):
        
        #start with full line
        fullLine = Line(Point(x_=x*gap,y_=0,z_=Z),Point(x_=x*gap,y_=bedWidth,z_=Z))
        inters = []

        #find intersections
        for line in perimeter:
            i = intersection(line,fullLine)
            if (i != None):
                inters.append(i)

        #sort by y to get matching pairs for internal lines
        inters.sort(key=lambda point: point.y)
        assert(len(inters)%2 == 0) #if not even, then perimeter was not manifold
        for i in range(int(round(len(inters)/2))):
            overlap = False;
            newLine = Line(inters[i*2],inters[i*2+1])
            for l in perimeter:
                if lineEqual(l,newLine):
                    overlap = True;
            if not overlap:
                infill.append(newLine)

    # for l in infill:
        # print("("+str(l.p0.x)+","+str(l.p0.y)+"),("+str(l.p1.x)+","+str(l.p1.y)+"), "+str(l.p0.z)+","+str(l.p1.z))

    return infill


# given a list of slices with the list of line segments
# to draw per slice, as a tuple with if the slice is a
# bottom or top, and a filename,
# write the G code to the given file 
def writeGcode(slices,filename):

    extrudeRate = 0.05
    f = open(filename[:-3] + "gcode",'w')

    #preamble
    f.write(";Start GCode\n")
    f.write("M109 S210.000000\n")
    f.write("G28 X0 Y0 Z0\n")
    f.write("G92 E0\n")
    f.write("G29\n")

    layer = 1;
    for s in slices:

        f.write(";Layer "+str(layer)+" of "+str(len(slices))+"\n")

        #fan
        if (layer == 2):
            f.write("M106 S127\n")
        if (layer == 3):
            f.write("M106 S255\n")
        
        f.write(";perimeter\n")
        for l in s.perimeter:
            #move to start of line
            f.write("G0 F2700 X"+str(l.p0.x)+" Y"+str(l.p0.y)+"Z "+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            f.write("G1 F900 X"+str(l.p1.x)+" Y"+str(l.p1.y)+"E "+str(dist*extrudeRate)+"\n")

        if len(s.support) > 0:
            f.write(";support\n")
        for l in s.support:
            #move to start of line
            f.write("G0 F2700 X"+str(l.p0.x)+" Y"+str(l.p0.y)+"Z "+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            f.write("G1 F900 X"+str(l.p1.x)+" Y"+str(l.p1.y)+"E "+str(dist*extrudeRate)+"\n")

        if len(s.infill) > 0:
            f.write(";infill\n")
        for l in s.infill:
            #move to start of line
            f.write("G0 F2700 X"+str(l.p0.x)+" Y"+str(l.p0.y)+"Z "+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            f.write("G1 F900 X"+str(l.p1.x)+" Y"+str(l.p1.y)+"E "+str(dist*extrudeRate)+"\n")

        
        layer+=1


    #postamble
    f.write(";End GCode\n")
    f.write("M104 S0\n")
    f.write("M140 S0\n")
    f.write("G91\n")
    f.write("G1 E-1 F300\n")
    f.write("G1 Z+0.5 E-5 X-20 Y-20 F2v700\n")
    f.write("G28 X0 Y0\n")
    f.write("M84\n")
    f.write("G90\n")

def main():
    filename = sys.argv[1]
    layerThickness = float(sys.argv[2])
    supportPercent = float(sys.argv[3])
    triangles = fileToTriangles(filename)

    slices = separateSlices(triangles, layerThickness)

    for s in slices:
        s.infill = infill(s.perimeter, supportPercent)

    writeGcode(slices,filename)
    

if __name__ == "__main__":
    main()
# pseudocode for computing brim of a single convex polyhedron base
# takes a listof(line segments) which are the base (bottom layer),
# a number of outlines, and an initial offset
# will also generate a skirt if offset is greater than 1 and number of outlines = 1
'''
def brim(base, outlines, offset):
    # remove interior lines in manifold
    for every line segment L in base
        if L is a duplicate, remove all L in base
    # compute centroid of polygonal manifold
    Cx = 1,Cy = 1, A = 1
    for every L in base:
        Cx = Cx * (L.p0.x + L.p1.x)
        Cy = Cy * (L.p0.y + L.p1.y)
        A = A * (L.p0.x * L.p1.y - L.p1.x * L.p0.y)
    A = A / 2
    Cx = Cx / (6 * A)
    Cy = Cy / (6 * A)

    # generate outlines
    brimlines = None
    for i from 1 to outlines:
        for L in base:
            L’ = L
            if L.p0.x > Cx:
                L’.p0.x += offset + extrudeWidth * i
            else:
                L’.p0.x -= offset + extrudeWidth * i
    if L.p0.y > Cy:
        L’.p0.y += offset + extrudeWidth * i
    else:
        L’.p0.y -= offset + extrudeWidth * i
    if L.p1.x > Cx:
        L’.p1.x += offset + extrudeWidth * i
    else:
        L’.p1.x -= offset + extrudeWidth * i
    if L.p1.y > Cy:
        L’.p1.y += offset + extrudeWidth * i
    else:
        L’p1.y -= offset + extrudeWidth * i
    brimlines = brimlines + {L’}
        return brimlines
'''


# pseudocode for computing a basic rectangular raft
# takes number of layers, an offset (how far raft extends from object), and an infill percentage
# and a listof(listof(line segments)) representing the object
# NOTE: raft goes UNDERNEATH object, all object / infill / support layers must be shifted up by
# height of raft (the brim, if used, is the same layer as raft and does not need to be elevated)
'''
def raft(layers, layernum, offset, infill):
	# compute enclosing rectangle on object
	x = 0, y = 0, _x = 0, _y = 0
	for L in layers
		for s in L
			if s.p0.x > x or s.p1.x > x:
				x = s.p0.x if s.p0.x > s.p1.x else s.p1.x
			if s.p0.x < _x or s.p1.x < _x:
				_x = s.p0.x if s.p0.x < s.p1.x else s.p1.x
			if s.p0.y > y or s.p1.y > y:
				y = s.p0.y if s.p0.y > s.p1.y else s.p1.y
			if s.p0.y < _y or s.p1.y < _y:
				_y = s.p0.y if s.p0.y < s.p1.y else s.p1.y
	x += offset, y += offset, _x -= offset, _y -= offset
	# compute number of lines and gaps
	xlen = x - _x, ylen = y - _y, A  =  xlen * ylen, O = A * infill
	xlines = floor(O / (extrudeWidth * xlen))
	xgap = (ylen - extrudeWidth * xlines) / xlines
	ylines = floor(O / (extrudeWidth * ylen))
	ygap = (xlen - extrudeWidth * ylines) / ylines
	# generate layers
	i = 0, z = least z of layers, R = listof(listof(line segments)), L = listof(line segments)
	switch = True
	while i < layers:
		if layers - i <= 1
			L = {	ls0(p0=(_x,_y,z),p1=(_x,y,z)) ,
				ls1(p0=(_x,y,z),p1=(x,y,z)) ,
				ls2(p0=(x,y,z),p1=(x,_y,z)) ,
				ls3(p0=(x,_y,z),p1=(_x,_y,z) }
		# at this point should call Aaron’s infill function at 100%  on L then combine
		# L and its infill to create a solid surface
		elif i % 2 == 1
			for k from 0 to xlines - 1
				if switch == True:
L  = {	ls0=(p0=(_x,_y+ygap * k,z),p1=(x,_y+ygap * k,z)) , 
ls1=(p0=(x,_y+ygap * k,z),p1=(x,_y+ygap*(k+1),z)) }
switch = False
				else:
					L = {	ls0=(p0=(x,_y+ygap * k,z),p1=(_x,_y+ygap* k,z)) ,
						ls1=(p0=(_x,_y+ygap *k,z)p1=(_x,_y+ygap * (k+1)))
					}
					switch = True
			
		elif i % 2 == 0
			for k from 0 to xlines -1
				if switch == True:
					L = {	ls0=(p0=(x-xgap * k, _y,z),p1=(x-xgap * k, y, z)) ,
ls1=(p0=(x-xgap * k,y,z),p1=(x-xgap * (k-1),y,z))
}
					switch = False
else:
L = {	ls0=(p0=(x-xgap * k,y,z),p1=(x-xgap * k,_y,z)) ,
	ls1=(p0=(x-xgap * k,_y,z),p1=(x-xgap * (k+1),_y,z)) 
}
switch = True
		R += {L}
switch = True
		z += sliceThickness
	return R
'''
