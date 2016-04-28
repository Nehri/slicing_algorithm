# -*- coding: utf-8 -*-

import math
import sys
import string
import copy

#printer specific constants, should be suplied as args
bedWidth = 150.0#mm
extrudeWidth = 0.71#mm
supportInfill = .5
delta = extrudeWidth/100.0 #delta for floating point comparison

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
        return "Point("+str(self.x)+","+str(self.y)+","+str(self.z)+")"
    def equals(self, p2):
        if close(self.x,p2.x) and close(self.y,p2.y) and close(self.z,p2.z):
            return True
        else:
            return False

def pointInLine(p, line):
    if close(p.x,line.p0.x) and close(p.y,line.p0.y) and close(p.z,line.p0.z):
        return True
    elif close(p.x,line.p1.x) and close(p.y,line.p1.y) and close(p.z,line.p1.z):
        return True
    else:
        return False

class Line:
    def __init__(self, p0_, p1_):
        self.p0 = p0_
        self.p1 = p1_

    def toString(self):
        return "Line("+self.p0.toString()+","+self.p1.toString()+")"
    def reverse(self):
        x_ = copy.copy(self.p0.x)
        y_ = copy.copy(self.p0.y)
        z_ = copy.copy(self.p0.z)
        self.p0.x = copy.copy(self.p1.x)
        self.p0.y = copy.copy(self.p1.y)
        self.p0.z = copy.copy(self.p1.z)
        self.p1.x = x_
        self.p1.y = y_
        self.p1.z = z_
        return self

#for floating point comparison
def close(f1,f2):
    comp = (max(f1,f2) - min(f1,f2))
    return (comp > -delta) and (comp < delta)

def lineEqual(L1,L2):
    if ((close(L1.p0.x, L2.p0.x) and close(L1.p0.y, L2.p0.y)
     and close(L1.p1.x, L2.p1.x) and close(L1.p1.y, L2.p1.y)) 
    or (close(L1.p0.x, L2.p1.x) and close(L1.p0.y, L2.p1.y) 
     and close(L1.p1.x, L2.p0.x) and close(L1.p1.y, L2.p0.y))):
        return True
    else:
        return False

class Triangle:
    def __init__(self, p0_, p1_, p2_, norm_):
        self.p0 = p0_
        self.p1 = p1_
        self.p2 = p2_
        self.norm = norm_
    def toString(self):
        return "Triangle("+self.p0.toString()+","+self.p1.toString()+","+self.p2.toString()+")"

def triangleEqual(T1,T2):
    if ((T1.p0.equals(T2.p0) and T1.p1.equals(T2.p1) and T1.p2.equals(T2.p2))
        or (T1.p0.equals(T2.p0) and T1.p1.equals(T2.p2) and T1.p2.equals(T2.p1))
        or (T1.p0.equals(T2.p1) and T1.p1.equals(T2.p0) and T1.p2.equals(T2.p2))
        or (T1.p0.equals(T2.p1) and T1.p1.equals(T2.p2) and T1.p2.equals(T2.p0))
        or (T1.p0.equals(T2.p2) and T1.p1.equals(T2.p0) and T1.p2.equals(T2.p1))
        or (T1.p0.equals(T2.p2) and T1.p1.equals(T2.p1) and T1.p2.equals(T2.p0))):
        return True
    else:
        return False

class Slice:
    def __init__(self, zValue_, perimeter_, isSurface_):
        self.zValue = zValue_
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
        points = list()
        for line in f:
            l_ = line.split(" ")
            l = [value for value in l_ if value != '']
            if counter == 6:
                counter = 0
                continue
            elif counter == 0:
                if l[0] == 'endsolid':
                    break
                points.insert(0, Point(float(l[2]), float(l[3]), float(l[4][:-1])))
            elif counter == 2:
                points.insert(0, Point(float(l[1]), float(l[2]), float(l[3])))
            elif counter == 3:
                points.insert(0, Point(float(l[1]), float(l[2]), float(l[3])))
            elif counter == 4:
                points.insert(0, Point(float(l[1]), float(l[2]), float(l[3])))
            counter += 1

        while points:
            triangles.insert(0, Triangle(points[2], points[1], points[0], points[3]))
            points = points[4:]

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
                return Point(x_=line.p0.x+t*slope.x, y_=line.p0.y+t*slope.y, z_=line.p0.z+t*slope.z)

            else: 
                return None
        else:
            return None


#helper for aboveTriangle
def sign(p1, p2, p3):
    return (p1.x-p3.x)*(p2.y-p3.y) - (p2.x-p3.x)*(p1.y-p3.y)

# given a point and a triangle, 
# returns True if that point is directly above that triangle,
# False otherwise
def aboveTriangle(point,triangle):
    if  (point.z > (triangle.p0.z-delta) and
        point.z > (triangle.p1.z-delta) and
        point.z > (triangle.p2.z-delta)):

        b1 = (sign(point, triangle.p0, triangle.p1) < 0.0)
        b2 = (sign(point, triangle.p1, triangle.p2) < 0.0)
        b3 = (sign(point, triangle.p2, triangle.p0) < 0.0)
        ret = ((b1 == b2) and (b2 == b3))
        return ret

    else:
        return False
        
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

            points_ = list(set([point1, point2, point3]))
            points = list()

            for point in points_:
                if point == None:
                    points_.remove(None)
                    break
            
            for i in range(0,len(points_)):
                j = i+1
                unique = True
                while j < len(points_):
                    if points_[i].equals(points_[j]):
                        unique = False
                    j+=1
                if unique:
                    points.insert(0,copy.deepcopy(points_[i]))

            if s <= (bounds[0]+layerThickness) or s >= (bounds[1]-layerThickness):
                currentSegmentSurface = True

            #if len(points) == 1:
                #currentSegment.append(Line(points[0], points[0]))
            if len(points) == 2:
                currentSegment.append(Line(points[0], points[1]))
            elif len(points) == 3:
                segment1 = Line(points[0], points[1])
                segment2 = Line(points[1], points[2])
                segment3 = Line(points[2], points[0])
                currentSegmentSurface = True

                currentSegment.append(segment1)
                currentSegment.append(segment2)
                currentSegment.append(segment3)
         
        segments.append(Slice(zValue_=s, perimeter_=copy.deepcopy(currentSegment),isSurface_=currentSegmentSurface))
        '''
        for line in currentSegment:
            print("appended "+line.toString())
        '''
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

        if ((intersect.x >= min(x1,x2)-delta) and (intersect.x <= max(x1,x2)+delta) and
            (intersect.y >= min(y1,y2)-delta) and (intersect.y <= max(y1,y2)+delta) and
            (intersect.x >= min(x3,x4)-delta) and (intersect.x <= max(x3,x4)+delta) and
            (intersect.y >= min(y3,y4)-delta) and (intersect.y <= max(y3,y4)+delta)):
            return intersect
        else:
            return None
        # return intersect
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

    numLines = int(round((bedWidth*percent)/extrudeWidth))
    gap = bedWidth/numLines
    infill = []

    #vertical lines
    for x in range(numLines):
        
        #start with full line

        fullLine = Line(Point((bedWidth/-2)+(x*gap),bedWidth/-2,Z),Point((bedWidth/-2)+(x*gap),bedWidth/2,Z))
        inters = []

        #find intersections without repeats
        for line in perimeter:
            sect = intersection(line,fullLine)
            if (sect != None):
                new = True
                for i in inters:
                    if close(i.y,sect.y):
                        new = False
                if new:
                    inters.append(copy.deepcopy(sect))

        #sort by y to get matching pairs for internal lines
        inters.sort(key=lambda point: point.y)
        
        if len(inters)%2 == 0: #if not even then something went wrong and its safer not to print
            for i in range(len(inters)):
                if i%2 != 0:
                    overlap = False;
                    newLine = Line(inters[i-1],inters[i])
                    for l in perimeter:
                        if lineEqual(l,newLine):
                            overlap = True;
                    if not overlap:
                        infill.append(newLine)
        '''
        else:
            print("Perimeter not manifold\n")
            print("fullLine: "+fullLine.toString())
            for line in perimeter:
                print(line.toString())
            print(" ")
            for p in inters:
                print(p.toString())
            print(" ")
        '''
    return infill


# given a list of line segments and a starting point,
# returns the location of the next line that connects to the points,
# or None if no point follows
def findNextPoint(point, lines):
    for i in range(0, len(lines)):
        line = lines[i]
        if pointInLine(point, line):
            return i
    return None

# given a slice with a list of line segments,
# returns a new slice free of duplicate or interior line segments
# and in order for optimized drawing
def cleanPerimeter(s):
    #for line in s:
        #if L is a duplicate and if every triangle containing L is on the slice, remove all L in base
    setPerimeter = copy.deepcopy(s.perimeter)
    
    i = 0
    while i < len(setPerimeter):
        j = i+1
        while j < len(setPerimeter):
            if lineEqual(setPerimeter[i],setPerimeter[j]):
                setPerimeter.remove(setPerimeter[j])
            else:
                j+=1
        i+=1


    '''    
    pathPerimeter = list()
    print("Perimetering")
    k = 0
    while setPerimeter:
        pathPerimeter.insert(0,copy.deepcopy(setPerimeter[0]))
        setPerimeter = setPerimeter[1:]
        while setPerimeter:
            print(len(pathPerimeter))
            loc = findNextPoint(pathPerimeter[k].p1, setPerimeter)
            if loc is None:
                #print(pathPerimeter[k].p1.toString())
                #for line in setPerimeter:
                #    print(line.toString())
                k+=1
                break
            if pathPerimeter[k].p1.equals(setPerimeter[loc].p0):
                pathPerimeter.insert(0,copy.deepcopy(setPerimeter[loc]))
            else:
                pathPerimeter.insert(0,copy.deepcopy(setPerimeter[loc].reverse()))
            setPerimeter.remove(setPerimeter[loc])
            k+=1
    
    for line in pathPerimeter:
        if line.p0.equals(line.p1):
            pathPerimeter.remove(line)
    finalPerimeter = [value for value in pathPerimeter if value != None]
    
    '''

    finalPerimeter = setPerimeter
    #need to order perimeter such that it is manifold
    return Slice(zValue_=s.zValue, perimeter_=finalPerimeter, isSurface_=s.isSurface)


# pseudocode for computing brim of a single convex polyhedron base
# takes a listof(line segments) which are the base (bottom layer),
# a number of outlines, and an initial offset
# will also generate a skirt if offset is greater than 1 and number of outlines = 1
def brim(base, numOutlines, offset):
    # compute centroid of polygonal manifold
    cx = 1
    cy = 1
    area = 1
    for line in base:
        cx = cx * (line.p0.x + line.p1.x)
        cy = cy * (line.p0.y + line.p1.y)
        area *= (line.p0.x * line.p1.y - line.p1.x * line.p0.y)
    area = area / 2
    cx = cx / (6 * area)
    cy = cy / (6 * area)

    # generate outlines
    brimlines = list()
    for i in range(1, numOutlines+1):
        for line in base:
            line_ = Line(Point(line.p0.x, line.p0.y, line.p0.z), Point(line.p1.x, line.p1.y, line.p1.z))

            if line.p0.x > cx:
                line_.p0.x += offset + extrudeWidth * i
            else:
                line_.p0.x -= offset + extrudeWidth * i
            if line.p0.y > cy:
                line_.p0.y += offset + extrudeWidth * i
            else:
                line_.p0.y -= offset + extrudeWidth * i
            if line.p1.x > cx:
                line_.p1.x += offset + extrudeWidth * i
            else:
                line_.p1.x -= offset + extrudeWidth * i
            if line.p1.y > cy:
                line_.p1.y += offset + extrudeWidth * i
            else:
                line_p1.y -= offset + extrudeWidth * i
            brimlines.append(line_)
    
    return brimlines


# pseudocode for computing a basic rectangular raft
# takes number of layers, an offset (how far raft extends from object), and an infill percentage
# and a listof(listof(line segments)) representing the object
# NOTE: raft goes UNDERNEATH object, all object / infill / support layers must be shifted up by
# height of raft (the brim, if used, is the same layer as raft and does not need to be elevated)
def raft(slices, numLayers, offset, infill, layerThickness):
    # compute enclosing rectangle on object
    x = 0
    y = 0
    x_ = 0
    y_ = 0
    for s in slices:
        for line in s.perimeter:
            if line.p0.x > x or line.p1.x > x:
                x = line.p0.x if line.p0.x > line.p1.x else line.p1.x
            if line.p0.x < x_ or line.p1.x < x_:
                x_ = line.p0.x if line.p0.x < line.p1.x else line.p1.x
            if line.p0.y > y or line.p1.y > y:
                y = line.p0.y if line.p0.y > line.p1.y else line.p1.y
            if line.p0.y < y_ or line.p1.y < y_:
                y_ = line.p0.y if line.p0.y < line.p1.y else line.p1.y

    x += offset
    y += offset
    x_ -= offset
    y_ -= offset
    # compute number of lines and gaps
    xlen = x - x_
    ylen = y - y_
    area  =  xlen * ylen
    totalArea = area * infill
    xlines = floor(totalArea / (extrudeWidth * xlen))
    xgap = (ylen - extrudeWidth * xlines) / xlines
    ylines = floor(totalArea / (extrudeWidth * ylen))
    ygap = (xlen - extrudeWidth * ylines) / ylines
    # generate layers
    i = 0
    z = slices[0].zValue
    allSegments = list()
    lines = list()
    switch = True
    while i < numLayers:
        if numLayers - i <= 1:
            lines_ = list(Line(p0_=Point(x_,y_,z),p1_=Point(x_,y,z)), Line(p0_=Point(x_,y,z),p1_=Point(x,y,z)), Line(p0_=Point(x,y,z),p1_=Point(x,y_,z)), Line(p0_=Point(x,y_,z),p1_=Point(x_,y_,z)))
            lines += lines_
            lines += infill(lines_, 1.0)
        # lines and its infill to create a solid surface

        elif i % 2 == 1:
            for k in range(0, xlines):
                if switch == True:
                    lines += list(Line(p0_=Point(x_,y_+ygap * k,z),p1_=Point(x,y_+ygap * k,z)), Line(p0_=Point(x,y_+ygap * k,z),p1_=Point(x,y_+ygap*(k+1),z)))
                    switch = False
                else:
                    lines += list(Line(p0_=Point(x,y_+ygap * k,z),p1_=Point(x_,y_+ygap* k,z)), Line(p0_=Point(x_,y_+ygap *k,z), p1_=Point(x_,y_+ygap * (k+1), z)))
                    switch = True
            
        elif i % 2 == 0:
            for k in range(0,xlines):
                if switch == True:
                    lines += list(Line(p0_=Point(x-xgap * k, y_,z),p1_=Point(x-xgap * k, y, z)), Line(p0_=Point(x-xgap * k,y,z),p1_=Point(x-xgap * (k-1),y,z)))
                    switch = False
        else:
            lines  += list(Line(p0_=Point(x-xgap * k,y,z),p1_=Point(x-xgap * k,y_,z)), Line(p0_=Point(x-xgap * k,y_,z),p1_=Point(x-xgap * (k+1),y_,z))) 
            switch = True
        allSegments.append(lines)
        switch = True
        z += layerThickness
    return allSegments

# given a list of triangles,
# returns a list of any downward-facing triangles
def downward(triangles):
    trianglesDown = list()
    for triangle in triangles:
        if triangle.norm.z < 0:
            trianglesDown.insert(0, copy.deepcopy(triangle))
    return trianglesDown

# given a downward-facing triangle and a list of all triangles,
# returns True if no triangles are in the way of the downward triangle
# and False if a triangle blocks the path directly downward
def supportNeeded(triangle, triangles, bottomZ):
    if (close(triangle.p0.z, bottomZ)
        and close(triangle.p1.z, bottomZ)
        and close(triangle.p2.z, bottomZ)):
        return False

    for tri in triangles:
        if (aboveTriangle(triangle.p0, tri) 
            or aboveTriangle(triangle.p1, tri) 
            or aboveTriangle(triangle.p2, tri)):
            return False

    return True


# given a triangle that requires support and a minimum Z value,
# returns the list of triangles required to form the support shape under that triangle
def generateSupportShape(triangle, bottomZ):
    triangleTop = copy.deepcopy(triangle)
    triangleBottom = Triangle(Point(triangleTop.p0.x, triangleTop.p0.y, bottomZ), 
                                Point(triangleTop.p1.x, triangleTop.p1.y, bottomZ), 
                                Point(triangleTop.p2.x, triangleTop.p2.y, bottomZ), None)
    newShape = [triangleTop]
    newShape.insert(0, triangleBottom)

    newShape.insert(0, Triangle(triangleTop.p0, triangleTop.p1, triangleBottom.p0, None))
    newShape.insert(0, Triangle(triangleTop.p1, triangleTop.p2, triangleBottom.p1, None))
    newShape.insert(0, Triangle(triangleTop.p2, triangleTop.p0, triangleBottom.p2, None))
    newShape.insert(0, Triangle(triangleBottom.p0, triangleBottom.p1, triangleTop.p1, None))
    newShape.insert(0, Triangle(triangleBottom.p1, triangleBottom.p2, triangleTop.p2, None))
    newShape.insert(0, Triangle(triangleBottom.p2, triangleBottom.p0, triangleTop.p0, None))

    i = 0
    while i < len(newShape):
        j = i+1
        while j < len(newShape):
            if triangleEqual(newShape[i],newShape[j]):
                newShape.remove(newShape[j])
            else:
                j+=1
        i+=1

    return newShape

# given a list of triangles
# returns a list of list of slices to draw supports for any downward-facing triangles
#returns a list of lists of slices
def generateSupports(triangles, layerThickness):

    bounds = findBoundaries(triangles)

    trianglesDown = downward(triangles)

    trianglesForSupport = list()
    for tri in trianglesDown:
        if supportNeeded(tri, triangles, bounds[0]):
            trianglesForSupport.insert(0,copy.deepcopy(tri))

    supportShapes = list()
    for triangle in trianglesForSupport:
        supportShapes.insert(0, generateSupportShape(triangle, bounds[0]))

    supportSlices = list()

    for shape in supportShapes:
        supportSlices.insert(0, separateSlices(shape, layerThickness))

    return supportSlices

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

    o = bedWidth/2 #origin
    layer = 1; #current layer/slice
    E = 0; #extrusion accumulator
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
            f.write("G0 F2700 X"+str(o+l.p0.x)+" Y"+str(o+l.p0.y)+" Z"+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            E += dist*extrudeRate
            f.write("G1 F900 X"+str(o+l.p1.x)+" Y"+str(o+l.p1.y)+" E"+str(E)+"\n")

        if len(s.support) > 0:
            f.write(";support\n")
        for l in s.support:
            #move to start of line
            f.write("G0 F2700 X"+str(o+l.p0.x)+" Y"+str(o+l.p0.y)+" Z"+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            E += dist*extrudeRate
            f.write("G1 F800 X"+str(o+l.p1.x)+" Y"+str(o+l.p1.y)+" E"+str(E)+"\n")

        if len(s.infill) > 0:
            f.write(";infill\n")
        for l in s.infill:
            #move to start of line
            f.write("G0 F2700 X"+str(o+l.p0.x)+" Y"+str(o+l.p0.y)+" Z"+str(l.p0.z)+"\n")
            #move to end while extruding
            dist = math.sqrt(pow(l.p1.x-l.p0.x,2) + pow(l.p1.y-l.p0.y,2))
            E += dist*extrudeRate
            f.write("G1 F900 X"+str(o+l.p1.x)+" Y"+str(o+l.p1.y)+" E"+str(E)+"\n")
        
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
    infillPercent = float(sys.argv[3])
    triangles = fileToTriangles(filename)

    slices_ = separateSlices(triangles, layerThickness)
    supportSlices = generateSupports(triangles, layerThickness)

    slices = list()

    for s in slices_:
        slices += [cleanPerimeter(s)]

    for s in slices:
        if s.isSurface:
            s.infill = infill(s.perimeter, 1)
        else:
            s.infill = infill(s.perimeter, infillPercent)

    for shape in supportSlices:
        for s in range(len(shape)):
            slices[s].support += infill(shape[s].perimeter,supportInfill)

    writeGcode(slices,filename)

if __name__ == "__main__":
    main()
