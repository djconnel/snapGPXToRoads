"""GPXUtils is a series of procedures useful for processing GPX files,
using the gpxp:y library.
version histry:
0.01: first release
0.02: additional methods
0.03: added arcFit
0.04: added calcDistances, cropCorners, smoothPoints
0.05: added minRadius, fixed bug w/ S/F croppingm, zigZag and loop checks, renamed gpxutils from GPXUtils
author: Daniel Connelly djconnel@gmail.com
"""

__version__ = 0.05

from math import pi, cos, sin, tan, sqrt, atan2, floor, ceil, exp
import gpxpy
import gpxpy.gpx
import sys
from sys import stderr
rEarth = 20037392 / pi

twopi   = 2 * pi
deg2rad = pi / 180
lat2y   = rEarth * deg2rad


def loadGPX(gpx_file):
  """
  Load GPX file and extract all track points.

  Args:
    gpx_file: Path to the GPX file

  Returns:
    List of (latitude, longitude) tuples
  """

  gpx = gpxpy.parse(gpx_file)

  points = []
  track = gpx.tracks[0]
  for segment in track.segments:
    for point in segment.points:
      points.append(point)

  stderr.write(f"Loaded {len(points)} points from GPX file\n")
  return points


def copyPoint(p):
  return gpxpy.gpx.GPXTrackPoint(
    latitude= p.latitude,
    longitude= p.longitude,
    elevation= p.elevation
  )

def reduceAngle(theta):
  theta -= 2 * pi * floor(0.5 + theta / (2 * pi))
  return theta

def reduceDegs(theta):
  theta -= 360 * floor(0.5 + theta / 360)
  return theta

def pointsAreClose(p1, p2, deltaDegs= 1e-6, deltaZ= 0.1):
  if (
    (p1.latitude is None) or
    (p2.latitude is None) or
    (p1.longitude is None) or
    (p2.longitude is None) or
    (abs(p1.latitude - p2.latitude) > deltaDegs) or
    ((deltaZ is not None ) and abs(p1.elevation - p2.elevation) > deltaZ)
  ):
    return False
  dLon = p1.longitude - p2.longitude
  dLon -= 360 * floor(dLon / 360 + 0.5)
  return (abs(dLon) < deltaDegs)

def pointsAreVeryClose(p1, p2, deltaDegs= 1e-8, deltaZ= 0.001):
  return pointsAreClose(p1, p2, deltaDegs = deltaDegs, deltaZ = deltaZ)

def pointsAreCoincident(p1, p2, deltaDegs= 1e-7):
  return pointsAreClose(p1, p2, deltaDegs= deltaDegs, deltaZ = None)

def pointdxdy(p1, p2):
  """
  Calculate the great circle distance in kilometers between two points
  on the earth (specified in decimal degrees)
  """
  # convert decimal degrees to radians
  lon1 = p1.longitude * deg2rad
  lat1 = p1.latitude * deg2rad
  lon2 = p2.longitude * deg2rad
  lat2 = p2.latitude * deg2rad

  # haversine formula
  dlon = lon2 - lon1
  dlon -= twopi * floor(0.5 + dlon / twopi)
  dlat = lat2 - lat1

  u = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
  v = sin(dlon) * cos(lat2)
  uv = sqrt(u ** 2 + v ** 2)
  if uv == 0:
    dx = 0
    dy = 0
  else:
    s = u / uv
    c = v / uv
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    d = 2 * rEarth * atan2( sqrt(a), sqrt(1 - a) )
    dx = d * c
    dy = d * s
  return (dx, dy)

def pointDistance(p1, p2):
  """
  Calculate the great circle distance in kilometers between two points
  on the earth (specified in decimal degrees)
  """
  # convert decimal degrees to radians
  lon1 = p1.longitude * deg2rad
  lat1 = p1.latitude * deg2rad
  lon2 = p2.longitude * deg2rad
  lat2 = p2.latitude * deg2rad

  # haversine formula
  dlon = lon2 - lon1
  dlon -= twopi * floor(0.5 + dlon / twopi)
  dlat = lat2 - lat1

  a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
  d = 2 * rEarth * atan2( sqrt(a), sqrt(1 - a) )
  return d

def deltaPosition(d1, d2, courseDistance=None, isLoop=False):
  delta = d2 - d1
  if isLoop and courseDistance is not None and courseDistance > 0:
    delta -= courseDistance * floor(0.5 + delta / courseDistance)
  return delta

def calcDistances(points, isLoop):
  if len(points) == 0:
    return []
  distances = [0]
  for i in range(len(points) - 1):
    distances.append(distances[-1] + pointDistance(points[i], points[i + 1]))
  courseDistance = distances[-1]
  if isLoop:
    courseDistance += pointDistance(points[-1], points[0])
    
  return distances, courseDistance

def reduceAngle(theta):
  theta -= twopi * floor(0.5 + theta / twopi)
  return theta;

def segmentDirection(p1, p2):
  dx, dy = pointdxdy(p1, p2)
  return atan2(dy,dx)

def averageAngles(theta1, theta2):
  return reduceAngle( theta1 + 0.5 * reduceAngle(theta2 - theta1) )

def pointDirection(p1, p2, p3):
  return averageAngles(segmentDirection(p1, p2), segmentDirection(p2, p3));

def calcDirections(points= [], isLoop = []):
  N = len(points)
  if N == 0:
    return points
  directions = [None] * N
  u = N - 1 if isLoop else 0
  v = 0
  w = 1
  dPrev = None
  while v < N:
    u = v
    w = v
    while pointsAreClose(points[u], points[v]):
      if (((u - 1) % N != w) if isLoop else u > 0):
        u = (u - 1) % N
      else:
        u = v
        break
    while pointsAreClose(points[w], points[v]):
      if (w + 1) % N != u if isLoop else w < N - 1 :
        w = (w + 1) % N
      else:
        w = v
        break

    d = dPrev if dPrev is not None else 0
    if u == v:
      if v != w:
        d = segmentDirection(points[v], points[w])
    else:
      if v == w:
        d = segmentDirection(points[u], points[v])
      else:
        d = pointDirection(points[u], points[v], points[w])
    if dPrev is not None:
      d = dPrev + reduceAngle(d - dPrev)
    dPrev = d
    directions[v] = d
    v += 1
  return directions


def distanceInRange(distance, start = None, stop = None, isLoop = None):
  """test whether the given distance is in the given range"""
  if isLoop and start is not None and stop is not None and stop < start:
    return start <= distance or stop >= distance
  else:
    return (start is None or start <= distance) and (stop is None or stop >= distance)

# add a vector to a point
# a single iteration will use the average cosine for the path rather than
# a starting cosine, for 1st order improvement
def addVectorToPoint(point, v):
  dx = v[0]
  dy = v[1]
  dz = v[2] if len(v) > 2 else 0
  pNew = copyPoint(point)
  dlat = dy / lat2y
  pNew.latitude += dlat
  if abs(pNew.latitude) > 90:
    sys.exit(0)
  c = cos(deg2rad * pNew.latitude + dlat)
  pNew.longitude += dx / c / lat2y
  pNew.longitude -= 360 * floor(0.5 + pNew.longitude / 360)
  pNew.elevation += dz
  return pNew

def pointDotCrossProduct(p1, p2, p3, p4):
  x1, y1 = pointdxdy(p1, p2)
  x2, y2 = pointdxdy(p3, p4)
  dot = x1 * x2 + y1 * y2
  cross = x1 * y2 - x2 * y1
  return (dot, cross)

def pointNormalizedDotCross(p1, p2, p3, p4):
  x1, y1 = pointdxdy(p1, p2)
  x2, y2 = pointdxdy(p3, p4)
  dot = x1 * x2 + y1 * y2
  cross = x1 * y2 - x2 * y1
  d = sqrt(dot ** 2 + cross ** 2)
  if d == 0:
    return (1, 0)
  else:
    return (dot / d, cross / d)

def pointdDot(p1, p2, p3, p4):
  dot, cross = pointDotCross(p1, p2, p3, p4)
  return dot

def pointCross(p1, p2, p3, p4):
  dot, cross = pointDotCross(p1, p2, p3, p4)
  return cross

def pointNormalizedDot(p1, p2, p3, p4):
  dot, cross = pointNormalizedDotCross(p1, p2, p3, p4)
  return dot

def pointNormalizedCross(p1, p2, p3, p4):
  dot, cross = pointNormalizedDotCross(p1, p2, p3, p4)
  return cross

def pointAngle(p1, p2, p3):
  dot, cross = pointNormalizedDotCross(p1, p2, p2, p3)
  angle = atan2(cross, dot)
  return angle

def xyPointOnLine(xy1, xy2, xy3):
  # x,y points
  x1, y1 = xy1
  x2, y2 = xy2
  x3, y3 = xy3
  if x1 == x2 and y1 == y2:
    return (None, None)
  f = ( (y3 - y1) * (y2 - y1) + (x3 - x1) * (x2 - x1) ) / ( (y2 - y1) ** 2 + (x2 - x1) ** 2 )
  d = sqrt( ( x1 - x3 + f * (x2 - x1) ) ** 2 + ( y1 - y3 + f * (y2 - y1) ) ** 2 )
  return ( f, d )

def removeDuplicatePoints(points):
  newPoints = [ points[0] ]
  for p in points:
    if not pointsAreVeryClose(p, newPoints[-1]):
      newPoints.append(p)
  return newPoints

def averageCoincidentPoints(points):
  """average coincident points, not wrapping around for loop courses"""
  newPoints = []
  i = 0
  while (i < len(points)):
    pi = points[i]
    newPoints.append(copyPoint(pi))
    j = i
    sumLat = points[i].latitude
    sumLon = points[i].longitude
    sumEle = points[i].elevation
    while ((j + 1 < len(points)) and pointsAreCoincident(pi, points[j + 1])):
      j += 1
      dLon = points[j].longitude - points[i].longitude
      dLon -= 360 * floor(0.5 + dLon / 360)
      sumLat += points[j].latitude
      sumLon += points[i].longitude + dLon
      sumEle += points[j].elevation
    if j > i:
      N = 1 + j - i
      newPoints[-1].longitude = sumLon / N
      newPoints[-1].latitude = sumLat / N
      newPoints[-1].elevation = sumEle / N
      newPoints[-1].longitude -= 360 * floor(0.5 + newPoints[-1].longitude / 360)
    i = j + 1
  return newPoints

def deltaPosition(d1, d2, courseDistance=None, isLoop=False):
  delta = d2 - d1
  if isLoop and courseDistance is not None and courseDistance > 0:
    delta -= courseDistance * floor(0.5 + delta / courseDistance)
  return delta

def interpolatePosition(d1, d2, f = 0.5, courseDistance = None, isLoop = False):
  dd = d2 - d1
  if isLoop:
    dd -= courseDistance * floor(0.5 + dd / courseDistance)
  d = d1 + dd * f
  if isLoop is courseDistance is not None:
    d -= courseDistance * floor(0.5 + d / courseDistance)
  return d

def applyUniformGradient(points, bow = None):
  """create a uniform gradient across the points"""
  distances = [0]
  for i in range(1, len(points)):
    distances.append(distances[-1] + pointDistance(points[i - 1], points[i]))
  if distances[-1] == 0:
    distances = [i for i in range(len(points))]
  z0 = points[0].elevation
  z1 = points[-1].elevation
  if bow == 0:
    bow = None
  for i in range(1, len(points) - 1):
    f = distances[i] / distances[-1]
    points[i].elevation = (1 - f) * z0 + f * z1
    if bow is not None:
      points[i].elevation += 4 * f * (1 - f) * bow

def interpolateLon(lon1, lon2, f = 0.5):
    dl = lon2 - lon1
    dl -= 360 * floor(0.5 + dl / 360)
    l = lon1 + f * dl
    l -= 360 * floor(0.5 + l / 360)
    return l

def interpolateLat(lat1, lat2, f = 0.5):
    dl = lat2 - lat1
    l = lat1 + f * dl
    return l

def interpolatePoint(p1, p2, f):
  p = copyPoint(p1)
  p.latitude = interpolateLat(p1.latitude, p2.latitude, f)
  p.longitude = interpolateLon(p1.longitude, p2.longitude, f)
  p.elevation = (1 - f) * p1.elevation + f * p2.elevation
  return p

def arcFit(p0, p1, p2, p3, rMax = None, maxRadius = 20, maxAngle = pi / 16):
  """fit an arc from p2 to p3 tangent to p1-p2 and p3-p4 with maximum radius maxRadius"""
  if rMax == 0:
    return []

  # calculate the intersection point: p1 -> d1
  d01  = pointDistance(p0, p1)
  d23  = pointDistance(p2, p3)
  xy01 = pointdxdy(p0, p1)
  xy12 = pointdxdy(p1, p2)
  xy23 = pointdxdy(p2, p3)
  
  c1 = xy01[0] / d01
  s1 = xy01[1] / d01
  c2 = xy23[0] / d23
  s2 = xy23[1] / d23
  cross = c1 * s2 - c2 * s1

  if abs(cross) < 1e-6:
    return [] 
  # sign positive: right turn
  # sign negative: left turn
  sign = 1 if cross < 0 else -1 if cross > 0 else 0

  # list of points between p1 and p2, exclusive
  points = []

  # c1, s1 are approaching intersection, so flip sign
  # c2, s2 are leaving intersection, so keep sign
  dot = c1 * c2 + s1 * s2
  # angle subtended by the corner: 0 => 180 degree turn, 180 => straight
  theta = reduceAngle( atan2(-cross, -dot) )
  # multiplication factor for translating radius to position
  fR  = abs(1 / tan(theta / 2))

  # distance to intercept from p1, p2
  # di1 = [ s2 (x2 - x1) - c2 (y2 - y1) ] / (c1 s2 - c2 s1)
  di1 = ( s2 * xy12[0] - c2 * xy12[1] ) / cross
  di2 = ( c1 * xy12[1] - s1 * xy12[0] ) / cross

  # if there's a negative value, bail (revisit this)
  if ((di1 < 0) or (di2 < 0)):
    return [] 

  # figure out radii from distance
  R = min(di1, di2) / fR
  if (maxRadius is not None) and (maxRadius > 0) and (R > maxRadius):
    R = maxRadius 

  a = di1 - R * fR
  b = di2 - R * fR

  # distance from p1 to new corner point
  if abs(a) > 0.01:
    f = a / (2 * R + a)
    dz = (p2.elevation - p1.elevation) * f
    p1 = addVectorToPoint(p1, [a * c1, a * s1, dz])
    points.append(p1)
    xy01 = pointdxdy(p0, p1) # new p1
    xy12 = pointdxdy(p1, p2) # new p1

  if abs(b) > 0.01:
    # distance from p2 to new corner point
    f = b / (2 * R + b)
    dz = (p1.elevation - p2.elevation) * f
    px = addVectorToPoint(p2, [-b * c2, -b * s2, dz])
    p2 = px
    xy12 = pointdxdy(p1, p2) # new p2
    xy23 = pointdxdy(p2, p3) # new p2
  else:
    px = None

  cx = sign * s1 * R
  cy = -sign * c1 * R

  theta1 = atan2(-cy, -cx)
  theta2 = atan2(xy12[1] - cy, xy12[0] - cx)

  dTheta = reduceAngle(theta2 - theta1)
  nPoints = int(abs(dTheta / maxAngle))

  z1 = p1.elevation
  z2 = p2.elevation
  deltaZ = z2 - z1 if (z1 is not None) and (z2 is not None) else 0
  for n in range(1, nPoints + 1):
    f = n / (nPoints + 1)
    theta = theta1 + dTheta * f
    dx = cx + R * cos(theta)
    dy = cy + R * sin(theta)
    dz = f * deltaZ
    points.append(addVectorToPoint(p1, [dx, dy, dz]))

  # if we need to add the extra point, do it here
  if px is not None:
    points.append(px)

  applyUniformGradient([p1] + points + [p2])

  return points

def UTurnCheck (p1, p2, p3, p4, maxCos = -0.98):
  d = pointNormalizedDot(p1, p2, p3, p4)
  return ((d is not None) and (d < maxCos))


def fixZigZags(points = [], isLoop= False):
  '''eliminate "zig-zags" from a list of GPX points'''

  N = len(points)
  if N < 3:
    return points

  for zigZagIter in range(10):
    dzigZag = 100
    UTurns = []
    iStart = 0 if isLoop else 1
    iEnd = len(points) - 1 if isLoop else len(points) - 2
    for i in range(iStart, iEnd + 1):
      if UTurnCheck(points[i - 1], points[i], points[i], points[(i + 1) % N], -0.9):
        UTurns.append(i)

    distances, courseLoop = calcDistances(points= points, isLoop= isLoop)
    zigZagCount = 0
    if UTurns:
      while len(UTurns) > 1:
        U1 = UTurns.pop(0)
        U2 = UTurns[0]
        p1 = points[U1]
        p2 = points[U2]
        d1 = distances[U1]
        d2 = distances[U2]
        if d2 - d1 < dzigZag:
          stderr.write(f'WARNING: zig-zag found on points : {U1} and {U2}, from {d1:g} m to {d2:g} m separated by {(d2 - d1):g} m\n')
          # repairing zig-zags...
          # zig-zags are two U-turns within a specified distance
          # p1  p2  ...   p3  p4
          # U-turn @ p2, and U-turn @ p3
          # 1. eliminate all points between p2 and p3
          # 2. as long as P3 has a U-turn, delete it... there will be a new P3
          # 3. as long as P2 has a U-turn, delete it...
          # 4. go back step 2 if we deleted any U-turns

          # eliminate points between
          u = U1
          v = U2 + 1
          while v < len(points) - 1 and UTurnCheck(points[u], points[v], points[v], points[v + 1]) :
            v += 1


          stderr.write(f'eliminating {v - u - 1} points\n')
          zigZagCount += 1
          pNew =  points[0:u + 1] + points[v: N]
          # if we ran out of points, something is wrong
          if len(pNew) < 2:
            stderr.write('repairing zig-zags eliminated entire route\n')
            sys.exit()

          points = pNew
          N = len(points)
          # adjust U-turn coordinates
          # We've eliminated the next Uturn, so pop it
          UTurns.pop(0)
          # adjust coordinates of remaining U-turns
          UTurns = [ uturn + u - v + 1 for uturn in UTurns]

          # get rid of obsolete U-turns
          while (UTurns and UTurns[0] < 0):
            UTurns.pop(0)
    if zigZagCount == 0:
      break
  return points

def findLoops(points= [], isLoop= False, loopDistance= 100, distances= None, courseDistance= None):
  N = len(points)
  if distances is None or courseDistance is None:
    distances, courseDistances = calcDistances(points= points, isLoop= isLoop)
  directions = calcDirections(points= points, isLoop= isLoop)
  u = 0
  v = 0
  loopAngle = 0.85 * twopi
  while (v < N - 1):
    p = points[u]
    d = distances[u]
    direction = directions[u]
    while (v < N - 1) and (distances[v + 1] < d + loopDistance):
      v += 1
    if (abs(direction - directions[v]) > loopAngle):
      while u + 1 < v and abs(directions[u + 1] - directions[v]) > loopAngle:
        u += 1
        while v - 1 > v and abs(directions[v - 1] - directions[u]) > loopAngle:
          v -= 1 
      stderr.write(f'WARNING: oop between distance: {distance[u] / 1000:.3f} km to {distance[v] / 1000:.3f} km\n')
      u = v
      next
    u += 1

def simplifyPoints(points, z0 = 0.2, r0 = 1):
  N = len(points)
  if r0 <= 0 or N < 3:
    return points

  xf, yf = pointdxdy(points[0], points[-1])

  if xf ** 2 + yf ** 2 < 10:
    iFurthest = None
    dFurthest = 0
    for i in range(1, N):
      d = pointDistance(points[0], points[i])
      if d > dFurthest:
        iFurthest = i
        dFurthest = d
    if iFurthest is not None:
      if iFurthest == N - 1:
        return [points[0], points[-1]]
      else:
        p1 = simplifyPoints(points = points[:iFurthest + 1], z0  = z0, r0 = r0)
        p2 = simplifyPoints(points = points[iFurthest:], z0  = z0, r0 = r0)
        return p1[1:] + p2
    else:
      return points

  zi = points[0].elevation
  dzf = points[-1].elevation - zi
  iMax = None
  scoreMax = 1  # only accept points if score is at least 1
  for i in range(1, N):
    x, y = pointdxdy(points[0], points[i])
    dz = points[i].elevation - zi
    # find the nearest point on the curve
    f, d = xyPointOnLine((0, 0), (xf, yf), (x, y))   # distance of the interpolated point
    ddz = dzf * f - dz                               # interpolated altitude difference
    score = (d / r0) ** 2 + (ddz / z0) ** 2
    if score > scoreMax:
      iMax = i
      scoreMax = score

  if iMax is not None:
    p1 = simplifyPoints(points= points[:iMax + 1], z0 = z0, r0 = r0)
    p2 = simplifyPoints(points= points[iMax:], z0 = z0, r0 = r0)
    return p1 + p2[1:]
  else:
    return [points[0], points[-1]]

def checkLoop(points):
  """check that the route thru the S/F is not a U-turn"""
  i1 = -1
  i2 = 0
  i3 = 1
  while ((pointDistance(points[i1], points[i2]) < 10) and (i1 > i2 - len(points) + 1)):
    i1 -= 1
  while ((pointDistance(points[i2], points[i3]) < 10) and (i3 < len(points) - 1)):
    i3 += 1
  isLoop = (
         (len(points) > 2) and
         (pointDistance(points[0], points[-1]) < 150) and
         (pointNormalizedDot(points[i1], points[i2], points[i2], points[i3]) > -0.1)
  )
  return isLoop

def fixLoopAltitudes(points):
  '''make sure altitudes are smooth at S/F for loop course'''
  ds, courseDistance = calcDistances(points= points, isLoop = True)

  f = (ds[-1] - ds[-2]) / courseDistance
  z = points[-2].elevation * (1 - f) + points[0].elevation * f
  dz = z - points[-1].elevation
  for i in range(len(points)):
    f = ds[i] / courseDistance
    points[i].elevation += dz * f

  
def cropCorners(points, cropping = None, isLoop = False, start = None, end = None):
  if cropping is None or cropping <= 0:
    return points

  # remove duplicate points at end (can re-add them later)
  removedPoints = []
  closedLoop = False
  if isLoop:
    while len(points) > 0 and pointsAreClose(points[0], points[-1]):
      closedLoop = True
      removedPoints.append(points.pop())
    removedPoints.reverse()

  N = len(points)
  if N == 0:
    return points

  ds, courseDistance = calcDistances(points= points, isLoop = isLoop)

  # find corner candidates
  cornerPoints = []
  for i in range(N):
    p = points[i]
    if not isLoop and (i == 1 or i == N - 1):
      continue
    j = (i - 1) % N
    while ((isLoop and j > 0) or ((not isLoop) and j != (i + 1) % N)) and pointsAreVeryClose(points[j], points[i]):
      j = (j - 1) % N
    k = (i + 1) % N
    while ((isLoop and k < N - 1) or ((not isLoop) and k != (i - 1) % N)) and pointsAreVeryClose(points[k], points[i]):
      k = (k + 1) % N
    p1 = points[j]
    p2 = points[k]
    d1 = ds[j]
    d2 = ds[k]
    d  = ds[i]
    dd1 = deltaPosition(d1, d, courseDistance= courseDistance, isLoop= isLoop)
    dd2 = deltaPosition(d, d2, courseDistance= courseDistance, isLoop= isLoop)
    if dd1 == 0 or dd2 == 0:
      stderr.write(f'coincident points found!\n')
      sys.exit()
    cross = pointNormalizedCross(p1, p, p, p2)
    if abs(cross) > 0.5:
      cornerPoints.append(i)

  # corners which are too close get pruned
  dMin = 2 * cropping
  cc = []
  for ic in range(len(cornerPoints)):
    d = ds[cornerPoints[ic]]
    if not isLoop and ic == 0:
      dPrev = 0
    else:
      dPrev = ds[cornerPoints[ic - 1]]
      if ic == 0:
        dPrev -= courseDistance
    if not isLoop and ic == len(cornerPoints) - 1:
      dNext = ds[-1]
    else:
      dNext = ds[cornerPoints[(ic + 1) % len(cornerPoints)]]
      if ic == len(cornerPoints) - 1:
        dNext += courseDistance

    # pass criteria:
    # 1. point is in limits defined by start and/or stop (if isLoop, start can be > stop)
    # 2. corner is sufficiently far from neighbor corners (or for point-to-point, course terminals)
    flag = True
    if isLoop and start is not None and end is not None and end < start:
      flag &= d >= end and d <= start
    else:
      if start is not None:
        flag &= d >= start
      if end is not None:
        flag &= d <= end
    flag &= dNext > d + dMin and dPrev <= d - dMin
    if flag:
      cc.append(cornerPoints[ic])
      
  cornerPoints = cc

  if len(cornerPoints) == 0:
    return points + removedPoints

  # keep track of corner points with pruning
  isCornerPoint = [False] * N
  for ic in cornerPoints:
    isCornerPoint[ic] = True

  skipPoint = [False] * N

  for i in range(N):
    if not isLoop and (i == 1 or i == N - 1):
      continue

    # this is corner to be cropped
    if isCornerPoint[i]:
      d0 = ds[i]
      dc1 = d0 - cropping + 0.01
      dc2 = d0 + cropping - 0.01

      # find the interval spanning the crop point
      j = (i - 1) % N
      while j != i if isLoop else 0:
        d1 = ds[j]
        while d1 > d0:
          d1 -= courseDistance
        if d1 < dc1:
          break
        skipPoint[j] = True
        j = (j - 1) % N

      j = (i + 1) % N
      while j != (i if isLoop else len(points) - 1):
        d1 = ds[j]
        while d1 < d0:
          d2 += courseDistance
        if d1 > dc2:
          break
        skipPoint[j] = True
        j = (j + 1) % N

  newPoints = []
  endPoints = []

  for i in range(N):
    p = points[i]
    d = ds[i]
    if skipPoint[i]:
      continue
    elif isCornerPoint[i]:
      j = (i - 1) % N
      while skipPoint[j]:
        j = (j - 1) % N
      p1a = points[j]
      d1a = ds[j]
      p1b = points[(j + 1) % N]
      d1b = ds[(j + 1) % N]
      k = (i + 1) % N
      while skipPoint[k]:
        k = (k + 1) % N
      p2a = points[k - 1]
      d2a = ds[k - 1]
      p2b = points[k]
      d2b = ds[k]
        
      # back up until we reach the crop distance
      dd1a = deltaPosition(d1a, d, courseDistance= courseDistance, isLoop= isLoop)
      dd1b = deltaPosition(d1b, d, courseDistance= courseDistance, isLoop= isLoop)
      dd2a = deltaPosition(d, d2a, courseDistance= courseDistance, isLoop= isLoop)
      dd2b = deltaPosition(d, d2b, courseDistance= courseDistance, isLoop= isLoop)
      
      f1 = (dd1a - cropping) / (dd1a - dd1b)
      f2 = (cropping - dd2a) / (dd2b - dd2a)
      addPoints = []
      if abs(dd1b - cropping) > 0.01:
        p1b = interpolatePoint(p1a, p1b, f1)
        if k < j:
          endPoints = [p1b]
        else:
          newPoints.append(p1b)
      if abs(dd2a - cropping) > 0.01:
        p2a = interpolatePoint(p2a, p2b, f2)
        addPoints.append(p2a)

      pArc = arcFit(p1a, p1b, p2a, p2b, maxRadius = cropping)
      
      if k < j:
        # this wraps around, so either first point is a corner point, or the last point
        removedPoints = []
        iMid = (len(pArc) - 1) // 2
        pArc1 = pArc[:iMid]
        pArc2 = pArc[iMid:]
        if i == 0:
          endPoints += pArc1
          newPoints = pArc2 + addPoints
        else:
          newPoints = pArc2 + addPoints + newPoints + pArc1
      else:
        newPoints += pArc + addPoints
    else:
      newPoints.append(p)

  newPoints += endPoints + removedPoints

  if closedLoop and not pointsAreClose(newPoints[-1], newPoints[0]):
    newPoints.append(copyPoint(newPoints[0]))
    
  stderr.write(f'points after cropping = {len(newPoints)}\n')
  return newPoints

def mergeRanges(points, ranges, isLoop, courseDistance = None):
  """
  Merge overlapping ranges into non-overlapping ranges.
  
  Args:
    ranges: List of lists, where each tuple is [start, end]
  
  Returns:
    List of non-overlapping ranges sorted by start position
  """
  # Handle empty input
  if not ranges:
    return []

  # handle wrap-around ranges: second element can overflow end of list
  N = len(points)
  if isLoop:
    for r in ranges:
      r[0] -= courseDistance * floor(r[0] / courseDistance)
      r[1] -= courseDistance * floor((r[1] - r[0]) / courseDistance)

  # Sort ranges by start position
  sorted_ranges = sorted(ranges, key=lambda x: x[0])
  
  # Initialize result with first range
  merged = [sorted_ranges[0]]

    # Iterate through remaining ranges
  for current_start, current_end in sorted_ranges[1:]:
    last_start, last_end = merged[-1]
    
    # Check if current range overlaps with last merged range
    if current_start <= last_end:
      # Merge by extending the end if necessary
      merged[-1] = [last_start, max(last_end, current_end)]
    else:
      # No overlap, add as new range
      merged.append([current_start, current_end])

  # handle loop: this is overkill, to be safe
  if isLoop:
    newRanges = []
    if courseDistance is None:
      distances, courseDistance = calcDistances(points, isLoop = isLoop)
    for r in merged:
      if r[1] > courseDistance:
        newRanges.append([0, r[1] - courseDistance])
        r[1] = courseDistance
    while len(newRanges) > 1:
      stop = max(newRanges[0][1], newRanges[1][1])
      newRanges = [0, stop] + newRanges[2:]
    if len(newRanges) > 0:
      while len(merged) > 0 and newRanges[0][1] > merged[0][0]:
        newRanges[0][1] = max(merged[0][1], newRanges[0][1])
        merged = merged[1:]
      merged = newRanges + merged

  return merged

def regionOverlap(r1, r2, isLoop= False, courseDistance= None):
  """return the overlap of two regiom, where regions are specified as a pair of positions"""
  if isLoop:
    if courseDistance is None:
      stderr.write(f'ERROR: region overlap called with isLoop = True but no courseDistance\n')
      sys.exit()
    # map coordinates to not wrap around
    r1[0] -= courseDistance * floor(r1[0] / courseDistance)
    r1[1] -= courseDistance * floor((r1[1] - r1[0]) / courseDistance)
    r2[0] -= courseDistance * floor(r2[0] / courseDistance)
    r2[1] -= courseDistance * floor((r2[1] - r2[0])/ courseDistance)
    
  start = max(r1[0], r2[0])
  stop = min(r1[1], r2[1])

  if isLoop:
    stop -= courseDistance * floor(stop / courseDistance)
    if start == stop:
      return []
    else:
      return [start, stop]
  else:
    if stop > start:
      return [start, stop]
    else:
      return None

def segmentIntercept(s12, s34):
  p1, p2 = s12
  p3, p4 = s34
  x1, y1 = 0, 0
  x2, y2 = pointdxdy(p1, p2)
  x3, y3 = pointdxdy(p1, p3)
  x4, y4 = pointdxdy(p1, p4)
  dx12, dy12 = x2, y2
  dx34, dy34 = pointdxdy(p3, p4)
  denom = dx34 * dy12 - dx12 * dy34
  a = ( (dx12 ** 2 + dy12 ** 2) * (dx34 ** 2 + dy34 ** 2) ) ** (1/2)
  if a == 0 or abs(denom) < 0.01 * a:
    return []
  f12 = (dx34 * (y3 - y1) - dy34 * (x3 - x1)) / denom

  if (f12 >= 0) and (f12 < 1):
    x = f12 * x2 + (1 - f12) * x1
    y = f12 * y2 + (1 - f12) * y1
    if abs (x3 - x4) > abs(y3 - y4):
      f23 =  (x - x3) / (x4 - x3)
    else:
      f23 = (y - y3) / (y4 - y3)
    if (f23 >= 0) and (f23 < 1):
      return f12, f23
  return []
    
def applySelectiveSpacing(points, spacing, ranges = None, isLoop = False, distances = None, courseDistance = None):
  """interpolate points to get a specified max point spacing from d1 to d2: wrap around for loops.
  It returns a list of points, a list of distances, and a map of indices to the original points"""

  # create distance field
  if distances is None or courseDistance is None:
    distances, courseDistance = calcDistances(points= points, isLoop = isLoop)

  N = len(points)
  
  if spacing is None or spacing <= 0 or spacing > courseDistance or N == 0:
    return points, distances, range(N)
    
  # merge ranges
  if ranges is not None:
    if len(ranges) == 0:
      return points, distances, range(N)
    ranges = mergeRanges(points= points, ranges= ranges, isLoop= isLoop, courseDistance = courseDistance)
  
  newPoints = []
  newDistances = []
  indexMap = []
  nRange = 0
  for i in range(N):
    # interval on points
    p1 = points[i]
    d1 = distances[i]
    j = (i + 1) % N
    p2 = points[j]
    d2 = distances[j]
    d12 = deltaPosition(d1, d2, isLoop = isLoop, courseDistance = courseDistance)
    newPoints.append(p1)
    newDistances.append(d1)
    indexMap.append(i)

    if not isLoop and i == N - 1:
      break

    # increment range or bail if no ranges remain
    if ranges is None:
      if isLoop or j > N:
        nPoints = int(d12 / spacing)
        for n in range(1, nPoints + 1):
          f = n / (nPoints + 1)
          d = d1 + f * d12
          newPoints.append(interpolatePoint(points[i], points[j], f))
          newDistances.append(d)
          indexMap.append(None)
    else:
      while ranges[nRange][1] < d1:
        nRange += 1
        if nRange >= len(ranges):
          newPoints += points[i + 1:]
          newDistances += distances[i + 1:]
          indexMap += range(i + 1, N)
          return newPoints, newDistances, indexMap
    
      # check for region overlap
      while nRange < len(ranges):
        overlap = regionOverlap([d1, d2], ranges[nRange], isLoop= isLoop, courseDistance= courseDistance)
        if not overlap:
          break
        # there is overlap, so add the points
        start, stop = overlap
        addPoints = []
        addDistances = []
        addIndices = []
        if abs(start - d1) < spacing / 3:
          start = d1
        if abs(start - d2) < spacing / 3:
          stop = d2
        l = deltaPosition(start, stop, isLoop= isLoop, courseDistance= courseDistance)
        nPoints = int(l / spacing)
        for n in range(nPoints + 1):
          if nPoints == 0:
            d = (start + stop) / 2
            f = 0.5
          else:
            d = start + n * l / nPoints
            f = (d - d1) / d12
          if d - d1 > spacing / 3 and (d2 - d) > spacing / 3:
            addPoints.append(interpolatePoint(p1, p2, f))
            addDistances.append(d)
            addIndices.append(None)

        newPoints += addPoints
        newDistances += addDistances
        indexMap += addIndices
        # if the range is in this segment, then jump to the next range
        if d2 > stop:
          nRange += 1
          if ranges is not None and nRange >= len(ranges):
            newPoints += points[i + 1:]
            newDistances += distances[i + 1:]
            indexMap += range(i + 1, N)
            return newPoints, newDistances, indexMap
        else:
          break
  return newPoints, newDistances, indexMap
  
def smoothPoints(points, smoothXY=0, smoothZ=0, isLoop=False, distances = None, courseDistance = None, refineMesh = True):
  '''smooth points with the given sigma values.  Loop courses wrap around.'''
  if smoothXY is None or smoothXY < 0:
    smoothXY = 0
  if smoothZ is None or smoothZ < 0:
    smoothZ = 0
  if smoothXY <= 0 and smoothZ <= 0:
    return points

  wList = [1]
  for u in range(1, 41):
    wList.append(exp(- u ** 2 / 200))

  N = len(points)

  # create distance field
  if distances is None or courseDistance is None:
    distances, courseDistance = calcDistances(points= points, isLoop = isLoop)

  #
  # refine the mesh: add mesh points proximate to points
  #

  if refineMesh:

    # do refinements in order of refinement
    sigmas = []
    if smoothXY > 0:
      sigmas.append(smoothXY)
    if smoothZ > 0:
      sigmas.append(smoothZ)
    if len(sigmas) > 1:
      if sigmas[0] < sigmas[1]:
        s = sigmas[0]
        sigmas[0] = sigmas[1]
        sigmas[1] = s
      # if the sigmas are similar, just use the broader sigma
      if sigmas[0] < 1.2 * sigmas[1]:
        sigmas = [sigmas[0]]

    newPoints = [[[], []] for p in points]
    newDistances = [[[], []] for p in points]

    # create refinement region in the zones around each point
    regions = []    # list: for each sigma, for each point, 2 refinement region, each with 2 bounds
    ns = 0
    for sigma in sigmas:
      regions.append([])
      r = 2 * sigma
      for i in range(N):
        d = distances[i]
        d1 = d - r
        if not isLoop and d1 < 0:
          d1 = 0
        if isLoop or i > 0:
          dMid = max(0, interpolatePosition(distances[i - 1], distances[i], 0.5, courseDistance, isLoop))   # midpoint between points
          if d1 < dMid:
            d1 = dMid
        d2 = d + r
        if not isLoop and d1 > courseDistance:
          d1 = courseDistance
        if isLoop or i < N - 1:
         dMid = interpolatePosition(distances[i], distances[(i + 1) % N], 0.5, courseDistance, isLoop)   # midpoint between points
         if d2 > dMid:
           d2 = dMid
        # refinement region for this spacing is from d1 to d2, retaining original point
        regions[ns].append([[d1, d], [d, d2]])

        # if there's an additional region, it gets split
        for n in range(ns):
          regions[n][i][0][1] = d1   # end of first region
          regions[n][i][1][0] = d2   # start of 2nd region
      ns += 1

    nNewPoints = 0
    for i in range(N):
      j = (i - 1) % N
      k = (i + 1) % N
      # sigma regions prior to point
      for ns in range(len(sigmas)):
        spacing = sigmas[ns] / 2
        d1 = regions[ns][i][0][0]
        d2 = regions[ns][i][0][1]
        NPoints = 1 + int((d2 - d1) / spacing)

        # avoid close points at the boundary
        n1 = 0
        dMid = max(0, interpolatePosition(distances[i], distances[j], 0.5, courseDistance, isLoop))   # midpoint between points
        if abs(d1 - dMid) < spacing / 4:
          n1 = 1

        dd = deltaPosition(distances[j], distances[i], courseDistance = courseDistance, isLoop = isLoop)
        for n in range(n1, NPoints):
          d = interpolatePosition(d1, d2, f= n / NPoints, courseDistance = courseDistance)
          f = deltaPosition(distances[j], d, courseDistance = courseDistance, isLoop = isLoop) / dd
          newPoints[i][0].append(interpolatePoint(points[j], points[i], f))
          newDistances[i][0].append(d)
          nNewPoints += 1
          
      # sigma regions after point
      for ns in reversed(range(len(sigmas))):
        spacing = sigmas[ns] / 2
        d1 = regions[ns][i][1][0]
        d2 = regions[ns][i][1][1]
        NPoints = 1 + int((d2 - d1) / spacing)
        dd = deltaPosition(distances[i], distances[k], courseDistance = courseDistance, isLoop = isLoop)
        for n in range(1, NPoints):
          d = interpolatePosition(d1, d2, f= n / NPoints, courseDistance = courseDistance)
          f = deltaPosition(distances[i], d, courseDistance = courseDistance, isLoop = isLoop) / dd
          newPoints[i][1].append(interpolatePoint(points[i], points[k], f))
          newDistances[i][1].append(d)
          nNewPoints += 1

        # revise possible duplicate points at the boundary
        if isLoop or i < N - 1:
          dMid = interpolatePosition(distances[i], distances[j] , 0.5, courseDistance, isLoop)   # midpoint between points
          # if we're close to the mid-point, use that instead
          if len(newPoints[i][1]) > 0 and abs(d2 - dMid) < spacing / 4:
            newPoints[i][1][-1] = interpolatePoint(points[i], points[j], 0.5)
            newDistances[i][1][-1] = dMid
        # end of point to point, don't put a point too close to the final point
        if not isLoop and i == N - 1 and abs(d2 - courseDistance) < spacing / 4 and len(newPoints[i][1]) > 0:
          newPoints[i][1].pop()
          newDistances[i][1].pop()
          nNewPoints -= 1

    # if we added new points, splice them in
    # before the first point goes at the end
    if nNewPoints > 0:
      pRefined = []
      dRefined = []
      for i in range(N):
        if i > 0:
          pRefined += newPoints[i][0]
          dRefined += newDistances[i][0]
        pRefined += [points[i]] + newPoints[i][1]
        dRefined += [distances[i]] + newDistances[i][1]
      if isLoop:
        pRefined += newPoints[0][0]
        dRefined += newDistances[0][0]
      points = pRefined
      distances = dRefined
      N = len(points)


  #
  # calculate measures of points
  #
  ms = []
  for i in range(len(distances)):
    j = i - 1 if i > 0 or isLoop else 0
    k = (i + 1) % N if i < N - 1 or isLoop else N - 1
    ms.append(deltaPosition(distances[j], distances[k], courseDistance = courseDistance, isLoop = isLoop))

  #
  # do smoothing on refined mesh
  #

  newLats = []
  newLons = []
  newEles = []
  if smoothXY > 0:
    uXYs = [d / smoothXY for d in distances]
    mXYs = []
    cdXY = courseDistance / smoothXY
    j = 0
    k = 0
    if isLoop:
      while ((j - 1) % N != i % N) and deltaPosition(uXYs[j], uXYs[i], cdXY, isLoop) < 3:
        j = (j - 1) % N

    lons = [points[0].longitude]   # avoid antimeridian problems
    for i in range(N - 1):
      dLon = points[i + 1].longitude - points[i].longitude
      dLon -= 360 * floor(0.5 + dLon / 360)
      lons.append(lons[-1] + dLon)

    for i in range(N):
      if isLoop:
        iLast = (i - 1) % N
      else:
        iLast = N - 1
      while (k % N != iLast) and deltaPosition(uXYs[i], uXYs[k], cdXY, isLoop) < 3:
        k = (k + 1) % N
      while (j % N != i) and deltaPosition(uXYs[(j + 1) % N], uXYs[i], cdXY, isLoop) > 3:
        j = (j + 1) % N
      # smooth from k to j
      u = j
      sumw = 0
      sumLat = 0
      sumLon = 0
      while True:
        d = deltaPosition(uXYs[u], uXYs[i], cdXY, isLoop)
        nw = int(10 * abs(d) + 0.5)
        w = ms[u] * wList[int(10 * abs(d) + 0.5)] if nw < len(wList) else 0
        sumw += w
        sumLat += w * points[u].latitude
        sumLon += w * lons[u]
        if u == k:
          break
        u = (u + 1) % N

      lat = sumLat / sumw
      lon = sumLon / sumw
      lon -= 360 * floor(0.5 + lon / 360)
      newLats.append(lat)
      newLons.append(lon)

  if (smoothZ > 0):
    uZs = [d / smoothZ for d in distances]
    cdZ = courseDistance / smoothZ
    j = 0
    k = 0
    i = 0
    if isLoop:
      while ((j - 1) % N != i % N) and deltaPosition(uZs[j], uZs[i], cdZ, isLoop) < 3:
        j = (j - 1) % N
    lons = [points[0].longitude]   # avoid antimeridian problems
    for i in range(N - 1):
      dLon = points[i + 1].longitude - points[i].longitude
      dLon -= 360 * floor(0.5 + dLon / 360)
      lons.append(lons[-1] + dLon)

    for i in range(N):
      if isLoop:
        iLast = (i - 1) % N
      else:
        iLast = N - 1
      while (k % N != iLast) and deltaPosition(uZs[i], uZs[k], cdZ, isLoop) < 3:
        k = (k + 1) % N
      while (j % N != i) and deltaPosition(uZs[(j + 1) % N], uZs[i], cdZ, isLoop) > 3:
        j = (j + 1) % N
      # smooth from k to j
      u = j
      sumw = 0
      sumEle = 0
      while True:
        d = deltaPosition(uZs[u], uZs[i], cdZ, isLoop)
        nw = int(10 * abs(d) + 0.5)
        w = wList[int(10 * abs(d) + 0.5)] if nw < len(wList) else 0
        sumw += w
        sumEle += w * points[u].elevation
        if u == k:
          break
        u = (u + 1) % N

      ele = sumEle / sumw
      newEles.append(ele)

  if len(newLons) > 0:
    for p, lat, lon  in zip(points, newLats, newLons):
      p.latitude = lat
      p.longitude = lon

  if len(newEles) > 0:
    for p, ele  in zip(points, newEles):
      p.elevation = ele

  return points

def calcCurvatures(points = [], isLoop = False, maxR = 10, distances = None, courseDistance = None):
  N = len(points)
  cs = []
  if courseDistance is None or distances is None:
    distances, courseDistance = calcDistances(points= points, isLoop= isLoop)
  for i in range(N):
    # get 3 points
    j = (i - 1) % N
    while j != (i + 1) % N and pointsAreVeryClose(points[j], points[i]):
      j = (j - 1) % N
    k = (i + 1) % N
    while k != (i - 1) % N and pointsAreVeryClose(points[k], points[i]):
      k = (k + 1) % N
    if not isLoop and (i == 0 or i == N - 1):
      cs.append(0)
    else:
      d1 = deltaPosition(distances[j], distances[i], courseDistance = courseDistance, isLoop = isLoop)
      d2 = deltaPosition(distances[i], distances[k], courseDistance = courseDistance, isLoop = isLoop)
      if maxR is not None:
        d1 = min(maxR, d1)
        d2 = min(maxR, d2)
      angle = pointAngle(points[j], points[i], points[k])
      c = 2 * angle / (d1 + d2)
      cs.append(c)
  return cs
      
def applyMinRadius(points = [], minRadius = None, isLoop = False, distances = None, courseDistance = None, kTransition= 0.2, lTransition= 5):
  """shift points to create a minimum radius in corners"""
  if len(points) < 3 or minRadius is None or minRadius <= 0:
    return points

  # optionally apply min radius in steps
  # steps = ceil(10 * sqrt(minRadius / 10))
  steps = 1

  for nStep in range(1, steps + 1):
    N = len(points)
    minR = minRadius * nStep / steps
  
    if courseDistance is None or distances is None:
      distances, courseDistance = calcDistances(points= points, isLoop= isLoop)

    curvatures = calcCurvatures(points, isLoop = isLoop, distances= distances, courseDistance = courseDistance)

    lShifts = []
    rShifts = []
    for i in range(N):
      c = curvatures[i]
      score = minR * c
      if score > 1:
        rShifts.append(minR - 1 / c)
      else:
        rShifts.append(0)
      if score < -1:
        lShifts.append(minR + 1 / c)
      else:
        lShifts.append(0)

    # refine the points around regions of shift
    ranges = []
    r = []
    for i in range(N):
      if rShifts[i] > 0 or lShifts[i] > 0:
        if r:
          r[-1] = distances[i]
        else:
          r = [distances[i], distances[i]]
      else:
        if r:
          r = [r[0] - 3 * lTransition, r[1] + 3 * lTransition]
          ranges.append(r)
          r = []

    if ranges:
      spacing = 1 + lTransition / 3
      points, distances, indexMap = applySelectiveSpacing(points = points, ranges = ranges, spacing = spacing, isLoop = isLoop, distances = distances, courseDistance = courseDistance)
      N = len(points)
      # adjust parameters to new points
      newLShifts = []
      newRShifts = []
      newCurvatures = []
      for i in range(len(indexMap)):
        if indexMap[i] is None:
          newLShifts.append(0)
          newRShifts.append(0)
          newCurvatures.append(0)
        else:
          newLShifts.append(lShifts[indexMap[i]])
          newRShifts.append(rShifts[indexMap[i]])
          newCurvatures.append(curvatures[indexMap[i]])
      lShifts = newLShifts
      rShifts = newRShifts
      curvatures = newCurvatures

    # apply shift transition forward
    for i in range(N):
      if isLoop or i > 0:
        j = i - 1
        dd = deltaPosition(distances[j], distances[i], courseDistance = courseDistance, isLoop = isLoop)
        rShiftMin = max(0, max(rShifts[j] - kTransition * dd, rShifts[j] * exp(-dd / lTransition)))
        if rShifts[i] < rShiftMin:
          rShifts[i] = rShiftMin
        lShiftMin = max(0, max(lShifts[j] - kTransition * dd, lShifts[j] * exp(-dd / lTransition)))
        if lShifts[i] < lShiftMin:
          lShifts[i] = lShiftMin
    for i in reversed(range(N)):
      if isLoop or i < N - 1:
        j = (i + 1) % N
        dd = deltaPosition(distances[i], distances[j], courseDistance = courseDistance, isLoop = isLoop)
        rShiftMin = max(0, max(rShifts[j] - kTransition * dd, rShifts[j] * exp(-dd / lTransition)))
        if rShifts[i] < rShiftMin:
          rShifts[i] = rShiftMin
        lShiftMin = max(0, max(lShifts[j] - kTransition * dd, lShifts[j] * exp(-dd / lTransition)))
        if lShifts[i] < lShiftMin:
          lShifts[i] = lShiftMin
    shifts = [rShifts[i] - lShifts[i] for i in range(N)]

    points = applyShift(points= points, shifts = shifts, curvatures = curvatures, distances = distances, courseDistance = courseDistance, isLoop = isLoop)

    distances = None
    courseDistance = None
    curvatures = None
    
  return points


def shiftPoint(point, direction, distance):
  """shift a point by a given distance in a given direction"""
  c = cos(direction)
  s = sin(direction)

  # lane shift, 90 degrees
  dx = s * distance
  dy = -c * distance
  dlng = dx / (lat2y * cos(deg2rad * point.latitude))
  dlat = dy / lat2y
  pNew = copyPoint(point)
  pNew.longitude += dlng
  pNew.latitude  += dlat
  return pNew


def shiftVertex (point, directions, distance):
  """ shift a vertex by a given distance in a given direction
  the vertex is the intercection of two lines, each of
  which are shifted"""

  dir1, dir2 = directions

  c1 = cos(dir1)
  s1 = sin(dir1)
  c2 = cos(dir2)
  s2 = sin(dir2)

  # lane shift, 90 degrees
  denom = c1 * s2 - c2 * s1
  if abs(denom) < 0.001:
    dx =  (s1 + s2) * distance / 2
    dy = -(c1 + c2) * distance / 2
  else:
    dx = (c1 - c2) / denom * distance
    dy = (s1 - s2) / denom * distance

  dlng = dx / (lat2y * cos(deg2rad * point.latitude))
  dlat = dy / lat2y
  pNew = copyPoint(point)
  pNew.longitude += dlng
  pNew.latitude  += dlat
  return pNew

def applyShift(points = [], shifts = None, curvatures = None, distances = None, courseDistance = None, isLoop = False):
  if shifts is None or len(shifts) == 0:
    return points

  if curvatures is None:
    curvatures = calcCurvatures(points, isLoop = isLoop, distances= distances, courseDistance = courseDistance)

  if distances is None or courseDistance is None:
    distances, courseDistance = calcDistances(points= points, isLoop= isLoop)

  N = len(points)

  newPoints = []

  # list of point vectors and directions
  dirs = []
  for i in range(N):
    j = (i + 1) % N
    while (pointsAreClose(points[i], points[j]) and i != (j + 1) % N):
      j = (j + 1) % N
    if i > j and not isLoop:
      dirs.append(dirs[-1])
    else:
      dx, dy = pointdxdy(points[i], points[j])
      dirs.append(atan2(dy, dx))

  # lane shift: to right, which means adding pi/2 to the direction
  for i in range(N):
    dir1 = dirs[i - 1] if i > 0 or isLoop else dirs[0]
    dir2 = dirs[i]

    # check for excessive shift: R = 1/curvature, shift > R, shift x curvature > 1
    t = curvatures[i] * shifts[i]
    if t < -0.75:
      shifts[i] = -0.75 / curvatures[i]
      # check to make sure lane shift isn't changing too rapidly to nearby points
      for direction in (-1, 1):
        j = (i + direction) % N
        while j != i and (isLoop or j > i):
          d = deltaPosition(distances[j], distances[i], courseDistance = courseDistance, isLoop = isLoop)
          dShift = shifts[j] - shifts[i]
          if (dShift > 0 and shifts[i] > 0) or (dShift < 0 and shifts[i] < 0):
            shifts[j] = shifts[i]
          elif 2 * abs(dShift) > abs(d):
            shifts[j] = shifts[i] + (d if dShift > 0 else -d) / 2
          else:
            break
          j = (j + direction) % N
    # for sharp turns repeat a point: there's no way to decide if it's an "inside" or "outside" sharp turn
    if (isLoop or (i > 0 and i < N - 1)) and abs(dir2 - dir1) > 0.99 * pi:
      pTurns = []
      for dir in [dir1, dir2]:
        pTurns.append(shiftPoint(point= points[i], direction= dir, distance = shifts[i]))

      # check if there's a knot.. if not use the doubled points
      fs = segmentIntercept([points[i - 1], pTurns[0]], [pTurns[1], points[(i + 1) % N]])
      if len(fs) == 0:
        newPoints += pTurns
        next

    newPoints.append(shiftVertex(point = points[i], directions = [dir1, dir2], distance = shifts[i]))

  return newPoints

