from math import pi, cos, sin, sqrt, atan2, floor
import gpxpy
import gpxpy.gpx
import sys
rEarth = 20037392 / pi


twopi = 2 * pi
deg2rad = pi / 180

def load_gpx(gpx_file):
  """
  Load GPX file and extract all track points.
  
  Args:
    filename: Path to the GPX file
    
  Returns:
    List of (latitude, longitude) tuples
  """

  gpx = gpxpy.parse(gpx_file)

  points = []
  for track in gpx.tracks:
    for segment in track.segments:
      for point in segment.points:
        points.append(point)
  
  sys.stderr.write(f"Loaded {len(points)} points from GPX file\n")
  return points


def copyPoint(p):
  return gpxpy.gpx.GPXTrackPoint(
    latitude= p.latitude,
    longitude= p.longitude,
    elevation= p.elevation
  )

def pointsAreClose(p1, p2, deltaDegs= 1e-6, deltaZ= 0.1):
  return ( \
    (p1.latitude is not None) and
    (p2.latitude is not None) and
    (p1.longitude is not None) and
    (p2.longitude is not None) and
    abs(p1.latitude - p2.latitude) <= deltaDegs and
    abs(p1.longitude - p2.longitude) <= deltaDegs and
    (
      (p1.elevation is None) or
      (p2.elevation is None) or
      (abs(p1.elevation - p2.elevation) <= deltaZ)
    )
  )

def pointsAreVeryClose(p1, p2, deltaDegs= 1e-8, deltaZ= 0.001):
  return pointsAreClose(p1, p2, deltaDegs = deltaDegs, deltaZ = deltaZ)

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
  d = sqrt( x1 ** 2 * x2 ** 2 + y1 ** 2 + y2 ** 2 + x1 ** 2 * y2 ** 2 + x2 ** 2 * y1 ** 2 )
  d = sqrt(dot ** 2 + cross ** 2)
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

def simplifyPoints(points, z0 = 0.2, r0 = 1):
  if r0 <= 0 or len(points) < 3:
    return points

  xf, yf = pointdxdy(points[0], points[-1])

  if xf ** 2 + yf ** 2 < 10:
    iFurthest = None
    dFurthest = 0
    for i in range(1, len(points)):
      d = pointDistance(points[0], points[i])
      if d > dFurthest:
        iFurthest = i
        dFurthest = d
    if iFurthest is not None:
      if iFurthest == len(points) - 1:
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
  for i in range(1, len(points)):
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
  # check that the route thru the S/F is not a U-turn
  i1 = -1;
  i2 = 0;
  i3 = 1;
  while ((pointDistance(points[i1], points[i2]) < 10) and (i1 > i2 - len(points) + 1)):
    i1 -= 1
  while ((pointDistance(points[i2], points[i3]) < 10) and (i3 < len(points) - 1)):
    i3 += 1
  isLoop = (
         (len(points) > 2) and
         (pointDistance(points[0], points[-1]) < 150) and
         (pointNormalizedDot(points[i1], points[i2], points[i2], points[i3]) > -0.1)
  )
  return isLoop;
