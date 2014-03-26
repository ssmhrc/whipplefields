import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from stellar_size import sample_sky

# call sample_sky directly
ra = np.arange(0, 360, 30)
dec = np.arange(-90, 90, 30)
pairs = list(product(ra, dec))
ra, dec = zip(*pairs)
counts = sample_sky(ra, dec)

def get_prev_run(subdir):
	ra, dec, counts = [], [], []
	for filename in os.listdir(subdir):
		if '.csv' in filename:
			df = pd.read_csv(subdir+'/'+filename)
			counts.append(df.shape[0])
			ra.append(int(filename.split('_')[0]))
			dec.append(int(filename.split('_')[1]))
	return ra, dec, counts

# get the input RA, Dec, and output counts from a previous run
subdir = 'sample_sky1'
ra, dec, counts = [np.array(i) for i in get_prev_run(subdir)]

# make a simple plot of the counts
plt.scatter(ra, dec, s=counts, alpha=0.5)
plt.xlabel('RA')
plt.ylabel('Dec')
for i in range(len(ra)):
	plt.text(ra[i]-8, dec[i]+10, str(counts[i]), fontsize=10)
plt.title('Candidate source number densities for Whipple')
plt.savefig('sample_sky1_bubble_plot.pdf')
plt.show()

# try a hexbin plot
plt.hexbin(ra, dec, counts, gridsize=5)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
plt.show()

# try some interpolation to make a heatmap
from scipy.interpolate import Rbf
xi = np.linspace(0,360,100)
yi = np.linspace(-90,90,100)
XI, YI = np.meshgrid(xi, yi)

rbfi = Rbf(ra, dec, counts, function='linear')
ZI = rbfi(XI, YI)
plt.pcolor(XI, YI, ZI)
plt.colorbar()
plt.show()

# from scipy.interpolate import NearestNDInterpolator
# interp = NearestNDInterpolator(np.array(zip(ra,dec)),np.array(counts))
# ZI = interp(XI, YI)
# plt.pcolor(XI, YI, ZI)
# plt.colorbar()
# plt.show()

# # try making a plot using an Aitoff map projection
# deg2rad = np.pi/180.
# xi = np.linspace(-np.pi, np.pi, 100)
# yi = np.linspace(-np.pi/2, np.pi/2, 100)
# XI, YI = np.meshgrid(xi, yi)
# rbfi = Rbf((ra - 180) * deg2rad, dec * deg2rad, counts, function='linear')
# ZI = rbfi(XI, YI)
# plt.subplot(111, projection="aitoff")
# plt.pcolor(XI, YI, ZI)
# plt.colorbar()
# ax = plt.gca()
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# plt.show()

# from scipy.interpolate import interp2d
# deg2rad = np.pi/180.
# z = counts.reshape(len(set(ra)), len(set(dec))).T
# intp = interp2d(np.arange(-180, 180, 30) * deg2rad, 
# 	np.arange(-90, 90, 30) * deg2rad, z)
# xi = np.linspace(-np.pi, np.pi, 100)
# yi = np.linspace(-np.pi/2, np.pi/2, 100)
# zi = intp(xi, yi)
# plt.subplot(111, projection="aitoff")
# plt.pcolor(zi)
# plt.colorbar()
# ax = plt.gca()
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# plt.show()

# intp = interp2d(ra , dec, counts) 
# xi = np.linspace(0, 360, 100)
# yi = np.linspace(-90, 90, 100)
# zi = intp(xi, yi)
# plt.pcolor(zi)
# plt.colorbar()
# plt.show()

# try linear interpolation over a regular grid
from scipy.interpolate import RectBivariateSpline
ra, dec, counts = [np.array(i) for i in get_prev_run(subdir)]
# convert to pandas dataframe to do nested 
df = pd.DataFrame(dict(ra=ra, dec=dec, counts=counts))
df_sorted = df.sort(['ra','dec'])
counts_sorted = np.array(df_sorted.counts.tolist())
x = np.arange(0, 360, 30)
y = np.arange(-90, 90, 30)
z = counts_sorted.reshape(x.size, y.size)
intp = RectBivariateSpline(x, y, z)
xi = np.linspace(0, 360, 200)
yi = np.linspace(-90, 90, 100)
zi = intp(xi, yi)
# add one more element to the x and y arrays to make pcolor happy
xil, yil = xi.tolist(), yi.tolist()
xil.append(2*xi[-1]-xi[-2])
yil.append(2*yi[-1]-yi[-2])
xia, yia = np.array(xil), np.array(yil)
# plt.pcolor(xia, yia, zi.T)
# plt.pcolormesh(xia, yia, zi.T)
plt.pcolormesh(zi.T)
xloc = np.arange(0, 201, 25)
xlab = map(lambda x: str(x), np.arange(0, 361, 45))
yloc = np.arange(0, 101, 10)
ylab = map(lambda x: str(x), np.arange(-90, 91, 18))
plt.xticks(xloc, xlab)
plt.yticks(yloc, ylab)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
plt.show()
# # try aitoff projection
# plt.subplot(111, projection="hammer")
# plt.pcolor((xia-180)*deg2rad, yia*deg2rad, zi.T)
# plt.pcolormesh((xia-180)*deg2rad, yia*deg2rad, zi.T)
# plt.colorbar()
# ax = plt.gca()
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# plt.show()
# # not working... only the lower left quadrant seems to come out intact




# =====================================================

# now repeat it with a finer grid of coordinates
ra = np.arange(0, 360, 5)
dec = np.arange(-90, 90, 5)
pairs = list(product(ra, dec))
ra, dec = zip(*pairs)
counts2 = sample_sky(ra, dec)

subdir = 'sample_sky2'
ra, dec, counts = [np.array(i) for i in get_prev_run(subdir)]

# make a hexbin plot
# -----------------------------------------------------
plt.figure(figsize=(9,5.5))
plt.hexbin(ra,dec,counts,gridsize=30)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
plt.savefig('sample_sky2_hexbin.pdf')
plt.show()
# -----------------------------------------------------

# make a linearly interpolated heatmap
# -----------------------------------------------------
# missing data at index 1755: insert missing RA/Dec and a nan in counts
ra, dec, counts = [i.tolist() for i in [ra, dec, counts]]
idx = 1755
ra.insert(idx,240)
dec.insert(idx,45)
counts.insert(idx,np.nan)
ra, dec, counts = [np.array(i) for i in [ra, dec, counts]]

df = pd.DataFrame(dict(ra=ra, dec=dec, counts=counts))
df_sorted = df.sort(['ra','dec'])
counts_sorted = np.array(df_sorted.counts.tolist())

x = np.arange(0, 360, 5)
y = np.arange(-90, 90, 5)
z = counts_sorted.reshape(x.size, y.size)

idx = np.argwhere(np.isnan(z)).flatten()
idxx, idxy = idx[0], idx[1]
# set the nan value equal to the average of its 8 neighbors
z[idxx, idxy] = np.nanmean(z[idxx-1:idxx+2,idxy-1:idxy+2])

intp = RectBivariateSpline(x, y, z)
xi = np.linspace(0, 360, 200)
yi = np.linspace(-90, 90, 100)
zi = intp(xi, yi)
# add one more element to the x and y arrays to make pcolor happy
xil, yil = xi.tolist(), yi.tolist()
xil.append(2*xi[-1]-xi[-2])
yil.append(2*yi[-1]-yi[-2])
xia, yia = np.array(xil), np.array(yil)
# plt.pcolor(xia, yia, zi.T)
# plt.pcolormesh(xia, yia, zi.T)
plt.figure(figsize=(9,5.5))
# plt.pcolor(zi.T)
plt.pcolormesh(zi.T)
xloc = np.arange(0, 201, 25)
# xloc = np.arange(0, 201, 25) * 2
xlab = map(lambda x: str(x), np.arange(0, 361, 45))
yloc = np.arange(0, 101, 10)
# yloc = np.arange(0, 101, 10) * 2
ylab = map(lambda x: str(x), np.arange(-90, 91, 18))
plt.xticks(xloc, xlab)
plt.yticks(yloc, ylab)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
plt.savefig('sample_sky2_linear_interp.pdf')
plt.show()
# -----------------------------------------------------


# # make a radial basis function smoothed heatmap
# # -----------------------------------------------------
# subdir = 'sample_sky2'
# ra, dec, counts = [np.array(i) for i in get_prev_run(subdir)]
# # missing data at index 1755: insert missing RA/Dec and a nan in counts
# ra, dec, counts = [i.tolist() for i in [ra, dec, counts]]
# idx = 1755
# ra.insert(idx,240)
# dec.insert(idx,45)
# counts.insert(idx,95.625)
# ra, dec, counts = [np.array(i) for i in [ra, dec, counts]]

# xi = np.linspace(0,360,200)
# yi = np.linspace(-90,90,100)
# XI, YI = np.meshgrid(xi, yi)

# rbfi = Rbf(ra, dec, counts, function='gaussian')
# ZI = rbfi(XI, YI)
# plt.pcolormesh(XI, YI, ZI)
# plt.colorbar()
# plt.show()
# # not working :(

# from scipy import interpolate
# df = pd.DataFrame(dict(ra=ra, dec=dec, counts=counts))
# df_sorted = df.sort(['ra','dec'])
# counts_sorted = np.array(df_sorted.counts.tolist())

# x, y = np.mgrid[0:360:72j,-90:90:36j]
# z = counts_sorted.reshape(72, 36)
# tck = interpolate.bisplrep(x, y, z, s=0.1)
# xi, yi = np.mgrid[0:360:288j,-90:90:72j]
# zi = interpolate.bisplev(xi[:,0],yi[0,:],tck)
# plt.pcolormesh(xi, yi, zi)
# plt.show()
# # not working :(
# # -----------------------------------------------------


# convert to galactic coordinates
# -----------------------------------------------------
subdir = 'sample_sky2'
ra, dec, counts = [np.array(i) for i in get_prev_run(subdir)]
ra, dec, counts = [i.tolist() for i in [ra, dec, counts]]
idx = 1755
ra.insert(idx,240)
dec.insert(idx,45)
counts.insert(idx,95.625)
ra, dec, counts = [np.array(i) for i in [ra, dec, counts]]

from astropy import coordinates as coord
from astropy import units as u
c = coord.ICRS(ra=ra, dec=dec, unit=(u.degree, u.degree))
galactic = c.transform_to(coord.GalacticCoordinates)
rad2deg = 180./np.pi
plt.figure(figsize=(9,5.5))
plt.hexbin(galactic.l * rad2deg, galactic.b * rad2deg, counts, gridsize=20)
plt.colorbar()
plt.xlabel(r'$l$',fontsize=18)
plt.ylabel(r'$b$',fontsize=18)
plt.title('Viable sources: counts per square degree')
plt.savefig('sample_sky2_hexbin_galactic.pdf')
plt.show()

# # now linearly interpolate
# df = pd.DataFrame(dict( l = 5*((galactic.l.value*rad2deg)/5).astype('int'), 
# 						b = 5*((galactic.b.value*rad2deg)/5).astype('int'), 
# 						counts = counts) )
# df_sorted = df.sort(['l','b'])
# counts_sorted = np.array(df_sorted.counts)

# x = np.arange(0, 360, 5)
# y = np.arange(-90, 90, 5)
# z = counts_sorted.reshape(x.size, y.size)

# intp = RectBivariateSpline(x, y, z)
# xi = np.linspace(0, 360, 200)
# yi = np.linspace(-90, 90, 100)
# zi = intp(xi, yi)
# # add one more element to the x and y arrays to make pcolor happy
# xil, yil = xi.tolist(), yi.tolist()
# xil.append(2*xi[-1]-xi[-2])
# yil.append(2*yi[-1]-yi[-2])
# xia, yia = np.array(xil), np.array(yil)
# # plt.pcolor(xia, yia, zi.T)
# # plt.pcolormesh(xia, yia, zi.T)
# plt.figure(figsize=(9,5.5))
# # plt.pcolor(zi.T)
# plt.pcolormesh(zi.T)
# xloc = np.arange(0, 201, 25)
# # xloc = np.arange(0, 201, 25) * 2
# xlab = map(lambda x: str(x), np.arange(0, 361, 45))
# yloc = np.arange(0, 101, 10)
# # yloc = np.arange(0, 101, 10) * 2
# ylab = map(lambda x: str(x), np.arange(-90, 91, 18))
# plt.xticks(xloc, xlab)
# plt.yticks(yloc, ylab)
# plt.colorbar()
# plt.xlabel(r'$l$', fontsize=18)
# plt.ylabel(r'$b$', fontsize=18)
# plt.title('Viable sources: counts per square degree')
# # plt.savefig('sample_sky2_linear_interp.pdf')
# plt.show()
# # not working :(


# convert to ecliptic coordinates
# -----------------------------------------------------

# def eq_to_ec(ra, dec):
# 	"""
# 	Converts equatorial coordinates to ecliptic, with ra, dec in degrees.
# 	"""
# 	from astropy.coordinates import spherical_to_cartesian
# 	from astropy.coordinates import cartesian_to_spherical
# 	deg2rad = np.pi/180.
# 	rad2deg = 180./np.pi
# 	# convert to cartesian
# 	x_eq, y_eq, z_eq = spherical_to_cartesian(1, ra * deg2rad, dec * deg2rad)
# 	# define obliquity of the ecliptic
# 	epsilon = 23.4 * deg2rad
# 	# define rotation matrix
# 	cos_eps = np.cos(epsilon)
# 	sin_eps = np.sin(epsilon)
# 	rotation = np.array([	[1, 0, 0], 
# 							[0, cos_eps, sin_eps], 
# 							[0, -sin_eps, cos_eps] ])
# 	c_eq = np.c_[x_eq, y_eq, z_eq].T
# 	x_ec, y_ec, z_ec = np.dot(rotation, c_eq)
# 	r, lam, beta = cartesian_to_spherical(x_ec, y_ec, z_ec)
# 	return lam * rad2deg, beta * rad2deg

# lam, beta = eq_to_ec(ra, dec)
# beta, lam = eq_to_ec(ra, dec)

from ephem import Equatorial, Ecliptic
lam, beta = [], []
for i, j in zip(ra, dec):
	eq = Equatorial(i * deg2rad, j * deg2rad)
	ec = Ecliptic(eq)
	lam.append(ec.lon * rad2deg)
	beta.append(ec.lat * rad2deg)

plt.figure(figsize=(9,5.5))
plt.hexbin(lam, beta, counts, gridsize=20)
plt.colorbar()
plt.xlabel(r'$\lambda$',fontsize=18)
plt.ylabel(r'$\beta$',fontsize=18)
plt.title('Viable sources: counts per square degree')
plt.savefig('sample_sky2_hexbin_ecliptic.pdf')
plt.show()


# try new approach to plotting using plt.scatter
# equatorial
plt.scatter(ra, dec, c=counts, marker='s', s=100, linewidths=0, alpha=0.5)
plt.show()
# galactic
l = galactic.l.value * rad2deg
b = galactic.b.value * rad2deg
plt.scatter(l, b, c=counts, marker='s', s=100, linewidths=0, alpha=0.5)
plt.show()
# ecliptic
plt.scatter(lam, beta, c=counts, marker='s', s=100, linewidths=0, alpha=0.5)
plt.show()

# now try the same as above but vary the marker size
s = np.zeros(counts.size)
for gal_lat in xrange(0,90,10):
	idx = np.argwhere( (gal_lat < np.abs(b)) & (np.abs(b) < gal_lat+10) )
	s[idx] = (gal_lat+50) * 1.25
plt.scatter(l, b, c=counts, marker='s', s=s, linewidths=0, alpha=0.5)
plt.show()

# would be better to vary the value for s by distance from the equatorial 
# poles in galactic coordinates
from ephem import Galactic
north_pole = Equatorial(0,90)
south_pole = Equatorial(0,-90)
l_np = Galactic(north_pole).lon * rad2deg
b_np = Galactic(north_pole).lat * rad2deg
l_sp = Galactic(south_pole).lon * rad2deg
b_sp = Galactic(south_pole).lat * rad2deg
s = np.zeros(counts.size)
def d_euclidian(x1, y1, x2, y2):
	return np.sqrt( (x1-x2)**2 + (y1-y2)**2 )
for i in range(len(s)):
	# set s[i] to the sum of the distance to each pole
	s[i] = ( d_euclidian(l[i], b[i], l_np, b_np) + \
		d_euclidian(l[i], b[i], l_sp, b_sp) )
plt.scatter(l, b, c=counts, marker='s', s=s/2, linewidths=0, alpha=0.5)
plt.colorbar()
plt.show()

# now try the same with Ecliptic
lam_np = Ecliptic(north_pole).lon * rad2deg
beta_np = Ecliptic(north_pole).lat * rad2deg
lam_sp = Ecliptic(south_pole).lon * rad2deg
beta_sp = Ecliptic(south_pole).lat * rad2deg
s = np.zeros(counts.size)
for i in range(len(s)):
	# set s[i] to the sum of the distance to each pole
	s[i] = ( d_euclidian(l[i], b[i], lam_np, beta_np) + \
		d_euclidian(l[i], b[i], lam_sp, beta_sp) )
plt.scatter(lam, beta, c=counts, marker='s', s=s/3, linewidths=0, alpha=0.5)
plt.colorbar()
plt.show()

# # try RBF on the converted data
# # not working, but I found the code below on a related post:
# # http://geoexamples.blogspot.com.es/2012/05/creating-grid-from-scattered-data-using.html
# def pointValue(x,y,power,smoothing,xv,yv,values):
#     nominator=0
#     denominator=0
#     for i in range(0,len(values)):
#         dist = np.sqrt((x-xv[i])*(x-xv[i])+(y-yv[i])*(y-yv[i])+smoothing*smoothing)
#         #If the point is really close to one of the data points, return the data point value to avoid singularities
#         if(dist<0.0000000001):
#             return values[i]
#         nominator=nominator+(values[i]/pow(dist,power))
#         denominator=denominator+(1/pow(dist,power))
#     #Return NODATA if the denominator is zero
#     if denominator > 0:
#         value = nominator/denominator
#     else:
#         value = -9999
#     return value

# def invDist(xv,yv,values,power=2,smoothing=0):
#     valuesGrid = np.zeros((len(yv),len(xv)))
#     for x in range(0,len(xv)):
#         for y in range(0,len(yv)):
#             valuesGrid[y][x] = pointValue(x,y,power,smoothing,xv,yv,values)
#     return valuesGrid

# xi = np.linspace(0,360,72)
# yi = np.linspace(-90,90,72)
# XI, YI = np.meshgrid(xi, yi)
# ZI = invDist(ra,dec,counts,power=1,smoothing=0)
# plt.pcolormesh(XI, YI, ZI)
# plt.show()
# # doesn't work 
# # runs slow because there are nested loops - needs vectorization


# try a narrow band around the plane of the ecliptic
# +/- 5 or 10 degrees from beta = 0 with a denser coordinate grid spacing
# start with coordinate grid in lambda, beta
lam = np.arange(0, 360, 2)
beta = np.arange(-10, 10, 1)
pairs = list(product(lam, beta))
# now convert from ecliptic to equatorial
ra, dec = [], []
for l, b in pairs:
	eq = Equatorial(Ecliptic(l, b))
	ra.append(eq.ra * rad2deg)
	dec.append(eq.dec * rad2deg)
counts = sample_sky(ra, dec)



# the previous run 'sample_sky2' was produced with a bug in cull_dataset.
# instead of correctly comparing the distance to the nearest neighbor of
# a source, it was comparing the neighbors to the wrong RA/Dec.
# recompute the counts from the VOTable xml files using the new corrected
# cull_dataset code
from stellar_size import cull_dataset
from astropy.io.votable import parse
vot_dir = 'sample_sky2'
outdir = vot_dir+'_fixed'
# os.mkdir(outdir)
vot_files = filter(lambda x: '.xml' in x, os.listdir(vot_dir))
# vot_paths = ['/'.join([vot_dir,i]) for i in vot_files]
ra, dec, counts = [], [], []
for v in vot_files:
	vot = parse(vot_dir+'/'+v)
	ra_i, dec_i = v.split('_')[:2]
	df = cull_dataset(outdir, ra_i, dec_i, vot.get_first_table())
	counts.append(df.shape[0])
	ra.append(float(ra_i))
	dec.append(float(dec_i))
	df.to_csv(outdir+'/'+v.replace('.xml','.csv'), index=False)

plt.figure(figsize=(9,5.5))
plt.hexbin(ra,dec,counts,gridsize=30)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
plt.savefig(outdir+'_hexbin.pdf')
plt.show()

deg2rad = np.pi/180.
rad2deg = 180./np.pi
from ephem import Equatorial, Ecliptic
lam, beta = [], []
for i, j in zip(ra, dec):
	eq = Equatorial(i * deg2rad, j * deg2rad)
	ec = Ecliptic(eq)
	lam.append(ec.lon * rad2deg)
	beta.append(ec.lat * rad2deg)

plt.figure(figsize=(9,5.5))
plt.hexbin(lam, beta, counts, gridsize=20)
plt.colorbar()
plt.xlabel(r'$\lambda$',fontsize=18)
plt.ylabel(r'$\beta$',fontsize=18)
plt.title('Viable sources: counts per square degree')
plt.savefig(outdir+'_hexbin_ecliptic.pdf')
plt.show()

# compare counts against previous (incorrect) results
subdir = 'sample_sky2'
ra_old, dec_old, counts_old = [np.array(i) for i in get_prev_run(subdir)]
df_old = pd.DataFrame(dict(ra=ra_old, dec=dec_old, counts=counts_old))
df_old_sorted = df_old.sort(['ra','dec'])
counts_old_sorted = np.array(df_old_sorted.counts.tolist())

df_new = pd.DataFrame(dict(ra=ra, dec=dec, counts=counts))
df_new_sorted = df_new.sort(['ra','dec'])
counts_new_sorted = np.array(df_new_sorted.counts.tolist())

counts_new_sorted.shape
# Out[1645]: 2591
counts_old_sorted.shape
# Out[1647]: (2591,)
(counts_new_sorted == counts_old_sorted).sum()/2591.
# Out[1646]: 0.99807024314936321

# so the fix to cull_dataset didn't change much at all -- about 0.2%


# got the new data back for the ecliptic plane
# -----------------------------------------------------
subdir = 'sample_sky_ecliptic_plane'

def get_prev_run(subdir):
	ra, dec, counts = [], [], []
	for filename in os.listdir(subdir):
		if '.csv' in filename:
			df = pd.read_csv(subdir+'/'+filename)
			counts.append(df.shape[0])
			ra.append(float(filename.split('_')[0]))
			dec.append(float(filename.split('_')[1]))
	return ra, dec, counts

ra, dec, counts = map(np.array,get_prev_run(subdir))

lam, beta = [], []
for i, j in zip(ra, dec):
	eq = Equatorial(i * deg2rad, j * deg2rad)
	ec = Ecliptic(eq)
	lam.append(ec.lon * rad2deg)
	beta.append(ec.lat * rad2deg)

# make a hexbin plot
plt.figure(figsize=(9,5.5))
plt.hexbin(ra,dec,counts,gridsize=30)
plt.hexbin(lam,beta,counts,gridsize=30)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
# plt.savefig('sample_sky2_hexbin.pdf')
plt.show()

# apparently the conversions using pyephem didn't work
# try using astropysics instead
from astropysics.coords import EclipticCoordinatesEquinox as Ecliptic
from astropysics.coords import EquatorialCoordinatesEquinox as Equatorial

lam = np.arange(0, 360, 2)
beta = np.arange(-10, 10, 1)
pairs = list(product(lam, beta))
# now convert from ecliptic to equatorial
ra, dec = [], []
for l, b in pairs:
	eq = Ecliptic(l, b).convert(Equatorial)
	ra.append(eq.ra.d)
	dec.append(eq.dec.d)

print pairs[0]
# (0, -10)
print Equatorial(ra[0],dec[0]).convert(Ecliptic)
# EclipticCoordinatesEquinox: lamb=0.000000,beta=-10.000000
counts = sample_sky(ra, dec)

# something is wrong with astropysics on ace7, cannot import

# found this implementation online at:
# http://www.astro.ucla.edu/~ianc/python/_modules/astrolib.html
def euler(ai, bi, select=1, fk4=False):
	"""
	NAME:
	    EULER
	PURPOSE:
	    Transform between Galactic, celestial, and ecliptic coordinates.
	EXPLANATION:
	    Use the procedure ASTRO to use this routine interactively

	CALLING SEQUENCE:
	     AO, BO = EULER(AI, BI, [SELECT=1, FK4=False])

	INPUTS:
	      AI - Input Longitude in DEGREES, scalar or vector.  If only two
	              parameters are supplied, then  AI and BI will be modified to
	              contain the output longitude and latitude.
	      BI - Input Latitude in DEGREES

	OPTIONAL INPUT:
	      SELECT - Integer (1-6) specifying type of coordinate transformation.

	     SELECT   From          To        |   SELECT      From            To
	      1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
	      2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
	      3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic

	     If not supplied as a parameter or keyword, then EULER will prompt for
	     the value of SELECT
	     Celestial coordinates (RA, Dec) should be given in equinox J2000
	     unless the /FK4 keyword is set.
	OUTPUTS:
	      AO - Output Longitude in DEGREES
	      BO - Output Latitude in DEGREES

	INPUT KEYWORD:
	      /FK4 - If this keyword is set and non-zero, then input and output
	            celestial and ecliptic coordinates should be given in equinox
	            B1950.
	      /SELECT  - The coordinate conversion integer (1-6) may alternatively be
	             specified as a keyword
	NOTES:
	      EULER was changed in December 1998 to use J2000 coordinates as the
	      default, ** and may be incompatible with earlier versions***.
	REVISION HISTORY:
	      Written W. Landsman,  February 1987
	      Adapted from Fortran by Daryl Yentis NRL
	      Converted to IDL V5.0   W. Landsman   September 1997
	      Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
	      Add option to specify SELECT as a keyword W. Landsman March 2003
	"""
	import numpy

	twopi = 2.0e0 * numpy.pi
	fourpi = 4.0e0 * numpy.pi

	#   J2000 coordinate conversions are based on the following constants
	#   (see the Hipparcos explanatory supplement).
	#  eps = 23.4392911111d              Obliquity of the ecliptic
	#  alphaG = 192.85948d               Right Ascension of Galactic North Pole
	#  deltaG = 27.12825d                Declination of Galactic North Pole
	#  lomega = 32.93192d                Galactic longitude of celestial equator
	#  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
	#  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
	#  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator

	if fk4:   
		equinox = '(B1950)'
		psi = numpy.array ([0.57595865315e0, 4.9261918136e0, 0.00000000000e0, 0.0000000000e0, 0.11129056012e0, 4.7005372834e0])
		stheta = numpy.array ([0.88781538514e0, -0.88781538514e0, 0.39788119938e0, -0.39788119938e0, 0.86766174755e0, -0.86766174755e0])
		ctheta = numpy.array([0.46019978478e0, 0.46019978478e0, 0.91743694670e0, 0.91743694670e0, 0.49715499774e0, 0.49715499774e0])
		phi = numpy.array([4.9261918136e0, 0.57595865315e0, 0.0000000000e0, 0.00000000000e0, 4.7005372834e0, 0.11129056012e0])
	else:   
		equinox = '(J2000)'
		psi = numpy.array([0.57477043300e0, 4.9368292465e0, 0.00000000000e0, 0.0000000000e0, 0.11142137093e0, 4.71279419371e0])
		stheta = numpy.array([0.88998808748e0, -0.88998808748e0, 0.39777715593e0, -0.39777715593e0, 0.86766622025e0, -0.86766622025e0])
		ctheta = numpy.array([0.45598377618e0, 0.45598377618e0, 0.91748206207e0, 0.91748206207e0, 0.49714719172e0, 0.49714719172e0])
		phi = numpy.array([4.9368292465e0, 0.57477043300e0, 0.0000000000e0, 0.00000000000e0, 4.71279419371e0, 0.11142137093e0])
	  
	i = select - 1
	a = numpy.deg2rad(ai) - phi[i]
	b = numpy.deg2rad(bi)
	sb = numpy.sin(b)
	cb = numpy.cos(b)
	cbsa = cb * numpy.sin(a)
	b = -stheta[i] * cbsa + ctheta[i] * sb
	bo = numpy.rad2deg(numpy.arcsin(numpy.minimum(b, 1.0)))
	del b
	a = numpy.arctan2(ctheta[i] * cbsa + stheta[i] * sb, cb * numpy.cos(a))
	del cb, cbsa, sb
	ao = numpy.rad2deg(((a + psi[i] + fourpi) % twopi) )

	return (ao, bo)


# finished running sample_sky on the ecliptic plane
subdir = 'ecliptic_plane'

def get_prev_run(subdir):
	ra, dec, counts = [], [], []
	for filename in os.listdir(subdir):
		if '.csv' in filename:
			df = pd.read_csv(subdir+'/'+filename)
			counts.append(df.shape[0])
			ra.append(float(filename.split('_')[0]))
			dec.append(float(filename.split('_')[1]))
	return ra, dec, counts

ra, dec, counts = map(np.array,get_prev_run(subdir))

# use astropysics to convert coords:
from astropysics.coords import EclipticCoordinatesEquinox as Ecliptic
from astropysics.coords import EquatorialCoordinatesEquinox as Equatorial
lam, beta = [], []
for i, j in zip(ra, dec):
	ec = Equatorial(i, j).convert(Ecliptic)
	lam.append(ec.lamb.d)
	beta.append(ec.beta.d)

# or can use euler:
# lam, beta = [], []
# for i, j in zip(ra, dec):
# 	l, b = euler(i, j, select=3)
# 	lam.append(l)
# 	beta.append(b)

# make hexbin plots
plt.figure(figsize=(9,5.5))
plt.hexbin(lam,beta,counts,gridsize=20)
plt.colorbar()
plt.xlabel(r'$\lambda$',fontsize=18)
plt.ylabel(r'$\beta$',fontsize=18)
plt.title('Viable sources: counts per square degree')
# plt.savefig('ecliptic_plane_hexbin_ecliptic.pdf')
plt.show()

plt.figure(figsize=(9,5.5))
plt.hexbin(ra,dec,counts,gridsize=20)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('Viable sources: counts per square degree')
# plt.savefig('ecliptic_plane_hexbin_equatorial.pdf')
plt.show()
