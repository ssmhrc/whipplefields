import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import urllib
from astropy.io.votable import parse, writeto
from scipy.spatial import cKDTree as KDT
from itertools import product


def get_url(ra, dec, radius, brightlim, faintlim):

	"""
	Returns a URL for a HTTP query of the USNO NOMAD catalog,
	centered on <ra>, <dec> (in decimal degrees), for the given 
	<radius>, bright magnitude limit <brighlitm> and faint 
	magnitude limit <faintlim>. URL will return data in XML/VOTable 
	format. See: 
	http://www.usno.navy.mil/USNO/astrometry/optical-IR-prod/icas/vo_nofs
	"""

	args = [ra, dec, radius, brightlim, faintlim]
	ra, dec, radius, brightlim, faintlim = [str(i) for i in args]
	url='http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=NOMAD'+\
	'&RA='+ra+'&DEC='+dec+'&SR='+radius+'&VERB=3&clr=V'+\
	'&bri='+brightlim+'&fai='+faintlim+\
	'&cftype=XML/VO&slf=ddd.dd/dd.ddd&skey=RA'
	return url


def get_votable(ra, dec, radius, brightlim, faintlim):

	"""
	Returns an astropy Table object containing the result of a
	query of the USNO NOMAD catalog, centered on <ra>, <dec> 
	(in decimal degrees), for the given <radius>, 
	bright magnitude limit <brightlim>, and faint magnitude limit
	<faintlim>.
	"""

	url = get_url(ra, dec, radius, brightlim, faintlim)
	xmlfile = urllib.urlretrieve(url)
	xmlpath = xmlfile[0]
	votable = parse(xmlpath)
	return votable


def get_ang_size_van_belle(table):

	"""
	Computes and returns a pandas DataFrame object containing 
	B and K band magnitudes for all sources in table with valid
	entries for both, and the resulting angular size in arcsec
	using the relations from the 1999 Van Belle paper:
	http://adsabs.harvard.edu/abs/1999PASP..111.1515V
	"""

	good = (table.array['B'] != 30.0) & (table.array['K'] != 30.0)
	b = table.array['B'][good]
	k = table.array['K'][good]
	group1 = (-0.6 < b-k) & (b-k < 2.0)		# main sequence
	group2 = (3.0 < b-k) & (b-k < 7.5)		# giants
	group3 = (9.0 < b-k) & (b-k < 16.0)		# variables
	theta1 = 10**(0.500 + 0.290 * (b-k)[group1] - 0.2 * b[group1])
	theta2 = 10**(0.648 + 0.220 * (b-k)[group2] - 0.2 * b[group2])
	theta3 = 10**(0.840 + 0.211 * (b-k)[group3] - 0.2 * b[group3])
	array1 = np.c_[b[group1], k[group1], theta1, np.ones(theta1.size)]
	array2 = np.c_[b[group2], k[group2], theta2, np.ones(theta2.size)+1]
	array3 = np.c_[b[group3], k[group3], theta3, np.ones(theta3.size)+2]
	df = pd.DataFrame(np.r_[array1, array2, array3], 
		columns=['B','K','theta','group'])
	return df


def get_ang_size_wang(table):

	"""
	Computes and returns a pandas DataFrame object containing 
	V and K band magnitudes for all sources in table with valid
	entries for both, and the resulting angular size in arcsec
	using the relations from the 2010 Wang et al. paper:
	http://adsabs.harvard.edu/abs/2010AJ....139.2003W
	"""

	good = (table.array['V'] != 30.0) & (table.array['K'] != 30.0)
	v = table.array['V'][good]
	k = table.array['K'][good]
	group1 = v-k <= 1.85	# main sequence
	group2 = 2.0 < v-k		# giants
	theta1 = 10**(0.453 + 0.246 * (v-k)[group1] - 0.2 * v[group1])
	theta2 = 10**(0.407 + 0.238 * (v-k)[group2] - 0.2 * v[group2])
	array1 = np.c_[v[group1], k[group1], theta1, np.ones(theta1.size)]
	array2 = np.c_[v[group2], k[group2], theta2, np.ones(theta2.size)+1]
	df = pd.DataFrame(np.r_[array1, array2], 
		columns=['V','K','theta','group'])
	return df


def plot_groupt_hist(df, theta_max):

	"""
	Plots a histogram of the angular sizes in the DataFrame objects
	returned by the get_size_* functions, with different colors
	for each group and a maximum angular size of <theta_max>
	"""

	for group in set(df.group):
		idx = (df.group == group) & (df.theta < theta_max)
		binwidth = theta_max/20.
		bins = np.arange(0,theta_max+binwidth,binwidth)
		plt.hist(df.theta[idx], bins=bins, alpha=0.5, normed=True,
			label='group '+str(group))
	plt.legend()
	plt.xlabel('angular size [arcsec]')
	plt.show()


def spherical_to_cartesian(ra, dec):

	"""
	Inputs in degrees.  Outputs x,y,z
	"""

	rar = np.radians(ra)
	decr = np.radians(dec)
	x = np.cos(rar) * np.cos(decr)
	y = np.sin(rar) * np.cos(decr)
	z = np.sin(decr)
	return x, y, z


def radec_to_coords(ra, dec):

	"""
	Helper function for constructing/querying k-d trees with coordinates
	in spherical geometry.
	Converts the input arrays from spherical coordinates to cartesion
	and populates a 3-dimensional array with the result.
	"""

	x, y, z = spherical_to_cartesian(ra, dec)
	coords = np.empty((x.size, 3))
	coords[:, 0] = x
	coords[:, 1] = y
	coords[:, 2] = z
	return coords


def great_circle_distance(ra1, dec1, ra2, dec2):

	"""
	Returns great circle distance.  Inputs in degrees.

	Uses vicenty distance formula - a bit slower than others, but
	numerically stable.
	"""	

	# terminology from the Vicenty formula - lambda and phi and
	# "standpoint" and "forepoint"
	lambs = np.radians(ra1)
	phis = np.radians(dec1)
	lambf = np.radians(ra2)
	phif = np.radians(dec2)

	dlamb = lambf - lambs

	numera = np.cos(phif) * np.sin(dlamb)
	numerb = np.cos(phis) * np.sin(phif) - np.sin(phis) * \
		np.cos(phif) * np.cos(dlamb)
	numer = np.hypot(numera, numerb)
	denom = np.sin(phis) * np.sin(phif) + np.cos(phis) * \
		np.cos(phif) * np.cos(dlamb)
	return np.degrees(np.arctan2(numer, denom))


def get_ang_size(df):

	"""
	Expects a DataFrame object containing V and K band magnitudes
	and calculates the resulting angular size in arcsec
	using the relations from the 2010 Wang et al. paper:
	http://adsabs.harvard.edu/abs/2010AJ....139.2003W
	"""

	group1 = df.V - df.K <= 1.85	# main sequence
	group2 = 2.0 < df.V - df.K		# giants
	theta1 = 10**(0.453 + 0.246 * (df.V - df.K)[group1] - 0.2 * df.V[group1])
	theta2 = 10**(0.407 + 0.238 * (df.V - df.K)[group2] - 0.2 * df.V[group2])
	theta = pd.Series(np.zeros(df.shape[0]))
	group = pd.Series(np.zeros(df.shape[0]))
	theta[group1] = theta1
	theta[group2] = theta2
	group[group1] = 1
	group[group2] = 2
	df['theta'] = theta
	df['group'] = group
	# throw away sources with 1.85 < V-K <= 2
	df = df[group1 | group2]
	return df


def log(*args):
	outdir = args[0]
	f = open(outdir+'/log.txt','a')
	if len(args) == 4:
		line = 'WARNING: no viable sources at position: '+\
			'{}, {}\n\t{}'.format(*args[1:])
	elif len(args) == 5:
		line = 'position ({}, {}): {} out of {} '.format(*args[1:])+\
			'are viable sources'
	f.write(line+'\n')


def cull_dataset(outdir, field_ra, field_dec, table):

	"""
	Efficiently finds all neighbors within 0.01 degrees using
	kdt.query_ball_point method to get points within radius d, where
	d is the cartesian distance equivalent to 0.01 degree separation
	resulting from calculation:

	ra1, ra2 = 0, 0.01
	dec1, dec2 = 0, 0
	c1 = spherical_to_cartesian(ra1, dec1)
	c2 = spherical_to_cartesian(ra2, dec2)
	d = np.sqrt(sum( [ (c1[i] - c2[i])**2 for i in range(3) ] ))

	If there are any neighbors within 0.01 degrees of a given source,
	and if any of these neighbors are brighter than 2 magnitudes fainter
	than the source, remove the source from the table. Also use the
	Wang et al. 2012 relations to get the angular size of each source and
	remove any sources with angular size greater than 0.01 arcsec.

	Returns a pandas DataFrame object containing the culled dataset.
	"""

	good = (table.array['V'] != 30.0) & (table.array['K'] != 30.0)
	arr = table.array[good]
	df = get_ang_size(pd.DataFrame.from_records(arr))

	# ignore sources with theta >= 0.01 arcsec
	df = df[df.theta < 0.01]
	if df.shape[0] == 0:
		return None
	ra, dec, Vmag = [np.array(i) for i in [df.RA, df.DEC, df.V]]
	kdt = KDT(radec_to_coords(ra, dec))
	d = 0.00017453292497790891
	no_neighbors = np.ones(df.shape[0])

	for i in range(df.shape[0]):

		coords = radec_to_coords(ra[i],dec[i])

		# skip the first returned index - this is the query point itself
		idx = kdt.query_ball_point(coords,d)[0][1:]

		if len(idx) < 1:
			continue

		ds = great_circle_distance(ra[i],dec[i],ra[idx],dec[idx])[0]
		Vmag_i = Vmag[i]
		Vmag_neighbors = Vmag[idx]

		# flag sources that have bright nearby neighbors as bad
		for Vmag_j in Vmag_neighbors:
			if Vmag_j - Vmag_i < 2:
				no_neighbors[i] = 0

	df = df[no_neighbors.astype('bool')]
	log(outdir, field_ra, field_dec, df.shape[0], arr.shape[0])
	return df


def sample_sky(ra, dec):

	"""
	Returns a list of the counts of viable sources per square degree
	for the input lists of RA and Dec.
	"""

	assert len(ra) == len(dec)
	deg_sq_radius = np.sqrt(1/np.pi)
	counts = []

	# if previous runs exist, increment outdir number by one
	outdir = 'sample_sky1'
	previous = filter(lambda x: outdir[:-1] in x, os.listdir('.'))
	if previous:
		previous.sort()
		spl = previous[-1].split(outdir[:-1])
		outdir = outdir[:-1] + str(int(spl[1])+1)
	os.mkdir(outdir)

	for i in range(len(ra)):
		try:
			vot = get_votable(ra[i], dec[i], deg_sq_radius, 10.0, 15.0)
		except:
			log(outdir, ra[i], dec[i], 'error on call to get_votable')
			continue
		df = cull_dataset(outdir, ra[i], dec[i], vot.get_first_table())
		if df is None:
			log(outdir, ra[i], dec[i], 'cull_dataset returned None')
			continue
		name = str(ra[i])+'_'+str(dec[i])
		vot_path = outdir+'/'+name+'_VOTable.xml'
		writeto(vot,vot_path)
		csv_path = outdir+'/'+name+'_culled.csv'
		df.to_csv(csv_path, index=False)
		counts.append(df.shape[0])
	return counts
