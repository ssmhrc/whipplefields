import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import urllib
from astropy.io.votable import parse

def get_url(ra,dec,radius,brightlim,faintlim):
	"""
	Returns a URL for a HTTP query of the USNO NOMAD catalog,
	centered on <ra>, <dec> (in decimal degrees), for the given 
	<radius>, bright magnitude limit <brighlitm> and faint 
	magnitude limit <faintlim>. URL will return data in XML/VOTable 
	format.
	"""
	args = [ra, dec, radius, brightlim, faintlim]
	ra, dec, radius, brightlim, faintlim = [str(i) for i in args]
	url='http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=NOMAD'+\
	'&RA='+ra+'&DEC='+dec+'&SR='+radius+'&VERB=3&clr=R2'+\
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
	array1 = np.c_[b[group1], k[group1], theta1, np.ones(theta1.size)]
	array2 = np.c_[b[group2], k[group2], theta2, np.ones(theta2.size)+1]
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

def save_table(table, outpath):
	"""
	Saves the data in <table> in CSV format to the path specified by 
	<outpath>.
	"""
	# df = pd.DataFrame(table.array.tolist(), 
	# 	columns=table.array.dtype.names)
	df = pd.DataFrame.from_records(table.array)
	df.to_csv(outpath, index=False)

# get data for 'murray1' field and calculate angular sizes using the
# relations from the Van Belle 1999 paper
votable = get_table(78.0, 40.0, 1.41, 10.0, 16.0)
table = votable.get_first_table()
save_table(table,'murray1_nomad_test.csv')
theta_vb = get_ang_size_van_belle(table)
plot_groupt_hist(theta_vb,0.3)
plot_groupt_hist(theta_vb,0.1)

# now calculate angular size using the relations from Wang et al. 2010 paper
theta_wang = get_ang_size_wang(table)
plot_groupt_hist(theta_wang,0.3)
plot_groupt_hist(theta_wang,0.1)
