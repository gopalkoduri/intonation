# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ##Intonation python module
# 
# Homepage: [https://github.com/gopalkoduri/intonation](https://github.com/gopalkoduri/intonation)
# 
# The intonation python module has broadly four classes:
# 
# * **Intervals**: defines an object which has a set of intervals, of course! It has basic functions that facilitate the easy access to these intervals. 
# * **Histogram**: defines a histogram object with methods used to find peaks using different methods and plot them. 
# * ** Pitch**: Given timestamps in seconds, and pitch values in cents, it defines a pitch object which has a number of methods which can be used to study the intervals
# * **Recording**: Given a pitch object, it defines a recording object which has methods to computer histogram, intonation profile and label sections of pitch contours with given set of intervals.

# <codecell>

%pylab inline

import intonation
print dir(intonation)

# <markdowncell>

# ##Load some data
# A sample file with pitch data, and another file with just-intonation intervals are included. 
# The pitch in the file given is in cents scale, normalized to tonic. If you don't have this, you should get it
# from [https://github.com/gopalkoduri/intonation](https://github.com/gopalkoduri/intonation), or 
# load your own data. 
# 
# Make sure the data is formatted as a numpy array of mx2 size where m is number of total points. 
# The first column should corresponds to time stamps in seconds and the second column should 
# correspond to the pitch value in cents. The given file is already formatted this way!

# <codecell>

import pickle
data = loadtxt("88d8196a-123a-4306-9856-4ce3faca14fc.txt")
intervals = pickle.load(file("ji-intervals.pickle"))

# <markdowncell>

# ##Have a look at the data, it always does good!

# <codecell>

plot(data[52000:53000, 0], data[52000:53000, 1])
ylim(200, 1000)

# <markdowncell>

# ##Load the data into a pitch object
# You can avail a number of different method on pitch object to study different aspects of intervals. Let's also look at what methods are 
# available and what they do.

# <codecell>

pitch_obj = intonation.Pitch(data[:, 0], data[:, 1])
help(pitch_obj)

# <markdowncell>

# ##Load the recording object
# Recording object takes the pitch object, and defines methods that access pitch data and functions defined over it, to create
# histogram and intonation profile of the corresponding recording. Load it and check the  methods available on it.

# <codecell>

rec_obj = intonation.Recording(pitch_obj)
help(rec_obj)

# <markdowncell>

# ##Compute intonation profile

# <codecell>

rec_obj.compute_hist(weight='duration')
rec_obj.histogram.get_peaks()
rec_obj.histogram.plot()
rec_obj.parametrize_peaks(intervals)

for peak_pos, parameters in rec_obj.intonation_profile.items():
    print "Peak position:", peak_pos
    print "Parameters:", parameters
    print "\n\n"

