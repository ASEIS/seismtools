List of Tools and Commands:
@userSnapshot.py: original version plotting plane displacements data.

@userSnapshot_v2.py: revised version plotting plane displacements data.
		- sample commands #1: python userSnapshot_v2.py planedisplacements.0 v 0.025 100 136 181 1000 1000 x yes yes multiple 5 no hot yes
		(component = x, unsigned data = yes, plot cumulative data = yes, multiple plots with 5 sec in difference;
		without bar limit; hot colormap; generate data files.)

		- sample commands #2: python userSnapshot_v2.py planedisplacements.0 v 0.025 100 136 181 1000 1000 yz yes single yes 0 1 something.cpt no
		(component = yz, plot cumulative data = yes, single plot; with bar limit 0 and 1; some.cpt colormap; does not generate data files.)



@dataToPlot: loads the 3-column data file and generates the graph.
		- sample commands: python datToPlot.py vel-x-mag-cum-50s.dat 136 181 single hot
		(data file name = "vel-x-mag-cum-50s.dat")

		- sample commands: python datToPlot.py vel-x-mag-cum-10s.dat 136 181 multiple 5 hot
		(data file to start with = "vel-x-mag-cum-10s.dat", search for files after this one)


@extractHer.py: loads the plane displacement data and generate 9-column hercules (station) file at given location.
		- sample commands: python extractHer.py planedisplacements2.0 0.1 50 136 181 1000 1000 60300 90700 test.her
		(x coordinate = 60300, y coordinate = 90700)


@response_spectra_map.py: original version generating the response spectra map; it loads file with user specified sapces.
		- samples commands: python response_spectra_map.py planedisplacements2.0 0.1 50 180000 135000 1000 1000 x d 20 5000 5000 no hot yes test.dat
		(period = 20, space = 5000)


@response_spectra_map2.py: revised version generating the response spectra map; it loads several grid points at one time. The number of grid points should be less than or euqal to number of points in a row. (For example: 181)
		- sample commands: python response_spectra_map2.py planedisplacements2.0 0.1 50 180000 135000 1000 1000 x d 20 181 n hot yes test.dat
		(period = 20, num_grid = 181)


@fourier_spectra_map.py: loads several grid points once and geneates fourier spectra map.
		- sample commands: python fourier_spectra_map.py planedisplacements2.0 0.1 50 180000 135000 1000 1000 x v 2 181 n hot yes test.dat
		(frequence = 2, num_grid = 181)

@userInput.py: holds all the userInput objects involved.

@planeData.py: holds the data objects used in userSnapshot_v2.py

@htools.py: holds all the shared functions.


