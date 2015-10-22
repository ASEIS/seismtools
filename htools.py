#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import sys

def gmtColormap(fileName):
      import colorsys
      """convert cpt file to color dictionary"""
      try:
          f = open(fileName)
      except:
        print "[ERROR]: cpt file not found."
        sys.exit()


      lines = f.readlines()
      f.close()

      x = []
      r = []
      g = []
      b = []
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
             pass
          else:
              x.append(float(ls[0]))
              r.append(float(ls[1]))
              g.append(float(ls[2]))
              b.append(float(ls[3]))
              xtemp = float(ls[4])
              rtemp = float(ls[5])
              gtemp = float(ls[6])
              btemp = float(ls[7])

      x.append(xtemp)
      r.append(rtemp)
      g.append(gtemp)
      b.append(btemp)

      nTable = len(r)
      x = np.array( x , float)
      r = np.array( r , float)
      g = np.array( g , float)
      b = np.array( b , float)
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "RGB":
          r = r/255.
          g = g/255.
          b = b/255.
      xNorm = (x - x[0])/(x[-1] - x[0])

      red = []
      blue = []
      green = []
      for i in range(len(x)):
          red.append([xNorm[i],r[i],r[i]])
          green.append([xNorm[i],g[i],g[i]])
          blue.append([xNorm[i],b[i],b[i]])
      colorDict = {"red":red, "green":green, "blue":blue}
      return colorDict

# cdict = gmtColormap('biasutti.cpt')
# my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
# pcolor(rand(10,10),cmap=my_cmap)
# colorbar()

def integrate(data, dt):
  newdata = cumtrapz(data, dx = dt, initial=0) + data[0]*dt/2.0
  return newdata

def derivative(data, dt):
  """compute derivative of an numpy."""
  newdata = np.insert(data, 0, 0)
  newdata = np.diff(newdata)/dt
  return newdata

# def plot(data, cmap):
#   im = plt.imshow(data, cmap = cmap)
#   plt.axis('off')
#   plt.gca().invert_yaxis()

#   plt.colorbar(im)
#   plt.xlabel('X')
#   plt.ylabel('Y')
#   # plt.suptitle('t = ' + (str)(index), fontsize=20)
#   plt.axis('scaled')
#   plt.show()
# # end of plot

def show_progress(i, num):
  # showing progress on terminal
  percent = float(i)/num
  hashes = '#'*int(round(percent*20))
  spaces = ' '*(20-len(hashes))
  sys.stdout.write("\rProgress: [{0}] {1}%".format(hashes+spaces, int(round(percent*100))))
  sys.stdout.flush()
# end of show_progress

def dis_to_acc(dis, deltaT):
  vel = derivative(dis, deltaT)
  acc = derivative(vel, deltaT)
  return acc
# end of dis_to_acc

def saveDat(filename, dimensionX, spaceX, spaceY, plotData):
  """print the response data in a separate file"""
  try:
    f = open(filename, 'w')
  except IOError, e:
    print e

  descriptor = '{:>12}'*2 + '{:>12.7f}' + '\n'
  x = np.arange(0, dimensionX+1, spaceX, dtype=np.int)
  for i in range(0, len(plotData)):
    y = np.empty(len(plotData[i]), dtype = np.int)
    y.fill(i*spaceY)
    values = plotData[i]
    for c0, c1, c2 in zip(x, y, values):
      f.write(descriptor.format(c0, c1, c2))
  f.close()
# end of saveDat

def plot(data, userInput):
  try:
    if userInput.barChoice == True:
      im = plt.imshow(data, vmin=userInput.barMin,
        vmax=userInput.barMax, cmap=userInput.colorMap)
    else:
      im = plt.imshow(data, cmap=userInput.colorMap)
  except AttributeError:
    im = plt.imshow(data, cmap=userInput.colorMap)

  plt.axis('off')
  plt.gca().invert_yaxis()

  plt.colorbar(im)
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.axis('scaled')
  plt.show()
# end of plot

def loadFile(fp, alongStrike, downDip, num_layers, x_coor, y_coor, size):
  """load the displacement data at given grid point"""
  dis = np.array([], float)
  base = alongStrike*downDip*3*num_layers # number of data in past layers
  index = base + (downDip*x_coor+y_coor)*3 # number of data in past layers + postion in current layer
  offset = index*8 # offset measured in byte

  try:
    dis = np.memmap(fp, np.float64, 'r', offset, (3*size)) # load three numbers for three orientations
    return dis
  except ValueError:
    print "[ERROR]: unable to load file."
    sys.exit()
# end of loadFile


