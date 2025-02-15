#!/usr/bin/python3

import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter
import matplotlib.pyplot as plt
import csv
import time

class Cycle:
  def __init__(self, fishes, x, y, weight):
    self.fishes = fishes
    self.x = x
    self.y = y
    self.weight = weight

cycles = []

with open('fish_positions.csv','r') as positions_csv:
  plots = csv.reader(positions_csv, delimiter=',')
  cycle_num = 0
  fishes = []
  x = []
  y = []
  weights = []
  for row in plots:
    if (plots.line_num == 1):
      continue
    if (int(row[0]) > cycle_num): 
      cycles.append(Cycle(fishes, x, y, weights))
      cycle_num = int(row[0])
      fishes = []
      x = []
      y = []
      weights = []
    
    fishes.append(int(row[2]))
    x.append(float(row[3]))
    y.append(float(row[4]))
    weights.append(float(row[5]) / 500.0)
  cycles.append(Cycle(fishes, x, y, weights))

plt.ioff()
figure = plt.figure()
lines_plotted = plt.plot([[0] for x in cycles], [[0] for y in cycles], "bo", fillstyle="none")
plt.xlim(-30, 30)
plt.ylim(-30, 30)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Fish positions')

def animation(frame):
  for i, line in enumerate(lines_plotted):
    line.set_xdata([frame[1].x])
    line.set_ydata([frame[1].y])
  plt.legend(lines_plotted,[f'Cycle {frame[0] + 1}'], loc = 'upper right')

anim_created = FuncAnimation(figure, animation, frames=enumerate(cycles),save_count=len(cycles), interval=25)
writervideo = FFMpegWriter(fps=10) 
anim_created.save('fish_positions.mp4', writer=writervideo) 
plt.close()
