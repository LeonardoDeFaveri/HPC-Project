import sys
import csv
import time

class TimeInfo:
  def __init__(self, world_size):
    self.world_size = world_size
    self.time_info = {}

def print_header(info):
  size = len(info.time_info) * 10
  fish_count = 'FISH_COUNT'.center(size)
  print(f'{" ":^10} | {fish_count}')
  print('WORLD_SIZE | ', end='')
  for key in info.time_info:
    print(f'{key:^10}', end='')
  print()
  print('-' * (size + 13))

def print_times(infos):
  print('Program execution times')
  print_header(infos[0])
  for info in infos:
    print(f'{info.world_size:^10} | ', end='')
    for f_count in info.time_info:
      print(f'{info.time_info[f_count]:^10}', end='')
    print()

def print_speedups(infos):
  print('Program speedups')
  info_s = infos[0]
  print_header(info_s)
  for info in infos:
    print(f'{info.world_size:^10} | ', end='')
    for f_count in info.time_info:
      print(f'{info_s.time_info[f_count] / info.time_info[f_count]:^ #10.5f}', end='')
    print()

def print_efficiency(infos):
  print('Program efficiency')
  info_s = infos[0]
  print_header(info_s)
  for info in infos:
    print(f'{info.world_size:^10} | ', end='')
    for f_count in info.time_info:
      eff = info_s.time_info[f_count] / (info.world_size * info.time_info[f_count])
      eff *= 100
      print(f'{eff:^ #10.3f}', end='')
    print()

infos = []

with open(sys.argv[1],'r') as time_results_csv:
  plots = csv.reader(time_results_csv, delimiter=',')
  process_count = 0

  info = None
  for row in plots:
    if (plots.line_num == 1):
      continue
    if (int(row[0]) > process_count):
      process_count = int(row[0])
      info = TimeInfo(int(row[0]))
      infos.append(info)
    info.time_info[int(row[1])] = float(row[2])

  print_times(infos)
  print_speedups(infos)
  print_efficiency(infos)

#plt.ioff()
#figure = plt.figure()
#lines_plotted = plt.plot([[0] for x in cycles], [[0] for y in cycles], "bo", fillstyle="none")
#plt.xlim(-30, 30)
#plt.ylim(-30, 30)
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Fish positions')
#
#def animation(frame):
#  for i, line in enumerate(lines_plotted):
#    line.set_xdata([frame[1].x])
#    line.set_ydata([frame[1].y])
#  plt.legend(lines_plotted,[f'Cycle {frame[0] + 1}'], loc = 'upper right')
#
#anim_created = FuncAnimation(figure, animation, frames=enumerate(cycles),save_count=len(cycles), interval=25)
#writervideo = FFMpegWriter(fps=10) 
#anim_created.save('fish_positions.mp4', writer=writervideo) 
#plt.close()