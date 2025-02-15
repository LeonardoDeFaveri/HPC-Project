import sys
import csv
import time

class InfoLine:
  def __init__(self, world_size):
    self.world_size = world_size
    self.info = {}

class Info:
  def __init__(self, time, mean, sd):
    self.time = time
    self.mean = mean
    self.sd = sd

def print_header(info):
  size = len(info.info) * 40
  fish_count = 'FISH_COUNT'.center(size)
  print(f'{" ":^10} | {fish_count}')
  print('WORLD_SIZE | ', end='')
  for key in info.info:
    print(f'{key:^30}', end='')
  print()
  print('-' * (size + 13))

def print_times(infos):
  print('Program execution times')
  print_header(infos[0])
  for info_line in infos:
    print(f'{info_line.world_size:^10} | ', end='')
    for f_count in info_line.info:
      data = info_line.info[f_count]
      print(f'({data.time:^10} : {data.mean:^10} : {data.sd:^10})', end='')
    print()

def print_speedups(infos):
  print('Program speedups')
  info_s = infos[0]
  print_header(info_s)
  for info_line in infos:
    print(f'{info_line.world_size:^10} | ', end='')
    for f_count in info_line.info:
      data = info_line.info[f_count]
      speedup = info_s.info[f_count].time / data.time
      print(f'({speedup:^10f} : {data.mean:^10} : {data.sd:^10})', end='')
    print()

def print_efficiency(infos):
  print('Program efficiency')
  info_s = infos[0]
  print_header(info_s)
  for info_line in infos:
    print(f'{info_line.world_size:^10} | ', end='')
    for f_count in info_line.info:
      data = info_line.info[f_count]
      eff = info_s.info[f_count].time / (info_line.world_size * data.time)
      eff *= 100
      print(f'({eff:^ #10.3f} : {data.mean:^10} : {data.sd:^10})', end='')
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
      info = InfoLine(int(row[0]))
      infos.append(info)
    info.info[int(row[1])] = Info(float(row[2]), float(row[3]), float(row[4]))

  print_times(infos)
  print_speedups(infos)
  print_efficiency(infos)
